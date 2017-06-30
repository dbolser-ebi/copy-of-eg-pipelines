=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::EGPipeline::FileDump::CheckSumChecking;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Copy;
use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'skip_file_match' => [],
  };
}

sub run {
  my ($self) = @_;
  
  my ($checking_dir, $new_files) = $self->create_checking_dir;
  $self->classify_files($checking_dir, $new_files);
  remove_tree($checking_dir);
  
  $self->copy_new;
  $self->link_last_release;
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id(
    {
      'new_files'     => $self->param('new_files'),
      'missing_files' => $self->param('missing_files'),
      'changed_files' => $self->param('changed_files'),
    }
  , 1);
}

sub create_checking_dir {
  my ($self) = @_;
  my $results_dir     = $self->param_required('results_dir');
  my $checksum_dir    = $self->param_required('checksum_dir');
  my $skip_file_match = $self->param_required('skip_file_match');
  
  # Create empty dir for comparing new files and old checksums
  my $checking_dir = catdir($checksum_dir, 'checking');
  remove_tree($checking_dir);
  make_path($checking_dir);
  
  # Keep a record of files which are not part of the checksumming;
  # these will need to be treated as new files.
  my @new_files;
  
  # Link to new files
  opendir(my $dh, $results_dir) || die "Failed to open '$results_dir': $!";
  my @file_names = grep { !/.md5/ && !/^\./ && -f "$results_dir/$_" } readdir($dh);
  closedir $dh;
  
  foreach my $file_name (@file_names) {
    my $do_checksum = 1;
    foreach my $match (@$skip_file_match) {
      if ($file_name =~ /\Q$match\E/) {
        $do_checksum = 0;
      }
    }
    if ($do_checksum) {
      symlink "$results_dir/$file_name", "$checking_dir/$file_name" or die $!;
    } else {
      push @new_files, $file_name;
    }
  }
  
  # Link to old checksums
  my $last_release_dir = catdir($checksum_dir, 'last_release');
  if (-e $last_release_dir) {
    opendir($dh, $last_release_dir) || die "Failed to open '$last_release_dir': $!";
    my @md5_files = grep { /.md5/ && -f "$last_release_dir/$_" } readdir($dh);
    closedir $dh;
    
    foreach my $md5_file (@md5_files) {
      symlink "$last_release_dir/$md5_file", "$checking_dir/$md5_file" or die $!;
    }
  }
  
  return ($checking_dir, \@new_files);
}

sub classify_files {
  my ($self, $checking_dir, $new_files) = @_;
  
  opendir(my $dh, $checking_dir) || die "Failed to open '$checking_dir': $!";
  my @files = grep { !/.md5/ && !/^\./ } readdir($dh);
  closedir $dh;
  
  opendir($dh, $checking_dir) || die "Failed to open '$checking_dir': $!";
  my @md5_files = grep { /.md5/ } readdir($dh);
  closedir $dh;
  
  my %files = map { $_ => 1 } @files;
  my %md5_files = map { $_ => 1 } @md5_files;
  
  my @missing_files = ();
  my @changed_files = ();
  
  foreach my $md5_file (sort @md5_files) {
    (my $file = $md5_file) =~ s/\.md5$//;
    
    if (-e "$checking_dir/$file") {
      chdir $checking_dir;
      my $cmd = "md5sum --check --quiet $md5_file 2>&1";
      if (`$cmd` =~ /FAILED/m) {
        push @changed_files, $file;
      }
    } else {
      push @missing_files, $file;
    }
  }
  
  foreach my $file (sort @files) {
    if (! -e "$checking_dir/$file.md5") {
      push @$new_files, $file;
    }
  }
  
  $self->param('new_files', $new_files);
  $self->param('missing_files', \@missing_files);
  $self->param('changed_files', \@changed_files);
}

sub copy_new {
  my ($self) = @_;
  my $results_dir     = $self->param_required('results_dir');
  my $checksum_dir    = $self->param_required('checksum_dir');
  my $release_date    = $self->param_required('release_date');
  my $skip_file_match = $self->param_required('skip_file_match');
  
  my $source_files = catdir($results_dir, '*.md5');
  my $target_dir = catdir($checksum_dir, $release_date);
  make_path($target_dir);
  
  foreach my $source_file (glob $source_files) {
    my $copy_checksum = 1;
    foreach my $match (@$skip_file_match) {
      if ($source_file =~ /\Q$match\E/) {
        $copy_checksum = 0;
      }
    }
    
    if ($copy_checksum) {
      copy ($source_file, $target_dir) or die $!;
    }
  }
}

sub link_last_release {
  my ($self) = @_;
  my $checksum_dir = $self->param_required('checksum_dir');
  my $release_date = $self->param_required('release_date');
  
  my $release_dir = catdir($checksum_dir, $release_date);
  my $last_release_dir = catdir($checksum_dir, 'last_release');
  
  unlink $last_release_dir or die $! if -e $last_release_dir;
  symlink $release_dir, $last_release_dir or die $!;
}

1;
