=head1 LICENSE

Copyright [1999-2014] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchRefSeq

=head1 DESCRIPTION

Download RefSeq files from ftp site and convert into Fasta format.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchRefSeq;

use strict;
use warnings;
use feature 'say';
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;
use File::Basename qw(fileparse);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use Net::FTP;
use URI;

sub param_defaults {
  return {
    'ebi_path'  => '/nfs/panda/ensemblgenomes/external/refseq',
    'ftp_uri'   => 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release',
    'data_type' => 'protein',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_dir            = $self->param_required('out_dir');
  my $taxonomic_division = $self->param_required('taxonomic_division');
  my $data_type          = $self->param_required('data_type');
  my $species            = $self->param_required('species');
  my $refseq_species     = $self->param('refseq_species');
  
  if (!-e $out_dir) {
    $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
    make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
  }
  
  if (! defined $refseq_species) {
    $refseq_species = $species;
    $self->param('refseq_species', $refseq_species);
  }
  
  my $refseq_files;
  if ($data_type eq 'protein') {
    $refseq_files = "$taxonomic_division.*.protein.faa.gz";
  } else {
    $refseq_files = "$taxonomic_division.*.rna.fna.gz";
  }
  my $output_file  = catdir($out_dir, "refseq_$refseq_species\_$data_type.fa");
  
  $self->param('refseq_files', $refseq_files);
  $self->param('output_file', $output_file);
}

sub run {
  my ($self) = @_;
  my $refseq_ebi_path    = $self->param_required('refseq_ebi_path');
  my $refseq_ftp_uri     = $self->param_required('refseq_ftp_uri');
  my $taxonomic_division = $self->param_required('taxonomic_division');
  my $data_type          = $self->param_required('data_type');
  my $out_dir            = $self->param_required('out_dir');
  my $refseq_files       = $self->param_required('refseq_files');
  my $output_file        = $self->param_required('output_file');
  my $refseq_species     = $self->param_required('refseq_species');
  
  my $ebi_files = catdir($refseq_ebi_path, $taxonomic_division, $refseq_files);
  
  # The trickiness here is that the refseq data is stored in multiple
  # files. We will iterate over them, appending to the output file;
  # so need to make sure any existing file is removed.
  if (-e $output_file) {
    unlink $output_file;
  }
  
  my @ebi_files = glob($ebi_files);
  my @local_files;
  
  if (scalar(@ebi_files)) {
    foreach my $ebi_file (@ebi_files) {
      my $filename = fileparse($ebi_file);
      my $local_file = catdir($out_dir, $filename);
      
      $self->fetch_ebi_file($ebi_file, $local_file);
      push @local_files, $local_file;
    }
  } else {
    my $uri = "$refseq_ftp_uri/$taxonomic_division";
    my $ftp = $self->get_ftp($uri);
    @local_files = @{ $self->fetch_ftp_file($ftp, $refseq_files, $out_dir) };
  }
  
  foreach my $local_file (@local_files) {
    $self->extract_species($local_file, $output_file, $refseq_species, $data_type);
  }
}

sub write_output {
  my ($self) = @_;
  
  my $output_id = {
    'db_fasta_file' => $self->param('output_file'),
  };
  $self->dataflow_output_id($output_id, 1);
}

sub fetch_ebi_file {
  my ($self, $file, $local_file) = @_;

  my $ebi_size = -s $file;
  my $ebi_mdtm = (stat $file)[9];
  
  if (-e $local_file) {
    my $local_size = -s $local_file;
    my $local_mdtm = (stat $local_file)[9];
    
    if ( ($ebi_size == $local_size) && ($ebi_mdtm == $local_mdtm) ) {
      $self->warning("Using existing file '$local_file' with matching timestamp.");
    } else {
      copy($file, $local_file) or $self->throw("Failed to get '$file': $@");
    }
  } else {
    copy($file, $local_file) or $self->throw("Failed to get '$file': $@");
  }
  
  # Set the local timestamp to match the remote one.
  utime $ebi_mdtm, $ebi_mdtm, $local_file;
  
  if (! -e $local_file) {
    $self->throw("Failed to copy file '$file'.");
  }
}

sub get_ftp {
  my ($self, $uri) = @_;
  
  my $ftp_uri = URI->new($uri);
  my $ftp_host = $ftp_uri->host;
  my $ftp_path = $ftp_uri->path;
  
  my $ftp = Net::FTP->new($ftp_host) or $self->throw("Failed to reach FTP host '$ftp_host': $@");
  $ftp->login or $self->throw(printf("Anonymous FTP login failed: %s.\n", $ftp->message));
  $ftp->cwd($ftp_path) or $self->throw(printf("Failed to change directory to '$ftp_path': %s.\n", $ftp->message));
  $ftp->binary();
  
  return $ftp;
}

sub fetch_ftp_file {
  my ($self, $ftp, $files_pattern, $out_dir) = @_;
  
  my @all_files = $ftp->ls();
  my @local_files;
  
  my @files = grep { $_ =~ /$files_pattern/ } @all_files;
  for my $file (@files) {
    my $remote_size = $ftp->size($file);
    my $remote_mdtm = $ftp->mdtm($file);
    
    my $local_file = catdir($out_dir, $file);
    
    if (-e $local_file) {
      my $local_size = -s $local_file;
      my $local_mdtm = (stat $local_file)[9];
      
      if ( ($remote_size == $local_size) && ($remote_mdtm == $local_mdtm) ) {
        $self->warning("Using existing file '$local_file' with matching timestamp.");
      } else {
        $ftp->get($file, $local_file) or $self->throw("Failed to get '$file': $@");
      }
    } else {
      $ftp->get($file, $local_file) or $self->throw("Failed to get '$file': $@");
    }
    
    # Set the local timestamp to match the remote one.
    utime $remote_mdtm, $remote_mdtm, $local_file;
    
    if (! -e $local_file) {
      $self->throw("Failed to download file '$file'.");
    } else {
      push @local_files, $local_file;
    }
  }
  
  return \@local_files;
}

sub extract_species {
  my ($self, $file, $output_file, $species, $data_type) = @_;
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>>'.$output_file,
    -format => 'fasta',
  );
  
  $species =~ s/_/ /g;
  $species =~ s/[A-Z]//g;
  $species = ucfirst($species);
  
  # Don't forget these are gzipped files...
  open(F, "gunzip -c $file |");
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );
  
  while (my $inseq = $seq_in->next_seq) {
    if ($self->species_match($inseq->desc, $species, $data_type)) {
      my $display_id = $inseq->display_id;
      $display_id =~ s/.*?([^\|]+)\|?$/$1/;
      $inseq->display_id($display_id);
      $seq_out->write_seq($inseq);
    }
  }
}

sub species_match {
  my ($self, $desc, $species, $data_type) = @_;
  
  if ($data_type eq 'protein') {
    return $desc =~ /\[$species\]/;
  } else {
    return $desc =~ /^$species/;
  }
}
1;
