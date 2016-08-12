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

Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::FetchUniprot

=head1 DESCRIPTION

Download UniProtKB from ftp site and convert into Fasta format.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::FetchUniprot;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use Net::FTP;
use URI;

sub param_defaults {
  return {
    'ebi_path'     => '/ebi/ftp/pub/databases/uniprot/current_release/knowledgebase/complete',
    'ftp_uri'      => 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete',
    'data_source'  => 'sprot',
    'file_varname' => 'uniprot_fasta_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_dir     = $self->param_required('out_dir');
  my $data_source = $self->param_required('data_source');
  
  if (!-e $out_dir) {
    $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
    make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
  }
  
  my $uniprot_file = "uniprot_$data_source.fasta.gz";
  
  $self->param('uniprot_file', $uniprot_file);
}

sub run {
  my ($self) = @_;
  my $ebi_path     = $self->param_required('ebi_path');
  my $ftp_uri      = $self->param_required('ftp_uri');
  my $out_dir      = $self->param_required('out_dir');
  my $uniprot_file = $self->param_required('uniprot_file');
  
  my $ebi_file   = catdir($ebi_path, $uniprot_file);
  my $local_file = catdir($out_dir,  $uniprot_file);
  
  if (-e $ebi_file) {
    $self->fetch_ebi_file($ebi_file, $local_file);
  } else {
    my $ftp = $self->get_ftp($ftp_uri);
    $self->fetch_ftp_file($ftp, $uniprot_file, $local_file);
  }
  
  $self->param('local_file', $local_file);
}

sub write_output {
  my ($self) = @_;
  
  my $output_id = {
    $self->param('file_varname') => $self->param('local_file'),
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
  
  my $ftp_uri  = URI->new($uri);
  my $ftp_host = $ftp_uri->host;
  my $ftp_path = $ftp_uri->path;
  
  my $ftp = Net::FTP->new($ftp_host) or $self->throw("Failed to reach FTP host '$ftp_host': $@");
  $ftp->login or $self->throw(printf("Anonymous FTP login failed: %s.\n", $ftp->message));
  $ftp->cwd($ftp_path) or $self->throw(printf("Failed to change directory to '$ftp_path': %s.\n", $ftp->message));
  $ftp->binary();
  
  return $ftp;
}

sub fetch_ftp_file {
  my ($self, $ftp, $file, $local_file) = @_;
  
  my $remote_size = $ftp->size($file);
  my $remote_mdtm = $ftp->mdtm($file);
  
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
  }
}

1;
