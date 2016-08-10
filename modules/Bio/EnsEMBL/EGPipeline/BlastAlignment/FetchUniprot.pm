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

Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot

=head1 DESCRIPTION

Download UniProtKB from ftp site and convert into Fasta format.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot;

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
    'ebi_path'        => '/ebi/ftp/pub/databases/uniprot/current_release/knowledgebase',
    'ftp_uri'         => 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase',
    'taxonomic_level' => 'species',
    'data_source'     => 'sprot',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_dir         = $self->param_required('out_dir');
  my $taxonomic_level = $self->param_required('taxonomic_level');
  my $data_source     = $self->param_required('data_source');
  my $species         = $self->param('species');
  my $uniprot_species = $self->param('uniprot_species');
  
  if (!-e $out_dir) {
    $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
    make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
  }
  
  if (! defined $uniprot_species) {
    if (! defined $species) {
      if ($taxonomic_level eq 'species') {
        $self->throw('-species or -uniprot_species parameter is required');
      }
    } else {
      $uniprot_species = $species;
      $self->param('uniprot_species', $uniprot_species);
    }
  }
  
  my ($uniprot_file, $output_file);
  
  if ($taxonomic_level eq 'species') {
    $self->param('sub_dir', 'complete');
    $uniprot_file = "uniprot_$data_source.fasta.gz";
    $output_file = catdir($out_dir, "$data_source\_$uniprot_species.fa");
  } else {
    $self->param('sub_dir', 'taxonomic_divisions');
    $uniprot_file = "uniprot_$data_source\_$taxonomic_level.dat.gz";
    $output_file = catdir($out_dir, "$data_source\_$taxonomic_level.fa");
  }
  
  $self->param('uniprot_file', $uniprot_file);
  $self->param('output_file', $output_file);
}

sub run {
  my ($self) = @_;
  my $ebi_path        = $self->param_required('ebi_path');
  my $ftp_uri         = $self->param_required('ftp_uri');
  my $sub_dir         = $self->param_required('sub_dir');
  my $out_dir         = $self->param_required('out_dir');
  my $uniprot_file    = $self->param_required('uniprot_file');
  my $output_file     = $self->param_required('output_file');
  my $taxonomic_level = $self->param_required('taxonomic_level');
  my $uniprot_species = $self->param('uniprot_species');
  
  my $ebi_file   = catdir($ebi_path, $sub_dir, $uniprot_file);
  my $local_file = catdir($out_dir, $uniprot_file);
  
  if (-e $ebi_file) {
    $self->fetch_ebi_file($ebi_file, $local_file);
  } else {
    my $uri = "$ftp_uri/$sub_dir";
    my $ftp = $self->get_ftp($uri);
    $self->fetch_ftp_file($ftp, $uniprot_file, $local_file);
  }
  
  if ($taxonomic_level eq 'species') {
    $self->extract_species($local_file, $output_file, $uniprot_species);
  } else {
    $self->convert_to_fasta($local_file, $output_file);
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

sub convert_to_fasta {
  my ($self, $file, $output_file) = @_;
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$output_file,
    -format => 'fasta',
  );
  
  my ($total, $counter) = (0, 0);
  
  # Don't forget these are gzipped files...
  open(F, "gunzip -c $file |");
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'swiss',
  );
  
  while (my $inseq = $seq_in->next_seq) {
    $total++;
    $counter++;
    if ($counter >= 100000) {
      $self->warning("Processed $total sequences.");
      $counter = 0;
    }
    $seq_out->write_seq($inseq);
  }
}

sub extract_species {
  my ($self, $file, $output_file, $species) = @_;
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$output_file,
    -format => 'fasta',
  );
  
  $species =~ s/_/ /g;
  $species =~ s/[A-Z]//g;
  $species = ucfirst($species);
  
  my ($total, $counter) = (0, 0);
  
  # Don't forget these are gzipped files...
  open(F, "gunzip -c $file |");
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );
  
  while (my $inseq = $seq_in->next_seq) {
    $total++;
    $counter++;
    if ($counter >= 100000) {
      $self->warning("Processed $total sequences.");
      $counter = 0;
    }
    
    if ($inseq->desc =~ /OS=$species/) {
      my $display_id = $inseq->display_id;
      $display_id =~ s/^[^\|]+\|([^\|]+).*/$1/;
      $inseq->display_id($display_id);
      $seq_out->write_seq($inseq);
    }
  }
}

1;
