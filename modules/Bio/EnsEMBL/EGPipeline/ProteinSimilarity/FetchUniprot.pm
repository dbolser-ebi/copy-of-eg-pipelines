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

Bio::EnsEMBL::EGPipeline::ProteinSimilarity::FetchUniprot

=head1 DESCRIPTION

Download UniProtKB from ftp site and convert into Fasta format.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::FetchUniprot;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;
use File::Path qw(make_path);
use Net::FTP;
use URI;

sub param_defaults {
  return {
    'ftp_uri'          => 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions',
    'taxonomic_levels' => [qw(fungi invertebrates plants)],
    'uniprot_sources'  => [qw(sprot trembl)],
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_dir = $self->param_required('out_dir');
  
  if (!-e $out_dir) {
    $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
    make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
  }
  
  my $db_fasta_file =
    "$out_dir/".
    join("_", @{$self->param('uniprot_sources')}).
    "_".
    join("_", @{$self->param('taxonomic_levels')}).
    ".fa";
  $self->param('db_fasta_file', $db_fasta_file);
}

sub write_output {
  my ($self) = @_;
  
  my $output_id = {
    'db_fasta_file' => $self->param('db_fasta_file'),
  };
  $self->dataflow_output_id($output_id, 1);
}

sub run {
  my ($self) = @_;
  
  my $ftp = $self->get_ftp();
  my $files = $self->fetch_uniprot_files($ftp);
  $self->convert_to_fasta($files);
}

sub get_ftp {
  my ($self) = @_;
  
  my $ftp_uri = URI->new($self->param('ftp_uri'));
  my $ftp_host = $ftp_uri->host;
  my $ftp_path = $ftp_uri->path;
  
  my $ftp = Net::FTP->new($ftp_host) or $self->throw("Failed to reach FTP host '$ftp_host': $@");
  $ftp->login or $self->throw(printf("Anonymous FTP login failed: %s.\n", $ftp->message));
  $ftp->cwd($ftp_path) or $self->throw(printf("Failed to change directory to '$ftp_path': %s.\n", $ftp->message));
  $ftp->binary();
  
  return $ftp;
}

sub fetch_uniprot_files {
  my ($self, $ftp) = @_;
  
  my $taxonomic_levels = $self->param('taxonomic_levels');
  my $uniprot_sources = $self->param('uniprot_sources');
  my $out_dir = $self->param('out_dir');
  
  my @files = ();
  foreach my $level (@$taxonomic_levels) {
    foreach my $source (@$uniprot_sources) {
      my $file = "uniprot_$source\_$level.dat.gz";
      my $local_file = "$out_dir/$file";
      my $remote_size = $ftp->size($file);
      my $remote_mdtm = $ftp->mdtm($file);
      
      if (-e $local_file) {
        my $local_size = -s $local_file;
        my $local_mdtm = (stat $local_file)[9];
        
        if ( ($remote_size == $local_size) && ($remote_mdtm == $local_mdtm) ) {
          $self->warning("Using existing file '$local_file' with matching timestamp.");
          push @files, $local_file;
          next;
        }
      }
      $ftp->get($file, $local_file) or $self->throw("Failed to get '$file': $@");
      
      # Set the local timestamp to match the remote one.
      utime $remote_mdtm, $remote_mdtm, $local_file;
      push @files, $local_file;
    }
  }
  
  $self->throw("No files for given taxonomic levels ("
    .join(", ", @$taxonomic_levels).") and UniProt sources ("
    .join(", ", @$uniprot_sources).").") unless scalar(@files) > 0;
  
  return \@files;
}

sub convert_to_fasta {
  my ($self, $files) = @_;
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$self->param('db_fasta_file'),
    -format => 'fasta',
  );
  
  foreach my $file (@$files) {
    my $counter = 0;
    
    # Don't forget these are gzipped files...
    open(F, "gunzip -c $file |");
    my $seq_in = Bio::SeqIO->new(
      -fh     => \*F,
      -format => 'swiss',
    );
  
    while (my $inseq = $seq_in->next_seq) {
      $counter++;
      if ($counter >= 100000) {
        $self->warning("Processed $counter sequences.");
        $counter = 0;
      }
      $seq_out->write_seq($inseq);
    }
  }
}

1;
