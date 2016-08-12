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

Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::ExtractSpeciesUniprot

=head1 DESCRIPTION

Extract UniProt entries for a particular species, from Fasta format file.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::ExtractSpeciesUniprot;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;

sub param_defaults {
  return {
    'file_varname' => 'uniprot_species_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $uniprot_fasta_file = $self->param_required('uniprot_fasta_file');
  my $species            = $self->param('species');
  my $uniprot_species    = $self->param('uniprot_species');
  
  if (! defined $uniprot_species) {
    if (! defined $species) {
      $self->throw('-species or -uniprot_species parameter is required');
    } else {
      $uniprot_species = $species;
      $self->param('uniprot_species', $uniprot_species);
    }
  }
  
  (my $output_file = $uniprot_fasta_file) =~ s/(\.\w+)\.gz$/_$uniprot_species$1/;
  
  $self->param('output_file', $output_file);
}

sub run {
  my ($self) = @_;
  my $uniprot_fasta_file = $self->param_required('uniprot_fasta_file');
  my $output_file        = $self->param_required('output_file');
  my $uniprot_species    = $self->param_required('uniprot_species');
  
  $self->extract_species($uniprot_fasta_file, $output_file, $uniprot_species);
}

sub write_output {
  my ($self) = @_;
  
  my $output_id = {
     $self->param('file_varname') => $self->param('output_file'),
  };
  $self->dataflow_output_id($output_id, 1);
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
