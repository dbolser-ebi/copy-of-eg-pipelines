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

Bio::EnsEMBL::EGPipeline::BlastAlignment::ExtractSpecies

=head1 DESCRIPTION

Extract entries for a particular species, from Fasta format file.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::ExtractSpecies;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;

sub param_defaults {
  return {
    'file_varname' => 'species_fasta_file',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $fasta_file     = $self->param_required('fasta_file');
  my $species        = $self->param('species');
  my $source_species = $self->param('source_species');

  if (! defined $source_species) {
    if (! defined $species) {
      $self->throw('-species or -source_species parameter is required');
    } else {
      $source_species = $species;
      $self->param('source_species', $source_species);
    }
  }

  (my $output_file = $fasta_file) =~ s/(\.\w+)$/_$source_species$1/;
  
  $source_species =~ s/_/ /g;
  $source_species =~ s/[A-Z]//g;
  $source_species = ucfirst($source_species);
  $self->param('formatted_species', $source_species);

  $self->param('output_file', $output_file);
}

sub run {
  my ($self) = @_;
  my $fasta_file        = $self->param_required('fasta_file');
  my $output_file       = $self->param_required('output_file');
  my $formatted_species = $self->param_required('formatted_species');
  my $data_source       = $self->param_required('data_source');
  my $data_type         = $self->param_required('data_type');
  
  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$output_file,
    -format => 'fasta',
  );

  open(F, $fasta_file);
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );

  while (my $inseq = $seq_in->next_seq) {
    if ($data_source eq 'uniprot') {
      $self->parse_uniprot($formatted_species, $seq_out, $inseq);
    } elsif ($data_source eq 'refseq') {
      $self->parse_refseq($formatted_species, $data_type, $seq_out, $inseq);
    }
  }
}

sub write_output {
  my ($self) = @_;
  my $output_file  = $self->param('output_file');
  my $file_varname = $self->param_required('file_varname');
  
  if (-s $output_file) {
    $self->dataflow_output_id({$file_varname => $output_file}, 1);
  } else {
    $self->warning("No data for $output_file.");
    $self->input_job->autoflow(0);
  }
}

sub parse_uniprot {
  my ($self, $species, $seq_out, $inseq) = @_;
  my $full_desc = $inseq->desc;
  
  # Some 'descriptions' have "OS=<species>" embedded within them,
  # presumably because UniProt records were used as labels, without
  # appropriate filtering.
  my @OSs = $full_desc =~ /OS=(\w+)/g;
  if (scalar(@OSs) == 2) { 
    $full_desc =~ s/\s+OS=.*?SV=\d+//;
  } elsif (scalar(@OSs) > 2) { 
    $self->throw("More than two 'OS=' sections in description: $full_desc");
  }
  
  if ($full_desc =~ /OS=$species/) {
    my ($desc, $version) = $full_desc =~ /^(.*)\s+OS=.*SV=(\d+)/;
    $inseq->desc(join('|', map { $_ || '' } ($primary_id, $desc, $version)));

    $seq_out->write_seq($inseq);
  }
}

sub parse_refseq {
  my ($self, $species, $data_type, $seq_out, $inseq) = @_;
  
  if ($self->species_match($inseq->desc, $species, $data_type)) {
    my ($primary_id, $version) = $inseq->display_id =~ /(\w+)\.(\d+)|[^\|]+$/;
    $inseq->display_id("$primary_id.$version");
    
    my ($desc) = $inseq->desc =~ /^(.*)\s+(?:\[|\().*/;
    $desc =~ s/PREDICTED:\s+//;
    $inseq->desc(join('|', map { $_ || '' } ("$primary_id.$version", $desc, $version)));

    $seq_out->write_seq($inseq);
  }
}

sub species_match {
  my ($self, $desc, $species, $data_type) = @_;

  if ($data_type eq 'pep') {
    return $desc =~ /\[$species\]/;
  } else {
    return $desc =~ /^$species/;
  }
}

1;
