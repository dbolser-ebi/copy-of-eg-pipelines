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

sub param_defaults {
  return {
    'file_varname' => 'species_fasta_file',
    'data_type'    => 'pep',
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

  $self->param('output_file', $output_file);
}

sub run {
  my ($self) = @_;
  my $fasta_file     = $self->param_required('fasta_file');
  my $output_file    = $self->param_required('output_file');
  my $source_species = $self->param_required('source_species');
  my $data_type      = $self->param_required('data_type');

  $self->extract_species($fasta_file, $output_file, $source_species, $data_type);
}

sub write_output {
  my ($self) = @_;
  my $output_file  = $self->param('output_file');
  my $file_varname = $self->param_required('file_varname')
  my $data_type    = $self->param_required('data_type');
  
  if (-s $output_file) {
    my $output_ids = {
      $file_varname => $output_file,
      data_type     => $data_type,
    };
    $self->dataflow_output_id($output_ids, 1);
  } else {
    my $source_species = $self->param('source_species');
    $self->warning("No $data_type data for $source_species.");
    $self->input_job->autoflow(0);
  }
}

1;
