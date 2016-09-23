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

Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast

=head1 DESCRIPTION

get the db_file from the database rather than the config

=head1 Author

Michael Paulini

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast');

sub fetch_runnable {
  my $self = shift @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  if ($self->param_required('blast_type') eq 'wu') {
    $ENV{BLASTMAT} = $self->param_required('blast_matrix');
  }
  $parameters{'TYPE'}     = $self->param_required('blast_type');
  $parameters{'DATABASE'} = $self->param('analysis')->db_file;
  $parameters{'PARSER'}   = $self->make_parser();
  $parameters{'FILTER'}   = undef;
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastEG->new
  (
    -query    => $self->param('query'),
    -program  => $self->param('program'),
    -analysis => $self->param('analysis'),
    -datadir  => $self->param('datadir'),
    -bindir   => $self->param('bindir'),
    -libdir   => $self->param('libdir'),
    -workdir  => $self->param('workdir'),
    %parameters,
  );
  
  if ($self->param_required('database_type') eq 'pep') {
    if ($self->param_required('query_type') eq 'pep') {
      $self->param('results_index', 'translation');
      $self->param('save_object_type', 'ProteinFeature');
    } else {
      $self->param('save_object_type', 'ProteinAlignFeature');
    }
  } else {
    $self->param('save_object_type', 'DnaAlignFeature');
  }
  
  return $runnable;
}

1;
