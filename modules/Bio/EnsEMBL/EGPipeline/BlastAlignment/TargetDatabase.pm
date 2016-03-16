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

Bio::EnsEMBL::EGPipeline::BlastAlignment::TargetDatabase

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::TargetDatabase;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub write_output {
  my ($self) = @_;
  my $blastp = $self->param_required('blastp');
  my $blastx = $self->param_required('blastx');
  
  if ($blastp) {
    $self->dataflow_output_id({}, 2);
  }
  if ($blastx) {
    $self->dataflow_output_id({}, 3);
  }
}

1;
