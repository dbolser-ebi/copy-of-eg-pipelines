=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::Xref::ExtractSpeciesFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub write_output {
  my ($self) = @_;
  my $logic_name = $self->param_required('logic_name');

  if ($logic_name =~ /sprot|trembl/) {
    $self->dataflow_output_id({}, 2);
  } elsif ($logic_name =~ /refseq/) {
    $self->dataflow_output_id({}, 3);
  }
}

1;
