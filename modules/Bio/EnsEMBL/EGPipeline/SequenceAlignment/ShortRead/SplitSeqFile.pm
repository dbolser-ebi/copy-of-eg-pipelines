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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SplitSeqFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit');

sub write_output {
  my ($self) = @_;
  
  foreach my $split_file (@{$self->param('split_files')}) {
    my $dataflow_output = {
      'seq_file_1' => $split_file,
      'sam_file'   => "$split_file.sam",
    };
    $self->dataflow_output_id($dataflow_output, 2);
  }
}

1;
