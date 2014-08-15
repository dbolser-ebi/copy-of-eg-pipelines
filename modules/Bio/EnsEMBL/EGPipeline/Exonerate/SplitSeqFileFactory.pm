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

package Bio::EnsEMBL::EGPipeline::Exonerate::SplitSeqFileFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {};
}

sub run {
  my ($self) = @_;
  my $species = $self->param_required('species');
  my $seq_file = $self->param('seq_file');
  my $seq_files = $self->param('seq_files');
  
  if (exists $$seq_files{$species}) {
    $self->param('fasta_file', $$seq_files{$species});
  } else {
    if (defined $seq_file) {
      $self->param('fasta_file', $seq_file);
    } else {
      $self->throw("No seq_file for $species.");
    }
  }
  
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({'fasta_file' => $self->param('fasta_file')}, 1);
  
}

1;

