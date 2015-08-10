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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::SeqFileFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Basename qw(fileparse);

sub write_output {
  my ($self) = @_;
  
  my $species      = $self->param_required('species');
  my $seq_file     = $self->param('seq_file');
  my $species_file = $self->param('species_file');
  my $merge_level  = lc($self->param_required('merge_level'));
  
  if (exists $$species_file{$species}) {
    push @$seq_file, $$species_file{$species};
  }
  
  foreach my $fasta_file (@$seq_file) {
    my $merge_id = 'all';
    if ($merge_level eq 'file') {
      $merge_id = fileparse($fasta_file);
    }
    my $dataflow_output = {
      'fasta_file' => $fasta_file,
      'merge_id'   => $merge_id,
    };
    
    $self->dataflow_output_id($dataflow_output, 1);
  }
}

1;
