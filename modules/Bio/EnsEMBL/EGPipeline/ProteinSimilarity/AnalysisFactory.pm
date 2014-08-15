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

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::AnalysisFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {};
}

sub run {
  my ($self) = @_;
  my $analyses = $self->param_required('analyses');
  my $db = $self->param_required('db');
  my $blastp = $self->param_required('blastp');
  my $blastx = $self->param_required('blastx');
  
  my %logic_names;
  if ($self->param_required('blastp')) {
    $logic_names{lc($db.'_blastp')} = 1;
  }
  if ($self->param_required('blastx')) {
    $logic_names{lc($db.'_blastx')} = 1;
  }
  
  my $filtered_analyses = [];
  foreach my $analysis (@{$analyses}) {
    $$analysis{'logic_name'} = lc($$analysis{'logic_name'});
    if (exists $logic_names{$$analysis{'logic_name'}}) {
      push @$filtered_analyses, $analysis;
    }
  }
  $self->param('filtered_analyses', $filtered_analyses);
  
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('filtered_analyses'), 1);
  
}

1;

