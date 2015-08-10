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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::FindRuns;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

sub run {
  my ($self) = @_;
  
  my $species       = $self->param_required('species');
  my $study_list    = $self->param('study');
  my $species_study = $self->param('species_study');
  my $merge_level   = lc($self->param_required('merge_level'));
  
  if (exists $$species_study{$species}) {
    push @$study_list, $$species_study{$species};
  }
  
  my @runs = ();
  foreach my $study_acc (@$study_list) {
    my $study_adaptor = get_adaptor('Study');
    foreach my $study (@{$study_adaptor->get_by_accession($study_acc)}) {
      foreach my $run (@{$study->runs()}) {
        my $experiment = $run->experiment;
        my $sample = $run->sample;
        
        my $merge_id;
        if ($merge_level eq 'study') {
          $merge_id = $study->accession;
        } elsif ($merge_level eq 'taxon') {
          $merge_id = $sample->taxon ? $sample->taxon->taxon_id : 'unknown';
        } elsif ($merge_level eq 'sample') {
          $merge_id = $sample->accession;
        } elsif ($merge_level eq 'experiment') {
          $merge_id = $experiment->accession;
        } else {
          $merge_id = 'all';
        }
        
        push @runs,
          {
            # Also create and pass output template, of the format <study>_<merge_id>_<ass.default>
            run_acc  => $run->accession,
            merge_id => $merge_id,
          };
      }
    }
  }
  
  $self->param('runs', \@runs);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $run (@{$self->param('runs')}) {
    $self->dataflow_output_id($run, 2);
  }
}

1;
