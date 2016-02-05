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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SRASeqFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::ProcessSRA');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

sub run {
  my ($self) = @_;
  
  my $work_dir    = $self->param_required('work_directory');
  my $study_id    = $self->param_required('study_id');
  my $merge_level = $self->param_required('merge_level');
  my $taxid       = $self->param('taxid');
  
  my @runs = ();
  
  my $study_adaptor = get_adaptor('Study');
  foreach my $study (@{$study_adaptor->get_by_accession($study_id)}) {
    RUN: foreach my $run (@{$study->runs()}) {

      # Skip this run if it is not the right taxid
      if (defined $taxid) {
          my $sample = $run->sample();
          my $run_taxid = $sample->taxon()->taxon_id();
          next RUN if not $run_taxid eq $taxid;
      }

      my $merge_id = $self->merge_id($merge_level, $study, $run);
      my ($seq_file_1, $seq_file_2, $sam_file) = $self->retrieve_files($work_dir, $run);
      
      push @runs,
        {
          seq_file_1  => $seq_file_1,
          seq_file_2  => $seq_file_2,
          sam_file    => $sam_file,
          merge_id    => $merge_id,
        };
    }
  }
  
  $self->param('runs', \@runs);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $dataflow_output (@{$self->param('runs')}) {
    $self->dataflow_output_id($dataflow_output, 2);
  }
}

1;
