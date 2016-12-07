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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::AnalysisFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::BlastAlignment::AnalysisFactory');

sub run {
  my ($self) = @_;
  my $analyses          = $self->param_required('analyses');
  my $logic_name_prefix = $self->param_required('logic_name_prefix');
  my $blastp            = $self->param_required('blastp');
  my $blastx            = $self->param_required('blastx');
  my $filter_top_x      = $self->param('filter_top_x');
  my $source_species    = $self->param('source_species');
  
  my $proteomic_analyses = [];
  my $genomic_analyses = [];
  
  foreach my $analysis (@{$analyses}) {
    my $db_file = $$analysis{db_file};
    $db_file =~ s!.*/([^/]+)$!$1!;
    $$analysis{'description'} .= " Source database: $db_file.";
    $$analysis{blast_db} = $$analysis{db_file};   
 
    if ($filter_top_x) {
      my $top_x;
      if ($$analysis{'program'} =~ /blastp$/) {
        $top_x = $self->param_required('blastp_top_x');
      } elsif ($$analysis{'program'} =~ /blastx$/) {
        $top_x = $self->param_required('blastx_top_x');
      }
      if ($top_x eq 1) {
        $$analysis{'description'} .= " Results are filtered to show the best unique hit, based on E-value.";
      } else {
        $$analysis{'description'} .= " Results are filtered to show the top $top_x unique hits, based on E-value.";
      }
    }
    
    if ($$analysis{'program'} =~ /blastp/ && $blastp) {
      push @$proteomic_analyses, $analysis;
    } elsif ($$analysis{'program'} =~ /blastx/ && $blastx) {
      push @$genomic_analyses, $analysis;
    }
  }
 
  $self->param('proteomic_analyses', $proteomic_analyses);
  $self->param('genomic_analyses', $genomic_analyses);
}

sub write_output {
  my ($self) = @_;
  my $blastp = $self->param_required('blastp');
  my $blastx = $self->param_required('blastx');
  
  if ($blastp) {
    $self->dataflow_output_id($self->param('proteomic_analyses'), 2);
  }
  
  if ($blastx) {  
    $self->dataflow_output_id($self->param('genomic_analyses'), 3);
  }
}

1;
