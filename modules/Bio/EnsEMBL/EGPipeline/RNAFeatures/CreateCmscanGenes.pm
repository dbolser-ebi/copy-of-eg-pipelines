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

Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateCmscanGenes

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateCmscanGenes;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::RNAFeatures::CreateGenes');

sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults},
    'evalue_threshold'     => 1e-6,
    'truncated'            => 0,
    'nonsignificant'       => 0,
    'bias_threshold'       => 0.3,
    'gene_overlap_sources' => ['mirbase_gene', 'trnascan_gene'],
  };
}

sub check_thresholds {
  my ($self, $feature) = @_;
  my $evalue_threshold = $self->param_required('evalue_threshold');
  my $truncated        = $self->param_required('truncated');
  my $nonsignificant   = $self->param_required('nonsignificant');
  my $bias_threshold   = $self->param_required('bias_threshold');
  
  my $evalue  = $feature->p_value;
  my $attribs = $feature->get_all_Attributes;
  
  my ($biotype, $cmscan_bias, $cmscan_significant, $cmscan_truncated);
  foreach my $attrib (@$attribs) {
    if ($attrib->code eq 'rna_gene_biotype') {
      $biotype = $attrib->value;
    }
    if ($attrib->code eq 'cmscan_bias') {
      $cmscan_bias = $attrib->value;
    }
    if ($attrib->code eq 'cmscan_significant') {
      $cmscan_significant = 1;
    }
    if ($attrib->code eq 'cmscan_truncated') {
      $cmscan_truncated = 1;
    }
  }
  
  my $pass_thresholds = 0;
  
  if (defined $evalue) {
    if ($evalue <= $evalue_threshold) {
      $pass_thresholds = 1;
    }
    
    if ($biotype eq 'misc_RNA') {
      $pass_thresholds = 0;
    }
    
    if (defined $cmscan_bias) {
      if ($cmscan_bias > $bias_threshold) {
        $pass_thresholds = 0;
      }
    }
    
    if (defined $cmscan_truncated) {
      if (! $truncated && $cmscan_truncated) {
        $pass_thresholds = 0;
      }
    }
    
    if (defined $cmscan_significant) {
      if (! $nonsignificant && ! $cmscan_significant) {
        $pass_thresholds = 0;
      }      
    }
  }
  
  return $pass_thresholds;
}

1;
