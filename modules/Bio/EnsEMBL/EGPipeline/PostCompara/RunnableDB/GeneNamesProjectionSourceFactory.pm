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

=head1 NAME

Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::GeneNamesProjectionSourceFactory;

=head1 DESCRIPTION

=head1 AUTHOR

ckong

=cut
package Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::GeneNamesProjectionSourceFactory;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::Base');

sub run{
    my ($self)  = @_;

    my $self = shift @_;
    my $species_list = $self->param('species_list');
    my $g_config = $self->param_required('g_config');
    my $final_species_list;

  if ($species_list)
  {
    $final_species_list=$species_list;
  }
  else
  {
    $final_species_list=$g_config;
  }
    foreach my $pair (sort (keys $final_species_list)){
       my $source                 = $final_species_list->{$pair}->{'source'};
       my $species                = $final_species_list->{$pair}->{'species'};
       my $antispecies            = $final_species_list->{$pair}->{'antispecies'};
       my $division               = $final_species_list->{$pair}->{'division'};
       my $run_all                = $final_species_list->{$pair}->{'run_all'};       
       my $method_link_type       = $final_species_list->{$pair}->{'method_link_type'};  
       my $homology_types_allowed = $final_species_list->{$pair}->{'homology_types_allowed'};
       my $percent_id_filter      = $final_species_list->{$pair}->{'percent_id_filter'};
       my $percent_cov_filter     = $final_species_list->{$pair}->{'percent_cov_filter'};
       my $taxon_filter           = $final_species_list->{$pair}->{'taxon_filter'};
       my $geneName_source        = $final_species_list->{$pair}->{'geneName_source'};
       my $geneDesc_rules         = $final_species_list->{$pair}->{'geneDesc_rules'};
       my $geneDesc_rules_target  = $final_species_list->{$pair}->{'geneDesc_rules_target'};

       # Remove source/target species from the hash
      delete $final_species_list->{$pair};
     $self->param('species_list', $final_species_list);
       $self->dataflow_output_id(
		{'source'      		  => $source, 
		 'species'     		  => $species, 
		 'antispecies' 		  => $antispecies, 
  		 'division'    	  	  => $division, 
		 'run_all' 		  => $run_all,
		 'method_link_type' 	  => $method_link_type,
                 'homology_types_allowed' => $homology_types_allowed,
  		 'percent_id_filter'      => $percent_id_filter,
		 'percent_cov_filter'     => $percent_cov_filter,
		 'taxon_filter'           => $taxon_filter,
		 'geneName_source'	  => $geneName_source,
		 'geneDesc_rules'	  => $geneDesc_rules,
		 'geneDesc_rules_target'  => $geneDesc_rules_target
		},2);
       $self->dataflow_output_id({'species_list'       => $self->param('species_list'),
                                 'species'                => $species,
                                 'source'                 => $source},1);
       last;
      }
return 0;
}

1;


