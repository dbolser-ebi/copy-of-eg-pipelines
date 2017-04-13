=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::FileDump::VEPSpeciesFactory;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory');

sub write_output {
  my ($self) = @_;
  my $core_dbas = $self->param('core_dbas');
  
  foreach my $species (sort keys %$core_dbas) {
    my $dba = $$core_dbas{$species};
    my $mca = $dba->get_MetaContainerAdaptor;
    
    my %species_properties;
    $species_properties{host}            = $dba->dbc->host;
    $species_properties{port}            = $dba->dbc->port;
    $species_properties{user}            = 'ensro';
    $species_properties{pass}            = '';
    $species_properties{group}           = 'core';
    $species_properties{dbname}          = $dba->dbc->dbname;
    $species_properties{is_multispecies} = $mca->is_multispecies;
    $species_properties{species}         = $species;
    $species_properties{species_id}      = $dba->species_id;
    $species_properties{assembly}        = $mca->single_value_by_key('assembly.default');
    $species_properties{regulation}      = undef;
    $species_properties{variation}       = undef;
    
    my $dbva = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'variation');
    if (defined $dbva) {
      $species_properties{variation} = $dbva->dbc->dbname;
    }
    
    $self->dataflow_output_id(\%species_properties, $self->param('core_flow'));
  }
}

1;
