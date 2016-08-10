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

Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref

=head1 Author

Naveen Kumar

=cut

package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  return {
    
  };
}

sub run {
 my ($self) = @_;

 my $dba  = $self->core_dba();
 my $pfa = $dba->get_adaptor('ProteinFeature');

 my @ProteinFeaturesArray = @{ $pfa->fetch_all_by_logic_name('aedes_blastp') };

 foreach my $ProtFeature (@ProteinFeaturesArray){
   my $dbe_adaptor = $self -> core_dba() -> get_adaptor('DBEntry');
   my $db_name = "RFAM_GENE";





}
1;

