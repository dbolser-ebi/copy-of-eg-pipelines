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

Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::MLSSJobFactory;

=head1 DESCRIPTION

=head1 AUTHOR

ckong

=cut
package Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::MLSSJobFactory;

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::SqlHelper;
use base ('Bio::EnsEMBL::Hive::Process');
use Bio::EnsEMBL::Utils::Exception qw(throw);

sub param_defaults {
    return {

           };
}

sub fetch_input {
    my ($self) = @_;

return 0;
}

sub run {
    my ($self)  = @_;

return 0;
}

sub write_output {
    my ($self)  = @_;

    my $compara = $self->param('compara');
    my $from_sp = $self->param('source');
    my $ml_type = $self->param('method_link_type');

    my $mlssa     = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'MethodLinkSpeciesSet');
    my $mlss_list = $mlssa->fetch_all_by_method_link_type($ml_type);

    foreach my $mlss (@$mlss_list){ 
       my $mlss_id = $mlss->dbID();
       my $gdbs    = $mlss->species_set_obj->genome_dbs();

       my @gdb_nm;    
   
       foreach my $gdb (@$gdbs){
          push @gdb_nm,$gdb->name();
       }
 
       # Dataflow only MLSS_ID containing the source species 
       $self->dataflow_output_id({'mlss_id' => $mlss_id, 'compara' => $compara ,'from_sp' => $from_sp }, 2) if(grep (/$from_sp/, @gdb_nm));         
   }
return 0;
}

1;
