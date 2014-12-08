
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

=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO

=head1 DESCRIPTION

Runnable that invokes LoadUniProtGO on a core database

=head1 Author

Dan Staines

=cut

use strict;
use warnings;

package Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO;
use base qw/Bio::EnsEMBL::Production::Pipeline::Base/;

use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::EGPipeline::Xref::UniProtGOLoader;

Log::Log4perl->easy_init($INFO);

sub run {
  my ($self)     = @_;
  my $dba        = $self->get_DBAdaptor;
  my $dbc        = $dba->dbc;
  my $species_id = $dba->species_id;
  
  my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
			  -USER   => $self->param('uniparc_user'),
			  -PASS   => $self->param('uniparc_pass'),
			  -HOST   => $self->param('uniparc_host'),
			  -PORT   => $self->param('uniparc_port'),
			  -DBNAME => $self->param('uniparc_dbname'),
			  -DRIVER => $self->param('uniparc_driver')
  );
  my $loader =
	Bio::EnsEMBL::EGPipeline::Xref::UniParcLoader->new(
										   -UNIPARC_DBA=>$uniparc_dba);

  my $logger = get_logger();
  $logger->info("Connecting to core database " . $dba->dbc()->dbname());
  $logger->info("Processing " . $dba->species());
  $loader->add_upis($dba);

} ## end sub run

1;
