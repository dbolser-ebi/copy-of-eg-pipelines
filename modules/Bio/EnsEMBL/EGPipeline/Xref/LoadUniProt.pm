
=head1 LICENSE

Copyright [1999-2015] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt

=head1 DESCRIPTION

Add UniProt xrefs to a core database.

=head1 Author

James Allen

=cut

use strict;
use warnings;

package Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::EGPipeline::Xref::UniProtLoader;

sub run {
  my ($self) = @_;
  
  my $uniparc_db  = $self->param_required('uniparc_db');
  my $uniprot_db  = $self->param_required('uniprot_db');
  my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniparc_db);
  my $uniprot_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniprot_db);
  my $loader      = Bio::EnsEMBL::EGPipeline::Xref::UniProtLoader->new
  (
    -UNIPARC_DBA  => $uniparc_dba,
    -UNIPROT_DBA  => $uniprot_dba,
    -REPLACE_ALL  => $self->param('replace_all'),
    -GENE_NAMES   => $self->param('gene_names'),
    -DESCRIPTIONS => $self->param('descriptions'),
  );
  
  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  $self->analysis_setup($dba);
  $self->external_db_reset($dba, 'Uniprot/SPTREMBL');
  $self->external_db_reset($dba, 'Uniprot/SWISSPROT');
  $loader->add_xrefs($dba);
  $self->external_db_update($dba, 'Uniprot/SPTREMBL');
  $self->external_db_update($dba, 'Uniprot/SWISSPROT');
}

1;
