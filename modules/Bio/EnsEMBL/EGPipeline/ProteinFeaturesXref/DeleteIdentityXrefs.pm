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

Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::DeleteIdentityXrefs

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::DeleteIdentityXrefs;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $external_db = $self->param_required('external_db');
  
  my $dba  = $self->core_dba();
  my $dbea = $dba->get_adaptor('DBEntry');
  
  my @dbes = @{ $dbea->fetch_all_by_source($external_db) };
  
  foreach my $dbe (@dbes) {
    $self->delete_identity_xref($dba, $dbe);
  }
}

sub delete_identity_xref {
  my ($self, $dba, $dbe) = @_;
  
  my $select_sql = '
    SELECT ox.object_xref_id FROM
      object_xref ox INNER JOIN
      identity_xref ix USING (object_xref_id)
    WHERE ox.xref_id = ?
  ';
  my $select_sth = $dba->dbc->prepare($select_sql);
  $select_sth->execute($dbe->dbID);
  
  my $delete_sql = '
    DELETE ox.*, ix.* FROM
      object_xref ox INNER JOIN
      identity_xref ix USING (object_xref_id)
    WHERE ox.object_xref_id = ?';
  my $delete_sth = $dba->dbc->prepare($delete_sql);
  
  while (my @row = $select_sth->fetchrow_array) {
    my ($ox_id) = @row;
    $delete_sth->execute($ox_id);
  }
}

1;
