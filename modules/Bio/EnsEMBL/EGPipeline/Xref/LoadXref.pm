
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

Bio::EnsEMBL::EGPipeline::Xref::LoadXref

=head1 DESCRIPTION

Generic class for adding xrefs to a core database.

=head1 Author

James Allen

=cut

use strict;
use warnings;

package Bio::EnsEMBL::EGPipeline::Xref::LoadXref;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup');
use Time::Piece;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'oracle_home' => '/sw/arch/dbtools/oracle/product/11.1.0.6.2/client',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  if (!exists $ENV{'ORACLE_HOME'}) {
    $ENV{'ORACLE_HOME'} = $self->param_required('oracle_home');
  }
  
  my $t = localtime;
  $self->param('timestamp', $t->datetime);
}

sub analysis_setup {
  my ($self, $dba) = @_;
  
  my $logic_name = $self->param_required('logic_name');
  my $aa         = $dba->get_adaptor('Analysis');
  my $analysis   = $aa->fetch_by_logic_name($logic_name);
  
  if ($self->param('production_lookup')) {
    $self->production_updates;
  }
  
  my $new_analysis = $self->create_analysis;
  $new_analysis->created($self->param('timestamp'));
  if (defined $analysis) {
    $new_analysis->adaptor($aa);
    $new_analysis->dbID($analysis->dbID);
    $aa->update($new_analysis) || $self->throw("Failed to update analysis '$logic_name'");
  } else {
    $aa->store($new_analysis) || $self->throw("Failed to store analysis '$logic_name'");
  }
}

sub external_db_reset {
  my ($self, $dba, $db_name) = @_;
  
  my $dbh = $dba->dbc->db_handle();
  my $sql = "UPDATE external_db SET db_release = NULL WHERE db_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($db_name) or $self->throw("Failed to execute ($db_name): $sql");
}

sub external_db_update {
  my ($self, $dba, $db_name) = @_;
  
  my $db_release = "EG Xref pipeline; ".$self->param('timestamp');
  my $dbh = $dba->dbc->db_handle();
  my $sql = "UPDATE external_db SET db_release = ? WHERE db_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($db_release, $db_name) or $self->throw("Failed to execute ($db_release, $db_name): $sql");
}

1;
