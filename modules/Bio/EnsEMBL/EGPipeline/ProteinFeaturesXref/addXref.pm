package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;

use Time::Piece;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $logic_name  = $self->param_required('logic_name');
  my $external_db = $self->param_required('external_db');
  
  my $dba  = $self->core_dba();
  my $pfa  = $dba->get_adaptor('ProteinFeature');
  my $aa   = $dba->get_adaptor('Analysis');
  my $dbea = $dba->get_adaptor('DBEntry');
  
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  my @ProteinFeaturesArray = @{ $pfa->fetch_all_by_logic_name($logic_name) };
  
  $self->external_db_reset($dba, $external_db);
  
  foreach my $feature (@ProteinFeaturesArray){
    $self->add_xref($dbea, $feature, $analysis, $external_db);
  }
  
  $self->external_db_update($dba, $external_db);
}

sub add_xref {
  my ($self, $dbea, $feature, $analysis, $external_db) = @_;
  my $hit_name = $feature->hseqname();
  
  my $xref = Bio::EnsEMBL::DBEntry->new
    (
      -dbname      => $external_db,
      -primary_id  => $hit_name,
      -display_id  => $hit_name,
    );
  
  $xref->analysis($analysis);
  $dbea->store($xref, $feature->seqname(), 'Translation');
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
  
  my $t = localtime;
  
  my $db_release = "EG Alignment Xref pipeline; ".$t->datetime;
  my $dbh = $dba->dbc->db_handle();
  my $sql = "UPDATE external_db SET db_release = ? WHERE db_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($db_release, $db_name) or $self->throw("Failed to execute ($db_release, $db_name): $sql");
}

1;
