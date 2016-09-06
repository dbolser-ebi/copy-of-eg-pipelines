package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref;

use strict;
use warnings;

use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;
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
  $self->delete_protein_features($dba,$logic_name);
}

sub add_xref {
  my ($self, $dbea, $feature, $analysis, $external_db) = @_;
 # my $external_db = $self->param_required('external_db');
  my $hit_name = $feature->hseqname();
  my $ignore_release = 1;
  my $ensembl_start = $feature->start();
  my $ensembl_end = $feature->end();
  my $xref_start = $feature->hstart();
  my $xref_end = $feature->hend();
  my $xref_identity = $feature->percent_id();
  my $evalue = $feature->p_value();
  my $score = $feature->score();
  my $coverage = $feature->coverage();
  my $ensembl_identity = ($xref_end - $xref_start + 1)  * ($xref_identity) / ($ensembl_end - $ensembl_start +      1); 

  my $xref = Bio::EnsEMBL::IdentityXref->new(
    -XREF_IDENTITY    => $xref_identity,
    -ENSEMBL_IDENTITY => $ensembl_identity,
    -EVALUE           => $evalue,
    -QUERY_START       => $xref_start,
    -QUERY_END         => $xref_end,
    -ENSEMBL_START    => $ensembl_start,
    -ENSEMBL_END      => $ensembl_end,
    -ADAPTOR          => $dbea,
    -PRIMARY_ID       => $hit_name,
    -DBNAME           => $external_db,
    -DISPLAY_ID       => $hit_name,
    -SCORE            => $score
  );

  $xref->analysis($analysis);
  $dbea->store($xref, $feature->seqname(), 'Translation', $ignore_release);

}

#sub add_xref {
#  my ($self, $dbea, $feature, $analysis, $external_db) = @_;
#  my $hit_name = $feature->hseqname();
  
#  my $xref = Bio::EnsEMBL::DBEntry->new
#    (
#      -dbname      => $external_db,
#      -primary_id  => $hit_name,
#      -display_id  => $hit_name,
#    );
#  $xref->analysis($analysis);
#  $dbea->store($xref, $feature->seqname(), 'Translation');
#}

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

sub delete_protein_features {
  my ($self, $dba, $logic_name) = @_;
  my $dbh = $dba->dbc->db_handle();
 
  my $sql = " DELETE protein_feature.* FROM protein_feature INNER JOIN analysis ON
  protein_feature.analysis_id = analysis.analysis_id WHERE analysis.logic_name = ?;"
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name) or $self->throw("Failed to execute ($logic_name): $sql");

}

1;
