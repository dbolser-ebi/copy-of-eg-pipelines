package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref;

use strict;
use warnings;

use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;
use List::Util qw(min max);
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
  
  my %protein_features;
  
  foreach my $feature (@ProteinFeaturesArray) {
    my $translation_id = $feature->seqname();
    my $hit_name       = $feature->hseqname();
    push @{$protein_features{$translation_id}{$hit_name}}, $feature;
  }
  
  $self->add_xrefs($dbea, $analysis, $external_db, \%protein_features);
    
  $self->external_db_update($dba, $external_db);
  $self->delete_protein_features($dba,$logic_name);
}

sub add_xrefs {
  my ($self, $dbea, $analysis, $external_db, $protein_features) = @_;
  my $ignore_release = 1;
  
  foreach my $translation_id (keys %$protein_features) {
    foreach my $hit_name (keys %{$$protein_features{$translation_id}}) {
      my $xref;
      
      my @features = @{$$protein_features{$translation_id}{$hit_name}};
      if (scalar(@features) == 1) {
        $self->add_xref($dbea, $features[0], $analysis, $external_db);
      } else {
        $self->add_compound_xref($dbea, $analysis, $external_db, $hit_name, \@features);
      }
      
      $dbea->store($xref, $translation_id, 'Translation', $ignore_release);
    }
  }
}

sub add_xref {
  my ($self, $dbea, $feature, $analysis, $external_db) = @_;
  
  my $hit_name = $feature->hseqname();
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
    -QUERY_START      => $xref_start,
    -QUERY_END        => $xref_end,
    -ENSEMBL_START    => $ensembl_start,
    -ENSEMBL_END      => $ensembl_end,
    -ADAPTOR          => $dbea,
    -PRIMARY_ID       => $hit_name,
    -DBNAME           => $external_db,
    -DISPLAY_ID       => $hit_name,
    -SCORE            => $score,
  );

  $xref->analysis($analysis);

  return $xref;
}

sub add_compound_xref {
  my ($self, $dbea, $analysis, $external_db, $hit_name, $features) = @_;
  
  my @ensembl_starts;
  my @ensembl_ends;
  my @xref_starts;
  my @xref_ends;
  my ($ensembl_length, $xref_length, $xref_match);
  
  foreach my $feature (@$features) {
    push @ensembl_starts, $feature->start;
    push @ensembl_ends,   $feature->end;
    push @xref_starts,    $feature->hstart;
    push @xref_ends,      $feature->hend;
    
    $ensembl_length += $feature->end  - $feature->start  + 1;
    $xref_length    += $feature->hend - $feature->hstart + 1;
    $xref_match     += ($feature->hend - $feature->hstart + 1) * $feature->percent_id;
  }
  
  my $ensembl_start = min(@ensembl_starts);
  my $ensembl_end = max(@ensembl_ends);
  my $xref_start = min(@xref_starts);
  my $xref_end = max(@xref_ends);
  my $xref_identity = sprintf("%.2f", $xref_match/$xref_length);
  my $ensembl_identity = sprintf("%.2f", $xref_match/$ensembl_length); 
  
  my $xref = Bio::EnsEMBL::IdentityXref->new(
    -XREF_IDENTITY    => $xref_identity,
    -ENSEMBL_IDENTITY => $ensembl_identity,
    -QUERY_START      => $xref_start,
    -QUERY_END        => $xref_end,
    -ENSEMBL_START    => $ensembl_start,
    -ENSEMBL_END      => $ensembl_end,
    -ADAPTOR          => $dbea,
    -PRIMARY_ID       => $hit_name,
    -DBNAME           => $external_db,
    -DISPLAY_ID       => $hit_name,
  );

  $xref->analysis($analysis);

  return $xref;
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

sub delete_protein_features {
  my ($self, $dba, $logic_name) = @_;
  my $dbh = $dba->dbc->db_handle();
 
  my $sql = "DELETE protein_feature.* FROM protein_feature INNER JOIN analysis ON protein_feature.analysis_id = analysis.analysis_id WHERE analysis.logic_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name) or $self->throw("Failed to execute ($logic_name): $sql");

}

1;
