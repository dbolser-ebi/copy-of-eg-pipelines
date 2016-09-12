package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref;

use strict;
use warnings;

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
  my $protein_features = $pfa->fetch_all_by_logic_name($logic_name);

  $self->external_db_reset($dba, $external_db);

  my %protein_features;

  foreach my $feature (@$protein_features) {
    my $translation_id = $feature->seqname();
    my $hit_name       = $feature->hseqname();
    push @{$protein_features{$translation_id}{$hit_name}}, $feature;
  }

  $self->add_xrefs($dbea, $analysis, $external_db, \%protein_features);

  $self->external_db_update($dba, $external_db);

  #$self->delete_protein_features($dba, $logic_name);
}

sub add_xrefs {
  my ($self, $dbea, $analysis, $external_db, $protein_features) = @_;

  foreach my $translation_id (keys %$protein_features) {
    foreach my $hit_name (keys %{$$protein_features{$translation_id}}) {
      my $xref;

      my @features = @{$$protein_features{$translation_id}{$hit_name}};
      if (scalar(@features) == 1) {
        $xref = $self->add_xref($analysis, $external_db, $hit_name, $features[0]);
      } else {
        $xref = $self->add_compound_xref($analysis, $external_db, $hit_name, \@features);
      }

      $dbea->store($xref, $translation_id, 'Translation');
    }
  }
}

sub add_xref {
  my ($self, $analysis, $external_db, $hit_name, $feature) = @_;

  my $ensembl_length   = $feature->end - $feature->start + 1;
  my $xref_length      = $feature->hend - $feature->hstart + 1;
  my $xref_identity    = $feature->percent_id;
  my $ensembl_identity = ($xref_length * $xref_identity) / $ensembl_length;

  my $xref = Bio::EnsEMBL::IdentityXref->new(
    -DBNAME           => $external_db,
    -PRIMARY_ID       => $hit_name,
    -DISPLAY_ID       => $hit_name,
    -ENSEMBL_START    => $feature->start,
    -ENSEMBL_END      => $feature->end,
    -QUERY_START      => $feature->hstart,
    -QUERY_END        => $feature->hend,
    -ENSEMBL_IDENTITY => $ensembl_identity,
    -XREF_IDENTITY    => $xref_identity,
    -EVALUE           => $feature->p_value,
    -SCORE            => $feature->score,
  );
  $xref->analysis($analysis);

  return $xref;
}

sub add_compound_xref {
  my ($self, $analysis, $external_db, $hit_name, $features) = @_;

  my @ensembl_starts;
  my @ensembl_ends;
  my @xref_starts;
  my @xref_ends;
  my ($ensembl_length, $xref_length, $matches);

  foreach my $feature (@$features) {
    push @ensembl_starts, $feature->start;
    push @ensembl_ends,   $feature->end;
    push @xref_starts,    $feature->hstart;
    push @xref_ends,      $feature->hend;

    $ensembl_length +=  $feature->end  - $feature->start  + 1;
    $xref_length    +=  $feature->hend - $feature->hstart + 1;
    $matches        += ($feature->hend - $feature->hstart + 1) * $feature->percent_id;
  }

  my $xref_identity    = sprintf("%.2f", $matches/$xref_length);
  my $ensembl_identity = sprintf("%.2f", $matches/$ensembl_length);

  my $xref = Bio::EnsEMBL::IdentityXref->new(
    -DBNAME           => $external_db,
    -PRIMARY_ID       => $hit_name,
    -DISPLAY_ID       => $hit_name,
    -ENSEMBL_START    => min(@ensembl_starts),
    -ENSEMBL_END      => max(@ensembl_ends),
    -QUERY_START      => min(@xref_starts),
    -QUERY_END        => max(@xref_ends),
    -ENSEMBL_IDENTITY => $ensembl_identity,
    -XREF_IDENTITY    => $xref_identity,
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

  my $sql =
    "DELETE protein_feature.* FROM protein_feature INNER JOIN ".
    "analysis ON protein_feature.analysis_id = analysis.analysis_id ".
    "WHERE analysis.logic_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name) or $self->throw("Failed to execute ($logic_name): $sql");
}

1;
