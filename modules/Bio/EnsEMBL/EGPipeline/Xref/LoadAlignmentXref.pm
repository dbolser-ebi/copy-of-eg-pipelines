package Bio::EnsEMBL::EGPipeline::Xref::LoadAlignmentXref;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;
use Bio::SeqIO;
use List::Util qw(min max);

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'database_type' => 'pep',
    'query_type'    => 'pep',
  };
}

sub run {
  my ($self) = @_;
  my $logic_name  = $self->param_required('logic_name');
  my $external_db = $self->param_required('external_db');

  my $dba = $self->core_dba();
  my $aa  = $dba->get_adaptor('Analysis');

  my $analysis = $aa->fetch_by_logic_name($logic_name);
  
  my ($obj_name, $table_name) = $self->feature_type();

  $self->external_db_reset($dba, $external_db);

  $self->add_xrefs($dba, $analysis, $external_db, $logic_name, $obj_name);

  $self->external_db_update($dba, $external_db);

  $self->delete_features($dba, $logic_name, $table_name);
}

sub feature_type {
  my ($self) = @_;
  my $database_type = $self->param_required('database_type');
  my $query_type    = $self->param_required('query_type');
  
  my ($obj_name, $table_name);
  
  if ($database_type eq 'pep') {
    if ($query_type eq 'pep') {
      $obj_name   = 'ProteinFeature';
      $table_name = 'protein_feature';
    } else {
      $obj_name   = 'ProteinAlignFeature';
      $table_name = 'protein_align_feature';
    }
  } else {
    $obj_name   = 'DnaAlignFeature';
    $table_name = 'dna_align_feature';
  }
  
  return ($obj_name, $table_name);
}

sub add_xrefs {
  my ($self, $dba, $analysis, $external_db, $logic_name, $obj_name) = @_;
  
  my $object_type = $obj_name eq 'ProteinFeature' ? 'Translation' : 'Transcript';
  
  my $xref_metadata = $self->xref_metadata();

  my $dbea = $dba->get_adaptor('DBEntry');

  my $features = $self->fetch_features($dba, $logic_name, $obj_name);

  foreach my $id (keys %$features) {
    foreach my $hit_name (keys %{$$features{$id}}) {
      my $xref;

      my @features = @{$$features{$id}{$hit_name}};
      if (scalar(@features) == 1) {
        $xref = $self->add_xref($analysis, $external_db, $hit_name, $features[0], $xref_metadata);
      } else {
        $xref = $self->add_compound_xref($analysis, $external_db, $hit_name, \@features, $xref_metadata);
      }

      $dbea->store($xref, $id, $object_type);
    }
  }
}

sub fetch_features {
  my ($self, $dba, $logic_name, $obj_name) = @_;
  
  my $fa = $dba->get_adaptor($obj_name);
  my $ta = $dba->get_adaptor('Transcript');

  my %features;
  
  my $all_features = $fa->fetch_all_by_logic_name($logic_name);
  foreach my $feature (@$all_features) {
    my @ids;
    
    if ($obj_name eq 'ProteinFeature') {
      push @ids, $feature->seqname; # Change this to translation_id, pending Ensembl fix.
    } else {
      my $nearby_transcripts = $ta->fetch_all_nearest_by_Feature(
        -FEATURE => $feature,
        -SAME_STRAND => 1,
        #-RANGE => 0,
        -LIMIT => 1,
      );
      
      foreach my $nearby_transcript (@$nearby_transcripts) {
        push @ids, $$nearby_transcript[0]->dbID;
      }
    }
    
    my $hit_name = $feature->hseqname;
    
    foreach my $id (@ids) {
      push @{$features{$id}{$hit_name}}, $feature;
    }
  }
  
  return \%features;
}

sub add_xref {
  my ($self, $analysis, $external_db, $hit_name, $feature, $xref_metadata) = @_;

  my $ensembl_length   = $feature->end - $feature->start + 1;
  my $xref_length      = $feature->hend - $feature->hstart + 1;
  my $xref_identity    = $feature->percent_id;
  my $ensembl_identity = ($xref_length * $xref_identity) / $ensembl_length;

  my $xref = Bio::EnsEMBL::IdentityXref->new(
    -DBNAME           => $external_db,
    -PRIMARY_ID       => $hit_name,
    -DISPLAY_ID       => $$xref_metadata{$hit_name}{display_id},
    -DESCRIPTION      => $$xref_metadata{$hit_name}{description},
    -VERSION          => $$xref_metadata{$hit_name}{version},
    -ENSEMBL_START    => $feature->start,
    -ENSEMBL_END      => $feature->end,
    -QUERY_START      => $feature->hstart,
    -QUERY_END        => $feature->hend,
    -ENSEMBL_IDENTITY => $ensembl_identity,
    -XREF_IDENTITY    => $xref_identity,
    -EVALUE           => $feature->p_value,
    -SCORE            => $feature->score,
		-INFO_TYPE        => 'SEQUENCE_MATCH',
  );
  $xref->analysis($analysis);

  return $xref;
}

sub add_compound_xref {
  my ($self, $analysis, $external_db, $hit_name, $features, $xref_metadata) = @_;

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
    -DISPLAY_ID       => $$xref_metadata{$hit_name}{display_id},
    -DESCRIPTION      => $$xref_metadata{$hit_name}{description},
    -VERSION          => $$xref_metadata{$hit_name}{version},
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

sub xref_metadata {
  my ($self) = @_;
  my $db_fasta_file = $self->param_required('db_fasta_file');

  my %xref_metadata;

  open(F, $db_fasta_file);
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );

  while (my $inseq = $seq_in->next_seq) {
    my ($secondary_id, $desc, $version) = split(/\|/, $inseq->desc);

    $xref_metadata{$inseq->display_id} = {
      display_id  => $secondary_id,
      description => $desc,
      version     => $version,
    };
  }

  return \%xref_metadata;
}

sub delete_features {
  my ($self, $dba, $logic_name, $table_name) = @_;
  my $dbh = $dba->dbc->db_handle();

  my $sql =
    "DELETE $table_name.* FROM $table_name INNER JOIN ".
    "analysis ON $table_name.analysis_id = analysis.analysis_id ".
    "WHERE analysis.logic_name = ?;";
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name) or $self->throw("Failed to execute ($logic_name): $sql");
}

1;
