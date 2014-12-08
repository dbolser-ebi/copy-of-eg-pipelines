
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

Bio::EnsEMBL::EGPipeline::Xref::UniProtPomBaseLoader

=head1 DESCRIPTION



=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::UniProtPomBaseLoader;
use base Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Digest::MD5;
use List::MoreUtils qw/uniq/;
use Data::Dumper;

=head1 CONSTRUCTOR
=head2 new
  Arg [-UNIPROT_DBA]  : 
       string - adaptor for UniProt Oracle database (e.g. SWPREAD)
  Arg [-REPLACE_ALL]    : 
       boolean - remove all UniProt references first

  Example    : $ldr = Bio::EnsEMBL::EGPipeline::Xref::UniProtGOLoader->new(...);
  Description: Creates a new loader object
  Returntype : Bio::EnsEMBL::EGPipeline::Xref::UniProtGOLoader
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ($proto, @args) = @_;
  my $self = $proto->SUPER::new(@args);
  ($self->{uniprot_dba}, $self->{replace_all})
	= rearrange(
				['UNIPROT_DBA',
				 'REPLACE_ALL'],
				@args);
  return $self;
}

=head1 METHODS
=head2 add_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Add xrefs to supplied core
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut

sub add_xrefs {
  my ($self, $dba) = @_;
  $self->{analysis} = $self->get_analysis($dba, 'xrefpombase');
  if (defined $self->{replace_all} && $self->{replace_all} == 1) {
	$self->remove_xrefs($dba);
  }
  # get translation_id,pombase
  my $translation_pombase = $self->get_translation_pombase($dba);
  $self->add_uniprot_xrefs($dba, $translation_pombase);
  return;
}

=head2 get_translation_pombase
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Get mapping of UniParc to translation
  Returntype : Hash of translation IDs to UniParc accessions
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub get_translation_pombase {
  my ($self, $dba) = @_;
  $self->logger()->info("Finding translation-pombase pairs");

  my $translation_pombase = {};
  $dba->dbc()->sql_helper()->execute_no_return(
	-SQL => q/
		select tl.translation_id,g.stable_id,g.gene_id
	from translation tl
	join transcript tr using (transcript_id)
	join gene g using (gene_id)
	join seq_region sr on (g.seq_region_id=sr.seq_region_id)
	join coord_system cs using (coord_system_id)
        where
        tr.biotype='protein_coding' and
	cs.species_id=?
	/,
	-CALLBACK => sub {
	  my ($tid, $pombase, $gene_id) = @{$_[0]};
	  $translation_pombase->{$tid} = {pombase => $pombase, gene_id => $gene_id};
	  return;
	},
	-PARAMS => [$dba->species_id()]);

  $self->logger()
	->info("Found " .
		   scalar(keys %$translation_pombase) . " translation-pombase pairs");

  return $translation_pombase;
} ## end sub get_translation_pombase

=head2 add_uniprot_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : hashref of translation ID to PomBase accessions
  Description: Add UniProt to specified translations
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub add_uniprot_xrefs {
  my ($self, $dba, $translation_pombase) = @_;
  my $ddba  = $dba->get_DBEntryAdaptor();
  my $gdba  = $dba->get_GeneAdaptor();
  my $tN    = 0;
  my $uN    = 0;
  $self->logger()->info("Adding UniProt xrefs");
  # hash of gene names, descriptions
  my $gene_attribs = {};
  while (my ($tid, $pombase) = each %$translation_pombase) {
	$tN++;
	$self->logger()->debug("Looking up entry for " . $pombase->{pombase});
	my $uniprots = $self->get_uniprot_for_pombase($pombase->{pombase});
	$uN +=
	  $self->store_uniprot_xrefs($ddba, $tid, $uniprots,
								 $pombase->{gene_id}, $gene_attribs);
	$self->logger()->info("Processed $tN translations ($uN xrefs)")
	  if ($tN % 1000 == 0);
  }
  $self->logger()->info("Stored $uN UniProt xrefs on $tN translations");
} ## end sub add_uniprot_xrefs

=head2 store_uniprot_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : Translation dbID
  Arg        : Bio::EnsEMBL::DBEntry for UniProt record
  Arg        : Corresponding gene dbID
  Arg        : Hash of gene-related identifiers
  Description: Add UniProt xrefs to specified translation
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub store_uniprot_xrefs {
  my ($self, $ddba, $tid, $uniprots, $gene_id, $gene_attribs) = @_;
  # remove existing uniprots for this translation first
  $ddba->dbc()->sql_helper()->execute_update(
	-SQL => q/
		delete ox.*,ix.*
		from object_xref ox
		join  xref x using (xref_id)
		join external_db d using (external_db_id)
		left join identity_xref ix using (object_xref_id)
		where
		d.db_name in ('Uniprot\/SWISSPROT','Uniprot\/SPTREMBL')
		and ox.ensembl_id = ?
		and ox.ensembl_object_type = 'Translation'
	/,
	-PARAMS => [$tid]);

  my $n = 0;
  return $n if scalar(@$uniprots) == 0;
  for my $uniprot (@$uniprots) {

	if (!$uniprot->{ac} || $uniprot->{ac} eq '') {
	  $self->logger()
		->warn(
"Empty $uniprot->{type} accession retrieved from UniProt for translation $tid");
	  next;
	}
	$self->logger()
	  ->debug(
		  "Storing $uniprot->{type} " . $uniprot->{ac} . " on translation $tid");
	my $dbentry =
	  $ddba->fetch_by_db_accession($uniprot->{type}, $uniprot->{ac});
	if (!defined $dbentry) {
	  $dbentry =
		Bio::EnsEMBL::DBEntry->new(
								-PRIMARY_ID  => $uniprot->{ac},
								-DISPLAY_ID  => $uniprot->{name},
								-VERSION     => $uniprot->{version},
								-DESCRIPTION => $uniprot->{description},
								-DBNAME      => $uniprot->{type});
	}
	$dbentry->analysis($self->{analysis});
	$ddba->store($dbentry, $tid, 'Translation');
	$n++;
  } ## end for my $uniprot (@$uniprots)
  return $n;
} ## end sub store_uniprot_xrefs


=head2 get_uniprot_for_pombase
  Arg        : Pombase identifier
  Description: Find UniProt accession for PomBase
  Returntype : hashref of matching UniProt records
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub get_uniprot_for_pombase {
  my ($self, $pombase) = @_;
  $self->logger()->debug("Finding UniProt accessions for $pombase");

  my $sql = q/SELECT d.accession from dbentry d 
  join gene g on ( d.dbentry_id = g.dbentry_id) 
  join gene_name gn on (g.gene_id = gn.gene_id) 
  join cv_gene_name_type cgnt on (gn.gene_name_type_id = cgnt.gene_name_type_id ) 
  WHERE gn.name=? and cgnt.type='ORFNames' and d.deleted='N' and d.entry_type in (0,1) and d.merge_status <> 'R'/;

  my @uniprot_acs = @{
	$self->{uniprot_dba}->dbc()->sql_helper()->execute_simple(
	  -SQL => $sql,
	  -PARAMS => [$pombase])
  };
  

  my $uniprots = [];

  for my $ac (@uniprot_acs) {
  $self->logger()->debug("Building xref for $ac");
	my $uniprot = {};
	$self->{uniprot_dba}->dbc()->sql_helper()->execute_no_return(
	  -SQL => q/
	  SELECT d.accession,
  d.name,
  REPLACE(NVL(sc1.text,sc3.text),'^'),
  d.entry_type, gn.name, cgnt.type,
  s.version
FROM SPTR.dbentry d
JOIN sequence s ON (s.dbentry_id=d.dbentry_id)
LEFT OUTER JOIN SPTR.dbentry_2_description dd
ON (dd.dbentry_id         = d.dbentry_id
AND dd.description_type_id=1)
LEFT OUTER JOIN SPTR.description_category dc1
ON (dd.dbentry_2_description_id=dc1.dbentry_2_description_id
AND dc1.category_type_id       =1)
LEFT OUTER JOIN SPTR.description_subcategory sc1
ON (dc1.category_id        = sc1.category_id
AND sc1.subcategory_type_id=1)
LEFT OUTER JOIN SPTR.description_category dc3
ON (dd.dbentry_2_description_id=dc3.dbentry_2_description_id
AND dc3.category_type_id       =3)
LEFT OUTER JOIN SPTR.description_subcategory sc3
ON (dc3.category_id        = sc3.category_id
AND sc3.subcategory_type_id=1)
left join gene g on ( d.dbentry_id = g.dbentry_id)
left join gene_name gn on (g.gene_id = gn.gene_id)
left join cv_gene_name_type cgnt on (gn.gene_name_type_id = cgnt.gene_name_type_id )
WHERE d.accession = ?
	/,
	  -PARAMS   => [$ac],
	  -CALLBACK => sub {
		my ($ac, $name, $des, $type, $gene_name, $gene_name_type, $version) =
		  @{$_[0]};
		if (defined $ac && $ac ne '') {
		  $uniprot->{ac} = $ac;
		}
		if (defined $des && $des ne '') {
		  $uniprot->{description} = $des;
		}
		if ( defined $version && $version ne '' ) {
		  $uniprot->{version} = $version;
		}
		if (defined $name && $name ne '') {
		  if (defined $uniprot->{name}) {
			$uniprot->{name} .= "; $name";
		  }
		  else {
			$uniprot->{name} = $name;
		  }
		}
		if (defined $type && $type ne '') {
		  $uniprot->{type} =
			$type == 0 ? "Uniprot/SWISSPROT" : "Uniprot/SPTREMBL";
		}
		if (defined $gene_name) {
		  if ($gene_name_type eq 'Name') {
			$uniprot->{gene_name} = $gene_name;
		  }
		  else {
			push @{$uniprot->{synonyms}}, $gene_name;
		  }
		}
		return;
	  });
	if (defined $uniprot->{ac}) {
	  push @$uniprots, $uniprot;
	}
  } ## end for my $ac (@uniprot_acs)
  return $uniprots;

} ## end sub get_uniprot_for_pombase

=head2 remove_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Remove existing UniProt cross-references from genome
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub remove_xrefs {
  my ($self, $dba) = @_;
  $self->logger()->info("Removing existing UniProt cross-references");
  $dba->dbc()->sql_helper()->execute_update(
	-SQL => q/
		delete ox.*, ix.*
		from object_xref ox
		join xref x using(xref_id)
		join external_db d using (external_db_id)
		join translation tl on (tl.translation_id=ox.ensembl_id and ox.ensembl_object_type = 'Translation')
		join transcript tr using (transcript_id)
		join seq_region sr using (seq_region_id)
		join coord_system cs using (coord_system_id)
		left join identity_xref ix using (object_xref_id)
		where
		d.db_name in ('Uniprot\/SWISSPROT','Uniprot\/SPTREMBL')
		and cs.species_id=?
	/,
	-PARAMS => [$dba->species_id()]);
  return;
}

1;
