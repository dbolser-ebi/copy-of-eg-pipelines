
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

Bio::EnsEMBL::EGPipeline::Xref::UniProtLoader

=head1 DESCRIPTION

Module for adding UniProt cross-references based on UPI identity. 
Note that where matches are found, existing cross-references for that translation are removed. 
This allows imperfect matches to be replaced with exact matches as they become available
in UniProt

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::UniProtLoader;
use base Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Digest::MD5;
use List::MoreUtils qw/uniq/;
use Data::Dumper;
use Bio::EnsEMBL::IdentityXref;

=head1 CONSTRUCTOR
=head2 new
  Arg [-UNIPROT_DBA]  : 
       string - adaptor for UniProt Oracle database (e.g. SWPREAD)
  Arg [-UNIPARC_DBA]  : 
       string - adaptor for UniParc Oracle database (e.g. UAPRO)
  Arg [-REPLACE_ALL]    : 
       boolean - remove all UniProt references first
  Arg [-GENE_NAMES]    : 
       boolean - add gene names from SwissProt
  Arg [-DESCRIPTIONS]    : 
       boolean - add descriptions from SwissProt

  Example    : $ldr = Bio::EnsEMBL::EGPipeline::Xref::UniProtGOLoader->new(...);
  Description: Creates a new loader object
  Returntype : Bio::EnsEMBL::EGPipeline::Xref::UniProtGOLoader
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ( $proto, @args ) = @_;
  my $self = $proto->SUPER::new(@args);
  ( $self->{uniparc_dba}, $self->{uniprot_dba},
	$self->{replace_all}, $self->{gene_names},
	$self->{descriptions} )
	= rearrange(
				 [ 'UNIPARC_DBA', 'UNIPROT_DBA',
				   'REPLACE_ALL', 'GENE_NAMES',
				   'DESCRIPTIONS' ],
				 @args );
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
  my ( $self, $dba ) = @_;
  $self->{analysis} = $self->get_analysis( $dba, 'xrefuniparc' );
  if ( defined $self->{replace_all} && $self->{replace_all} == 1 ) {
	$self->remove_xrefs($dba);
  }
  $self->add_uniprot_xrefs( $dba,
							$self->get_translation_upis($dba),
							$self->get_translation_lengths($dba) );
  return;
}

=head2 get_translation_upis
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Get mapping of UniParc to translation
  Returntype : Hash of translation IDs to UniParc accessions
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub get_translation_upis {
  my ( $self, $dba ) = @_;
  $self->logger()->info("Finding translation-UPI pairs");

  my $translation_upis = {};
  $dba->dbc()->sql_helper()->execute_no_return(
	-SQL => q/
		select tl.translation_id,upix.dbprimary_acc,tr.gene_id
	from translation tl
	join transcript tr using (transcript_id)
	join seq_region sr using (seq_region_id)
	join coord_system cs using (coord_system_id)
	join object_xref upiox on (upiox.ensembl_object_type='Translation' and upiox.ensembl_id=tl.translation_id)
	join xref upix using (xref_id) 
	join external_db upie using (external_db_id) 
        where
        tr.biotype='protein_coding' and
	cs.species_id=? and 
	upie.db_name='UniParc'
	/,
	-CALLBACK => sub {
	  my ( $tid, $upi, $gene_id ) = @{ $_[0] };
	  $translation_upis->{$tid} = { upi => $upi, gene_id => $gene_id };
	  return;
	},
	-PARAMS => [ $dba->species_id() ] );

  $self->logger()
	->info( "Found " .
		 scalar( keys %$translation_upis ) . " translation-UPI pairs" );

  return $translation_upis;
} ## end sub get_translation_upis

sub get_translation_lengths {
  my ( $self, $dba ) = @_;

  $self->logger()->info("Finding translation lengths");

  my $translation_lengths = {};
  $dba->dbc()->sql_helper()->execute_no_return(
	-SQL => q/
		select tl.translation_id,ta.value
	from translation tl
	join translation_attrib ta using (translation_id)
	join attrib_type a using (attrib_type_id)
	join transcript tr using (transcript_id)
	join seq_region sr using (seq_region_id)
	join coord_system cs using (coord_system_id)
        where
        tr.biotype='protein_coding' and
	cs.species_id=? and 
	a.code='NumResidues'
	/,
	-CALLBACK => sub {
	  my ( $tid, $len ) = @{ $_[0] };
	  $translation_lens->{$tid} = $len;
	  return;
	},
	-PARAMS => [ $dba->species_id() ] );

  $self->logger()
	->info( "Found " .
		   scalar( keys %$translation_lens ) . " translation lengths" );

  return $translation_lens;

} ## end sub get_translation_lengths

=head2 add_uniprot_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : hashref of translation ID to UniParc accessions
  Description: Add UniProt to specified translations
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub add_uniprot_xrefs {
  my ( $self, $dba, $translation_upis, $translation_lengths ) = @_;
  my $taxid = $dba->get_MetaContainer()->get_taxonomy_id();
  my $ddba  = $dba->get_DBEntryAdaptor();
  my $gdba  = $dba->get_GeneAdaptor();
  my $tN    = 0;
  my $uN    = 0;
  $self->logger()->info("Adding UniProt xrefs");
  # hash of gene names, descriptions
  my $gene_attribs = {};
  while ( my ( $tid, $upi ) = each %$translation_upis ) {

	$tN++;
	$self->logger()->debug( "Looking up entry for " . $upi->{upi} );
	my $uniprots = $self->get_uniprot_for_upi( $taxid, $upi->{upi} );
	$uN +=
	  $self->store_uniprot_xrefs( $ddba,
								  $tid,
								  $uniprots,
								  $upi->{gene_id},
								  $gene_attribs,
								  $translation_lengths->{$tid} );
	$self->logger()->info("Processed $tN translations ($uN xrefs)")
	  if ( $tN % 1000 == 0 );
  }

  if ( defined $self->{gene_names} && $self->{gene_names} == 1 ) {
	$self->set_gene_names( $ddba, $gene_attribs );
  }
  if ( defined $self->{descriptions} && $self->{descriptions} == 1 ) {
	$self->set_descriptions( $gdba, $gene_attribs );
  }
  $self->logger()->info("Stored $uN UniProt xrefs on $tN translations");
} ## end sub add_uniprot_xrefs

=head2 store_uniprot_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : Translation dbID
  Arg        : Bio::EnsEMBL::DBEntry for UniProt record
  Arg        : Corresponding gene dbID
  Arg        : Hash of gene-related identifiers
  Arg        : Length of translation
  Description: Add UniProt xrefs to specified translation
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub store_uniprot_xrefs {
  my ( $self, $ddba, $tid, $uniprots, $gene_id, $gene_attribs, $len ) =
	@_;
  my $n = 0;
  return $n if scalar(@$uniprots) == 0;

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
	-PARAMS => [$tid] );

  for my $uniprot (@$uniprots) {

	if ( !$uniprot->{ac} || $uniprot->{ac} eq '' ) {
	  $self->logger()
		->warn(
"Empty $uniprot->{type} accession retrieved from UniProt for translation $tid"
		);
	  next;
	}
	$self->logger()
	  ->debug( "Storing $uniprot->{type} " .
			   $uniprot->{ac} . " on translation $tid" );

	my $dbentry =
	  Bio::EnsEMBL::IdentityXref->new(
								-DBNAME      => $uniprot->{type},
								-PRIMARY_ID  => $uniprot->{ac},
								-DISPLAY_ID  => $uniprot->{name},
								-DESCRIPTION => $uniprot->{description},
								-VERSION     => $uniprot->{version},
								-XREF_IDENTITY    => 100,
								-ENSEMBL_IDENTITY => 100,
								-EVALUE           => 0,
								-SCORE            => 100,
								-CIGAR_LINE       => $len . 'M',
								-XREF_START       => 1,
								-XREF_END         => $len,
								-ENSEMBL_START    => 1,
								-ENSEMBL_END      => $len,
								-ANALYSIS         => $self->{analysis}
	  );

	$ddba->store( $dbentry, $tid, 'Translation' );
	# track names and descriptions
	if ( defined $uniprot->{description} &&
		 $uniprot->{type} eq 'Uniprot/SWISSPROT' )
	{
	  push @{ $gene_attribs->{descriptions}->{$gene_id}
		  ->{ $uniprot->{type} }->{ $uniprot->{description} } },
		'[Source:' . $uniprot->{type} . ';Acc:' . $uniprot->{ac} . ']';
	}
	if ( defined $uniprot->{gene_name} &&
		 $uniprot->{type} eq 'Uniprot/SWISSPROT' )
	{
	  $gene_attribs->{gene_names}->{$gene_id}->{ $uniprot->{type} }
		->{ $uniprot->{gene_name} } += 1;
	  if ( defined $uniprot->{synonyms} ) {
		for my $synonym ( @{ $uniprot->{synonyms} } ) {
		  push(
			  @{ $gene_attribs->{synonyms}->{ $uniprot->{gene_name} } },
			  $synonym );
		}
	  }
	}
	$n++;
  } ## end for my $uniprot (@$uniprots)
  return $n;
} ## end sub store_uniprot_xrefs

=head2 set_descriptions
  Arg        : Bio::EnsEMBL::DBSQL::GeneAdaptor for core database to write to
  Arg        : Hash of gene-related identifiers
  Description: Add descriptions to names based on matched UniProt records
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub set_descriptions {
  my ( $self, $gdba, $gene_attribs ) = @_;
  my $nDes = 0;
  while ( my ( $gid, $descrs ) =
		  each %{ $gene_attribs->{descriptions} } )
  {
	# work out the best description

	my $candidates = $descrs->{'Uniprot/SWISSPROT'};
	if ( !defined $candidates ) {
	  $candidates = $descrs->{'Uniprot/SPTREMBL'};
	}
	my @descs = sort {
	  scalar @{ $candidates->{$b} } <=> scalar @{ $candidates->{$a} }
	} keys %{$candidates};
	if ( scalar @descs > 0 ) {
	  my $description = $descs[0];
	  if ( defined $description ) {
		my $src = $candidates->{$description};
		if ( defined $src && scalar @{$src} > 0 ) {
		  $description .= " " . $src->[0];
		}
		# store description for gene
		$self->logger()
		  ->debug(
				 "Setting gene $gene_id description to '$description'");
		$nDes++;
		$gdba->dbc()->sql_helper()->execute_update(
			   -SQL => q/update gene set description=? where gene_id=?/,
			   -PARAMS => [ $description, $gid ] );
	  }
	}
  } ## end while ( my ( $gid, $descrs...))
  $self->logger()->info("Stored $nDes descriptions");
  return;
} ## end sub set_descriptions

=head2 set_gene_names
  Arg        : Bio::EnsEMBL::DBSQL::DBEntryAdaptor for core database to write to
  Arg        : Hash of gene-related identifiers
  Description: Add gene names to names based on matched UniProt records
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub set_gene_names {
  my ( $self, $ddba, $gene_attribs ) = @_;
  my $nNames = 0;
  while ( my ( $gid, $names ) = each %{ $gene_attribs->{gene_names} } )
  {
	# work out the best name
	my $candidates = $names->{'Uniprot/SWISSPROT'};
	if ( !defined $candidates ) {
	  $candidates = $names->{'Uniprot/SPTREMBL'};
	}
	my @gene_names =
	  sort { $candidates->{$b} <=> $candidates->{$a} }
	  keys %{$candidates};
	if ( scalar @gene_names > 0 ) {
	  my $gene_name = $gene_names[0];
	  if ( defined $gene_name ) {
		# create dbentry
		my $gd =
		  $ddba->fetch_by_db_accession( 'Uniprot_gn', $gene_name );
		if ( !defined $gd ) {
		  $gd =
			Bio::EnsEMBL::DBEntry->new( -PRIMARY_ID => $gene_name,
										-DISPLAY_ID => $gene_name,
										-DBNAME     => 'Uniprot_gn' );
		  # synonyms for this name
		  my $synonyms = $gene_attribs->{synonyms}->{$gene_name};
		  if ( defined $synonyms ) {
			for my $synonym ( uniq @$synonyms ) {
			  $gd->add_synonym($synonym);
			}
		  }
		  $ddba->store($gd);
		  if ( !defined $gd ) {
			$self->logger()
			  ->warn(
				   "Could not store xref for gene name " . $gene_name );
			next;
		  }
		}
		# set as display_xref_id
		$self->logger()
		  ->debug("Setting gene $gene_id name to '$gene_name'");
		$nNames++;
		$ddba->dbc()->sql_helper()->execute_update(
		   -SQL => q/update gene set display_xref_id=? where gene_id=?/,
		   -PARAMS => [ $gd->dbID(), $gid ] );
	  } ## end if ( defined $gene_name)
	} ## end if ( scalar @gene_names...)
  } ## end while ( my ( $gid, $names...))
  $self->logger()->info("Stored $nNames gene names");
  return;
} ## end sub set_gene_names

=head2 get_uniprot_for_upi
  Arg        : Taxonomy ID
  Arg        : UniParc identifier
  Description: Find UniProt accession for UPI and taxonomy
  Returntype : hashref of matching UniProt records
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub get_uniprot_for_upi {
  my ( $self, $taxid, $upi ) = @_;
  my @uniprot_acs = @{
	$self->{uniparc_dba}->dbc()->sql_helper()->execute_simple(
	  -SQL => q/
	  select ac
	  from uniparc.xref x 
	  join uniparc.cv_database d on (d.id=x.dbid)
	  where upi=? and taxid=? and descr in ('TREMBL','SWISSPROT') and deleted='N'/,
	  -PARAMS => [ $upi, $taxid ] ) };

  my $uniprots = [];

  for my $ac (@uniprot_acs) {
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
		my ( $ac, $name, $des, $type, $gene_name, $gene_name_type,
			 $version )
		  = @{ $_[0] };
		if ( defined $ac && $ac ne '' ) {
		  $uniprot->{ac} = $ac;
		}
		if ( defined $des && $des ne '' ) {
		  $uniprot->{description} = $des;
		}
		if ( defined $version && $version ne '' ) {
		  $uniprot->{version} = $version;
		}
		if ( defined $name && $name ne '' ) {
		  if ( defined $uniprot->{name} ) {
			$uniprot->{name} .= "; $name";
		  }
		  else {
			$uniprot->{name} = $name;
		  }
		}
		if ( defined $type && $type ne '' ) {
		  $uniprot->{type} =
			$type == 0 ? "Uniprot/SWISSPROT" : "Uniprot/SPTREMBL";
		}
		if ( defined $gene_name ) {
		  if ( $gene_name_type eq 'Name' ) {
			$uniprot->{gene_name} = $gene_name;
		  }
		  else {
			push @{ $uniprot->{synonyms} }, $gene_name;
		  }
		}
		return;
	  } );
	if ( defined $uniprot->{ac} ) {
	  push @$uniprots, $uniprot;
	}
  } ## end for my $ac (@uniprot_acs)
  return $uniprots;

} ## end sub get_uniprot_for_upi

=head2 remove_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Remove existing UniProt cross-references from genome
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub remove_xrefs {
  my ( $self, $dba ) = @_;
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
	-PARAMS => [ $dba->species_id() ] );
  return;
}

1;
