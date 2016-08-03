
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

Bio::EnsEMBL::EGPipeline::Xref::UniProtXrefLoader

=head1 DESCRIPTION

Loader that adds new xrefs to translations based on existing UniProt cross-references

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::UniProtXrefLoader;
use base Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Digest::MD5;
use Data::Dumper;

=head1 CONSTRUCTOR
=head2 new
  Arg [-UNIPROT_DBA]  : 
       string - adaptor for UniProt Oracle database (e.g. SWPREAD)
  Arg [-DBNAMES]    : 
       array - array of database names to process (default is ArrayExpress, PDB, EMBL, ChEMBL)

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
  ($self->{uniprot_dba}, $self->{dbnames}) =
	rearrange(['UNIPROT_DBA', 'DBNAMES'], @args);
  if (!defined $self->{dbnames}) {
	$self->{dbnames} = qw/ArrayExpress ChEMBL EMBL MEROPS PDB/;
  }
  $self->{dbnames} = {%hash = map { $_ => 1 } @{$self->{dbnames}}};
  return $self;
}

=head1 METHODS
=head2 load_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Add xrefs to supplied core
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut
sub load_xrefs {
  my ($self, $dba) = @_;
  $self->{analysis} = $self->get_analysis($dba, 'xrefuniprot');
  # get translation_id,UniProt xref
  my $translation_uniprot = $self->get_translation_uniprot($dba);
  $self->logger()
	->info("Found " .
		   scalar(keys %$translation_uniprot) .
		   " translations with UniProt entries");
  $self->add_xrefs($dba, $translation_uniprot);
  $self->logger()->info("Finished loading xrefs");
  return;
}

=head2 add_xrefs
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : hashref of translation ID to UniProt accessions
  Description: Add xrefs to specified translations
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub add_xrefs {
  my ($self, $dba, $translation_uniprot) = @_;
  my $ddba = $dba->get_DBEntryAdaptor();
  my $tN   = 0;
  my $uN   = 0;
  $self->logger()->info("Adding xrefs");
  while (my ($tid, $uniprot) = each %$translation_uniprot) {
	$tN++;
	$uN += $self->store_xref($ddba, $tid, $uniprot);
	$uN += $n;
	$self->logger()->info("Processed $tN translations ($uN xrefs)")
	  if ($tN % 1000 == 0);
  }
  $self->logger()->info("Stored $uN xrefs on $tN translations");
  return;
}

=head2 store_xref
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : Translation dbID
  Arg        : Bio::EnsEMBL::DBEntry for UniProt record
  Description: Add xrefs to specified translation
  Returntype : number of xrefs attached
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub store_xref {
  my ($self, $ddba, $tid, $uniprot) = @_;
  # get xrefs we're interested in
  my @xrefs = @{$self->get_xrefs_for_uniprot($uniprot->primary_id())};
  
  ## Dump some info about what we're (potentially) ignoring...
  for( @xrefs ){
      $self->logger->debug(
          join("\t", $tid, $_->{DBNAME}, $_->{PRIMARY_ID})
      );
  }
  
  ## Trim out all but the list of what we want
  @xrefs = grep { defined $self->{dbnames}{$_->{DBNAME}} }
       @xrefs;
  
  # special rules for ENA - we don't want genomic references where we have the CDS
  @xrefs = 
	map {
	  if ($_->{DBNAME} eq 'EMBL' &&
		  defined $_->{QUATERNARY_ID} &&
		  defined $_->{SECONDARY_ID}  &&
		  $_->{QUATERNARY_ID} eq 'Genomic_DNA')
	  {
		$_->{PRIMARY_ID} = $_->{SECONDARY_ID};
		$_->{DBNAME} = 'protein_id';
	  }
	  $_;
	  }
	  @xrefs;	

  my $n = 0;
  for my $xref (@xrefs) {
	$n++;
	$self->logger()->debug("Attaching ".$xref->{DBNAME}.":".$xref->{PRIMARY_ID}." to translation ".$tid);
	$ddba->dbc()->sql_helper()->execute_update(
	  -SQL => q/delete ox.* from object_xref ox 
	join xref x using (xref_id) 
	join external_db e using (external_db_id) 
	where ox.ensembl_id=? and ox.ensembl_object_type='Translation'
	and e.db_name=? and x.dbprimary_acc=?/,
	  -PARAMS => [$tid, $xref->{DBNAME}, $xref->{PRIMARY_ID}]);
	my $dbentry =
	  Bio::EnsEMBL::DBEntry->new(-DBNAME      => $xref->{DBNAME},
								 -PRIMARY_ID  => $xref->{PRIMARY_ID},
								 -DISPLAY_ID  => $xref->{PRIMARY_ID});
	$dbentry->analysis($self->{analysis});
	$ddba->store($dbentry, $tid, 'Translation');
  }
  return $n;
} ## end sub store_xref

=head2 get_xrefs_for_uniprot
  Arg        : UniProt accession
  Description: Find xrefs from UniProt for given accession
  Returntype : Array of hashref of xrefs
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub get_xrefs_for_uniprot {
  my ($self, $ac) = @_;
  $self->logger()->debug("Getting xrefs for $ac");
  my $xrefs = $self->{uniprot_dba}->dbc()->sql_helper()->execute(
	-USE_HASHREFS => 1,
	-SQL          => q/
	    select 
		abbreviation as dbname, primary_id, secondary_id, note, quaternary_id
		from dbentry d,
		dbentry_2_database dd,
		 database_name db
		where d.dbentry_id = dd.dbentry_id
		and
		db.database_id=dd.database_id
		and
		d.accession=?/,
	-PARAMS => [$ac]);
  return $xrefs;
} ## end sub get_xrefs_for_uniprot

1;
