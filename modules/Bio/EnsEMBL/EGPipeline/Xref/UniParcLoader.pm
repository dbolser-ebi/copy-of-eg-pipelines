
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

Bio::EnsEMBL::EGPipeline::Xref::UniParcLoader

=head1 DESCRIPTION

Loader that adds UPI xrefs to translations based on checksums

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::UniParcLoader;
use base Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Digest::MD5;

=head1 CONSTRUCTOR
=head2 new
  Arg [-UNIPARC]  : 
       Bio::EnsEMBL::DBSQL::DBAdaptor - UniParc database
  Example    : $loader = Bio::EnsEMBL::EGPipeline::Xref::UniParcLoader->new(...);
  Description: Creates a new loader
  Returntype : Bio::EnsEMBL::EGPipeline::Xref::UniParcLoader
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut
sub new {
  my ($proto, @args) = @_;
  my $self = $proto->SUPER::new(@args);
  ($self->{uniparc_dba}) =
	rearrange(['UNIPARC_DBA'], @args);
  return $self;
}
=head1 METHODS
=head2 add_upis
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor (core)
  Arg        : (optional) set to preserve old UPIs
  Description: Adds UPIs to supplied core
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut
sub add_upis {
  my ($self, $dba, $preserve_old) = @_;
  if (!defined $preserve_old) {
	$self->remove_upis($dba);
  }
  # now work through transcripts
  my $translationN = 0;
  my $upiN         = 0;
  my $ddba         = $dba->get_DBEntryAdaptor();
  $self->logger()
	->info("Processing translations for " . $dba->species());
  my @translations =
	@{$dba->get_TranscriptAdaptor->fetch_all_by_biotype(
													 'protein_coding')};
  my $analysis = $self->get_analysis($dba, 'xrefchecksum');
  for my $transcript (@translations) {
	$translationN++;
	$self->logger()
	  ->info("Processed $translationN/" .
			 scalar @translations . " translations")
	  if ($translationN % 1000 == 0);
	$upiN +=
	  $self->add_upi($ddba, $transcript->translation(), $analysis);
  }
  $self->logger()
	->info("Stored UPIs for $upiN of $translationN translations");
  return;
} ## end sub add_upis

=head2 remove_upis
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor (core)
  Description: Remove existing UPIs from supplied core
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub remove_upis {
  my ($self, $dba) = @_;
  $self->logger()->info("Removing existing UniParc cross-references");
  $dba->dbc()->sql_helper()->execute_update(
	-SQL => q/
		delete ox.*
		from object_xref ox,
		xref x,
		external_db d,
		translation tl,
		transcript tr,
		seq_region sr,
		coord_system cs
		where
		d.db_name='UniParc'
		and d.external_db_id=x.external_db_id
		and x.xref_id = ox.xref_id
		and ox.ensembl_id = tl.translation_id
		and ox.ensembl_object_type = 'Translation'
		and tl.transcript_id=tr.transcript_id
		and tr.seq_region_id=sr.seq_region_id
		and sr.coord_system_id=cs.coord_system_id
		and cs.species_id=?
	/,
	-PARAMS => [$dba->species_id()]);
  return;
} ## end sub remove_upis

=head2 add_upi
  Arg        : Bio::EnsEMBL::DBSQL::DBEntryAdaptor (core)
  Arg        : Bio::EnsEMBL::Translation - translation to analysis
  Arg        : Bio::EnsEMBL::Analysis - analysis to use
  Description: Find and add UPI for supplied translation
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub add_upi {
  my ($self, $ddba, $translation, $analysis) = @_;
  my $stored = 0;
  $self->logger()
	->debug("Finding UPI for " . $translation->stable_id());
  my $hash = $self->md5_checksum($translation);
  my @upis = @{
	$self->{uniparc_dba}->dbc()->sql_helper->execute_simple(
				 -SQL => q/select upi from uniparc.protein where md5=?/,
				 -PARAMS => [$hash])};
  my $nUpis = scalar(@upis);
  if ($nUpis == 0) {
	$self->logger()
	  ->warn(
		   "No UPI found for translation " . $translation->stable_id());
  }
  elsif ($nUpis == 1) {
	$stored = 1;
	$self->logger()
	  ->debug("UPI $upis[0] found for translation " .
			  $translation->stable_id() . " - storing...");
	my $dbentry =
	  Bio::EnsEMBL::DBEntry->new(-PRIMARY_ID => $upis[0],
								 -DISPLAY_ID => $upis[0],
								 -DBNAME     => 'UniParc',
								 -INFO_TYPE  => 'CHECKSUM');
	$dbentry->analysis($analysis);
	$ddba->store($dbentry, $translation->dbID(), 'Translation');
  }
  else {
	$self->logger()
	  ->warn("Multiple UPIs found for translation " .
			 $translation->stable_id());
  }
  return $stored;
} ## end sub add_upi

=head2 md5_checksum
  Arg        : Bio::EnsEMBL::Translation
  Description: Calculate checksum for supplied sequence
  Returntype : string
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub md5_checksum {
  my ($self, $sequence) = @_;
  my $digest = Digest::MD5->new();
  $digest->add($sequence->seq());
  return uc($digest->hexdigest());
}

1;
