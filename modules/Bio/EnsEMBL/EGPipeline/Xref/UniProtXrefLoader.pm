
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

=cut

package Bio::EnsEMBL::EGPipeline::Xref::UniProtXrefLoader;
use base Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Digest::MD5;
use Data::Dumper;

sub new {
  my ($proto, @args) = @_;
  my $self = $proto->SUPER::new(@args);
  ($self->{uniprot_dba}, $self->{dbnames}) =
	rearrange(['UNIPROT_DBA', 'DBNAMES'], @args);
  if (!defined $self->{dbnames}) {
	$self->{dbnames} = qw/ArrayExpress PDB/;
  }
  $self->{dbnames} = {%hash = map { $_ => 1 } @{$self->{dbnames}}};
  return $self;
}

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

sub get_translation_uniprot {
  my ($self, $dba) = @_;
  # get hash of xrefs by translation ID
  my $translation_accs = {};
  my $dbea             = $dba->get_DBEntryAdaptor();
  $dba->dbc()->sql_helper()->execute_no_return(
	-SQL => q/
		select tl.translation_id,unix.xref_id
		from
	translation tl
	join transcript tr using (transcript_id)
	join seq_region sr using (seq_region_id)
	join coord_system cs using (coord_system_id)
	join object_xref uniox on (uniox.ensembl_object_type='Translation' and uniox.ensembl_id=tl.translation_id)
	join xref unix using (xref_id) 
	join external_db unie using (external_db_id) 
        where
	cs.species_id=? and 
	unie.db_name in ('Uniprot\/SWISSPROT','Uniprot\/SPTREMBL')
	/,
	-CALLBACK => sub {
	  my ($tid, $xid) = @{$_[0]};
	  $translation_accs->{$tid} = $dbea->fetch_by_dbID($xid);
	  return;
	},
	-PARAMS => [$dba->species_id()]);
  return $translation_accs;
} ## end sub get_translation_uniprot

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
}

sub store_xref {
  my ($self, $ddba, $tid, $uniprot) = @_;
  my @xrefs = grep { defined $self->{dbnames}{$_->{dbname}} }
	@{$self->get_xrefs_for_uniprot($uniprot->primary_id())};
  my $n = 0;
  for my $xref (@xrefs) {
	$n++;
	my $des;
	if (defined $xref->{secondary_id} &&
		$xref->{secondary_id} ne $xref->{primary_id})
	{
	  $des = $xref->{secondary_id};
	}
	$ddba->dbc()->sql_helper()->execute_update(
	  -SQL => q/delete ox.* from object_xref ox 
	join xref x using (xref_id) 
	join external_db e using (external_db_id) 
	where ox.ensembl_id=? and ox.ensembl_object_type='Translation'
	and e.db_name=? and x.dbprimary_acc=?/,
	  -PARAMS => [$tid, $xref->{dbname}, $xref->{primary_id}]);
	my $dbentry =
	  Bio::EnsEMBL::DBEntry->new(-DBNAME      => $xref->{dbname},
								 -PRIMARY_ID  => $xref->{primary_id},
								 -DISPLAY_ID  => $xref->{primary_id},
								 -DESCRIPTION => $des);
	$dbentry->analysis($self->{analysis});
	$ddba->store($dbentry, $tid, 'Translation');
  }
  return $n;
} ## end sub store_xref

sub get_xrefs_for_uniprot {
  my ($self, $ac) = @_;
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
}

1;
