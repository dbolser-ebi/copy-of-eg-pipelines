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

Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO

=head1 DESCRIPTION

Runnable that invokes LoadUniProtGO on a core database

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::BacteriaUniProtLoader;
use base Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Digest::MD5;
use List::MoreUtils qw/uniq natatime/;
use Data::Dumper;

my @whitelist = qw/
AGD
ArrayExpress
BindingDB
BioCyc
BRENDA
BuruList
ChEMBL
CYGD
DIP
DisProt
DrugBank
EchoBASE
EcoGene
GlycoSuiteDB
GO
IntAct
IntEnz
LegioList
ListiList
MEROPS
MIM
MINT
MypuList
Orphanet
PDB
PeroxiBase
PharmGKB
PhotoList
PRIDE
PseudoCAP
REBASE
SGD
SMR
STRING
SubtiList
TAIR
TIGR
TRANSFAC
TubercuList
UniGene
UniPathway/;

=head1 CONSTRUCTOR
=head2 new
  Arg [-UNIPROT_DBA]  : 
       string - adaptor for UniProt Oracle database (e.g. SWPREAD)
  Arg [-REPLACE_ALL]    : 
       boolean - remove all GO references first
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
  (  $self->{uniprot_dba} )
	= rearrange(
				 [ 'UNIPROT_DBA' ],
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
  $self->{analysis} = $self->get_analysis( $dba, 'xrefpid' );
  # get translation_id,PIDs where UniProt not set and PID is set
  my $translation_pids = $self->get_translation_pids($dba);
  $self->add_uniprot_xrefs( $dba, $translation_pids );
  return;
}

=head2 get_translation_uniprot
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Get mapping of UniParc to translation
  Returntype : Hash of translation IDs to UniParc accessions
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub get_translation_pids {
  my ( $self, $dba ) = @_;
  $self->logger()->info("Finding translation-PID pairs");

  my $translation_pids = {};
  $dba->dbc()->sql_helper()->execute_no_return(
	-SQL => q/
		select tl.translation_id,upix.dbprimary_acc
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
	upie.db_name='protein_id'
        and translation_id not in (
                select ensembl_id from object_xref sox 
                join xref sx using (xref_id) 
                where external_db_id in (2200,2000) and 
                ensembl_object_type='Translation')
	/,
	-CALLBACK => sub {
	  my ( $tid, $pid ) = @{ $_[0] };
	  $translation_pids->{$tid} = $pid;
	  return;
	},
	-PARAMS => [ $dba->species_id() ] );

  $self->logger()
	->info( "Found " .
		 scalar( keys %$translation_pids ) . " translation-PID pairs" );

  return $translation_pids;
} ## end sub get_translation_upis

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
  my ( $self, $dba, $translation_pids ) = @_;

  my $taxid = $dba->get_MetaContainer()->get_taxonomy_id();
  my $ddba  = $dba->get_DBEntryAdaptor();
  my $gdba  = $dba->get_GeneAdaptor();
  my $tN    = 0;
  my $uN    = 0;
  $self->logger()->info("Adding UniProt xrefs");
  # hash of gene names, descriptions
  my $uniprots_for_pids = $self->get_uniprot_for_pids( [values %$translation_pids] );
  while ( my ( $tid, $pid ) = each %$translation_pids ) {
	$tN++;
	$self->logger()->debug( "Looking up entry for $tid/$pid" );        
        my $uniprots = $uniprots_for_pids->{$pid};
	$uN +=
            $self->store_uniprot_xrefs( $ddba, $tid, $uniprots);
	$self->logger()->info("Processed $tN translations ($uN xrefs)")
	  if ( $tN % 1000 == 0 );
  }
  $self->logger()->info("Stored $uN UniProt xrefs on $tN translations");
  return;
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
  my ( $self, $ddba, $tid, $uniprots ) = @_;
  my $n = 0;
  my $m = 0;
  for my $uniprot(values %$uniprots) { 
      if ( !$uniprot->{acc} || $uniprot->{acc} eq '' ) {
	  $self->logger()
              ->warn(
              "Empty $uniprot->{type} accession retrieved from UniProt for translation $tid"
              );
	  next;
      }
      $self->logger()
          ->debug( "Storing $uniprot->{type} " .
                   $uniprot->{acc} . " on translation $tid" );
      my $dbentry =
	  $ddba->fetch_by_db_accession( $uniprot->{type}, $uniprot->{ac} );
      if ( !defined $dbentry ) {
	  $dbentry = Bio::EnsEMBL::DBEntry->new(
              -PRIMARY_ID  => $uniprot->{acc},
              -DISPLAY_ID  => $uniprot->{name},
              -DESCRIPTION => $uniprot->{description},
              -VERSION     => $uniprot->{version},
              -DBNAME      => $uniprot->{type} );
      }
      $dbentry->analysis( $self->{analysis} );
#	$ddba->store( $dbentry, $tid, 'Translation' );
      my @xrefs = grep { defined $self->{dbnames}{$_->{dbname}} }
      @{$uniprot->{xrefs}};
      
      for my $xref (@xrefs) {
          $m++;
          $self->logger()->debug("Attaching ".$xref->{dbname}.":".$xref->{primary_id}." to translation ".$tid);
          if($xref->{dbname} eq 'GO') {
              	my $go_xref =
                    Bio::EnsEMBL::OntologyXref->new(-dbname     => 'GO',
                                                    -primary_id => $xref->{primary_id},
                                                    -DISPLAY_ID => $xref->{primary_id});
                $xref->analysis($self->{analysis});
                my $linkage_type = $xref->{note1};
                if ($linkage_type) {
                    $linkage_type =~ s/:.*//;
                    $go_xref->add_linkage_type($linkage_type, $dbentry);
                }
                #$ddba->store($go_xref, $tid, 'Translation');
          } else {
              $ddba->dbc()->sql_helper()->execute_update(
                  -SQL => q/delete ox.* from object_xref ox 
	join xref x using (xref_id) 
	join external_db e using (external_db_id) 
	where ox.ensembl_id=? and ox.ensembl_object_type='Translation'
	and e.db_name=? and x.dbprimary_acc=?/,
                  -PARAMS => [$tid, $xref->{dbname}, $xref->{primary_id}]);
              my $dbentry =
                  Bio::EnsEMBL::DBEntry->new(-dbname      => $xref->{dbname},
                                             -primary_id  => $xref->{primary_id},
                                             -DISPLAY_ID  => $xref->{primary_id});
              $dbentry->analysis($self->{analysis});
              #$ddba->store($dbentry, $tid, 'Translation');
          }
      }

      $n++;
      

  }

#  for my $uniprot (@$uniprots) {
#
#  } ## end for my $uniprot (@$uniprots)
  return $n;
} ## end sub store_uniprot_xrefs


=head2 get_uniprot_for_pid
  Arg        : Taxonomy ID
  Arg        : UniParc identifier
  Description: Find UniProt accession for PID and taxonomy
  Returntype : hashref of matching UniProt records
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub get_uniprot_for_pids {
    my ( $self, $pids ) = @_;
    my $whitelist_str = join ',', map {"'$_'"} @whitelist;
  my $it = natatime 500, @$pids;
    my $uniprots = {};
  while (my @pids_sub = $it->())
  {
    my $pids_str =join ',', map {"'$_'"} @pids_sub; 
    my $sql = qq/
	  SELECT p.protein_id, d.accession,
  d.name,  REPLACE(NVL(sc1.text,sc3.text),'^'),
  d.entry_type,  s.version,
		db.abbreviation as dbname, ddd.primary_id, ddd.secondary_id, ddd.note, ddd.quaternary_id
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
JOIN embl_protein_id p on (p.dbentry_id=d.dbentry_id)
left JOIN dbentry_2_database ddd ON (ddd.dbentry_id=d.dbentry_id)
left JOIN database_name db ON (ddd.database_id=db.database_id)
WHERE d.entry_type in (0,1) and p.protein_id in ($pids_str)
and db.abbreviation in ($whitelist_str)
/;
    $self->{uniprot_dba}->dbc()->sql_helper()->execute_no_return(
        -SQL => $sql,
        -CALLBACK => sub {
            my ( $pid, $acc, $name, $des, $type, $version, $dbname, $primary_id, $secondary_id, $note, $quaternary_id )
                = @{ $_[0] };
            my $uniprot = $uniprots->{$pid}->{$acc};
            if(!defined $uniprot) {
                $uniprot = {
                    acc=>$acc,
                    description=>$des,
                    name=>$name,
                    version=>$version,
                    xrefs=>[]
                };                
                if ( defined $type && $type ne '' ) {
                    $uniprot->{type} =
                        $type == 0 ? "Uniprot/SWISSPROT" : "Uniprot/SPTREMBL";                    
                }
                $uniprots->{$pid}->{$acc} = $uniprot;
            }
            push @{$uniprot->{xrefs}}, {dbname=>$dbname, primary_id=>$primary_id,secondary_id=>$secondary_id,note=>$note, quaternary_id=>$quaternary_id}; 
            return;
        } );
  }  
    return $uniprots;
} ## end sub get_uniprot_for_pid

1;
