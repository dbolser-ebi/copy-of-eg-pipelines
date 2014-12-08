# Copyright [1999-2014] EMBL-European Bioinformatics Institute
# and Wellcome Trust Sanger Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::OntologyXref;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Bio::EnsEMBL::Utils::CliHelper;
use Carp;
use Data::Dumper;
use Log::Log4perl qw(:easy);

my $logger = get_logger();

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = [ @{ $cli_helper->get_dba_opts() },
			  @{ $cli_helper->get_dba_opts('ont') },
			  @{ $cli_helper->get_dba_opts('tax') } ];

push( @{$optsd}, "file:s" );
push( @{$optsd}, "verbose" );
push( @{$optsd}, "write" );
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&pod2usage );
if ( $opts->{verbose} ) {
  Log::Log4perl->easy_init($DEBUG);
}
else {
  Log::Log4perl->easy_init($INFO);
}

if ( defined $opts->{dbname} ) {
  $logger->info( "Loading " . $opts->{dbname} );
  Bio::EnsEMBL::DBSQL::DBAdaptor->new( -USER   => $opts->{user},
									   -PASS   => $opts->{pass},
									   -HOST   => $opts->{host},
									   -PORT   => $opts->{port},
									   -DBNAME => $opts->{dbname} );
}
else {
  $logger->info("Loading registry");
  Bio::EnsEMBL::Registry->load_registry_from_db( -USER => $opts->{user},
												 -PASS => $opts->{pass},
												 -HOST => $opts->{host},
												 -PORT => $opts->{port}
  );
}
$logger->info("Loading helper");
my $lookup =
  Bio::EnsEMBL::LookUp::LocalLookUp->new( -SKIP_CONTIGS => 1,
										  -NO_CACHE     => 1 );

my ($ont_dba_details) =
  @{ $cli_helper->get_dba_args_for_opts( $opts, 1, 'ont' ) };
my $ont_dba =
  new Bio::EnsEMBL::DBSQL::OntologyDBAdaptor(%$ont_dba_details);

# Get lists of term names for host, phenotype and condition

my $ont_adaptor = $ont_dba->get_adaptor("OntologyTerm");

my $ident_root_term = $ont_adaptor->fetch_by_accession("PHI:0");

my @ident_terms =
  @{ $ont_adaptor->fetch_all_by_ancestor_term($ident_root_term) };

my $pheno_root_term = $ont_adaptor->fetch_by_accession("PHI:001");

my @pheno_terms =
  @{ $ont_adaptor->fetch_all_by_ancestor_term($pheno_root_term) };

my $cond_root_term = $ont_adaptor->fetch_by_accession("PHI:002");

my @cond_terms =
  @{ $ont_adaptor->fetch_all_by_ancestor_term($cond_root_term) };

#TODO replace with taxonomy adaptor call
my %ncbi_host = ( 'wheat'                  => '4565',
				  'arabidopsis thaliana'   => '3702',
				  'thale cress'            => '3702',
				  'barley'                 => '112509',
				  'maize'                  => '4577',
				  'oat'                    => '4498',
				  'rye'                    => '39947',
				  'rice'                   => '39947',
				  'tomato'                 => '4081',
				  'cabbage'                => '51351',
				  'melon'                  => '3656',
				  'brassica juncea'        => '3707',
				  'canola'                 => '138011',
				  'rape seed oil'          => '3708',
				  'american chestnut tree' => '134033',
				  'apple'                  => '3750',
				  'pea'                    => '3888',
				  'chickpea'               => '3827',
				  'nicotiana benthamiana'  => '4100',
				  'tobacco'                => '4097',
				  'solanum huancabambense' => '856898',
				  'solanum microdontum'    => '73574',
				  'soybean'                => '3847',
				  'parsley'                => '4043',
				  'bean'                   => '3885',
				  'gerbera'                => '13546',
				  'grape'                  => '29760', );

my %ncbi_host_ids = reverse %ncbi_host;

# Get all phibase rows in the cvs file that are also in the core db

my $phibase_file = $opts->{file};

open( my $INP, "<", $phibase_file ) or
  croak "Could not open $phibase_file for reading";

my %host_voc;
my %condition_voc;
my %phenotype_voc;

my $updated_dbcs = {};

# nested hash of phibase xrefs by species-translation_id-phiId
my $phibase_xrefs = {};
# hash of dbas by species
my $dbas_by_genome = {};

open my $ok_outfile, ">",
  "mapped_phibase.txt" || croak "Could not open mapped_phibase.txt";
open my $fail_outfile, ">",
  "unmapped_phibase.txt" || croak "Could not open unmapped_phibase.txt";
my $n = 0;
my $x = 0;
LINE: while ( my $line = <$INP> ) {

  next if $line =~ m/^PHI-base accession no.*/;

  chomp $line;

  my $msg     = '';
  my $success = 0;
  my $found   = 0;

  # 0 PHI-base accession no,
  # 1 PHI-base accession,
  # 2 Obsolete PHI accession,
  # 3 DB_Type
  # 4 Accession,
  # 5 Obsolete EMBL accession,
  # 6 Locus ID,
  # 7 AA sequence #no EMBL#,
  # 8 NT sequence #no EMBL#,
  # 9 Associated strain,
  # 10 Gene name,
  # 11 Genome location,
  # 12 Multiple mutation,
  # 13 Pathogen NCBI Taxonomy ID,
  # 14 Pathogen species,
  # 15 Strain,
  # 16 Disease name,
  # 17 Monocot/Dicot plant,
  # 18 Host NCBI Taxonomy ID,
  # 19 Experimental host,
  # 20 Function,
  # 21 GO annotation,
  # 22 Database,
  # 23 Pathway,
  # 24 Phenotype of mutant,
  # 25 Mating defect prior to penetration,
  # 26 Pre-penetration defect,
  # 27 Penetration defect,
  # 29 Post-penetration defect,
  # 29 Vegetative spores,
  # 30 Sexual spores,
  # 31 In vitro growth,
  # 32 Spore germination,
  # 33 Essential gene (Lethal knockout),
  # 34 Inducer,
  # 35 CAS,
  # 36 Host target,
  # 37 Host response,
  # 38 Experimental evidence,
  # 39 Species Expert,
  # 40 Entered by,
  # 41 Manual (M) or textmining (T),
  # 42 Literature ID,
  # 43 Literature source,
  # 44 DOI,
  # 45 Full citation,
  # 46 Author email,
  # 47 Comments,
  # 48 Reference
  # 49 Year published
  my @cols = split( "\t", $line );

  my $phibase_id      = $cols[1];
  my $db_name         = $cols[3];
  my $acc             = $cols[4];
  my $locus           = $cols[6];
  my $gene_name       = $cols[10];
  my $tax_id          = $cols[13];
  my $species_name    = $cols[14];
  my $host_ids        = $cols[18];
  my $host_names      = $cols[19];
  my $phenotype_name  = $cols[24];
  my $condition_names = $cols[38];
  my $literature_ids  = $cols[42];
  my $dois            = $cols[44];

  for my $var ( $phibase_id, $acc, $tax_id, $phenotype_name ) {
	if ( !defined $var ) {
	  $success = 0;
	  $msg     = "Cannot parse line";
	  next LINE;
	}
  }

  $locus =~ s/\..*//;
  $acc       = rm_sp($acc);
  $locus     = rm_sp($locus);
  $gene_name = rm_sp($gene_name);
  # get dbadaptor based on tax ID
  $logger->debug(
"Processing entry: $phibase_id annotation on $db_name:$acc (gene $locus $gene_name) from $species_name ($tax_id)"
  );
  $logger->debug("Getting DBAs for species '$species_name' ($tax_id)");

  my $dbas;
  if ( defined $tax_id ) {
	$dbas = $lookup->get_all_by_taxon_branch($tax_id);
  }

  if ( !defined $dbas || scalar(@$dbas) == 0 ) {
	$msg = "No DBA for for taxon $tax_id (name $species_name)";
	$logger->warn($msg);
  }

  my $dba;

  for $dba ( @{$dbas} ) {

	$dbas_by_genome->{ $dba->species() } = $dba;

	$logger->debug( "Found DBA species " . $dba->species() .
			  " in " . $dba->dbc()->dbname() . "/" . $dba->species_id );
	my @transcripts_ids;
	my @gene_ids;

	my $translation =
	  find_translation( $dba, $acc, $locus, $gene_name );

	if ( !defined $translation ) {
	  $msg =
"Failed to find translation for $db_name:$acc (gene $locus $gene_name) from $species_name ($tax_id) in "
		. $dba->species()
		. " jn " . $dba->dbc()->dbname() . "/" . $dba->species_id();
	  $logger->warn($msg);
	  next;
	}
	$found = 1;
	$logger->debug( "Found translation " .
			   $translation->stable_id() . "/" . $translation->dbID() );

	my $translation_ass = {};

	my @phibase_dbentries;
	my @phibase_types;
	$logger->debug("Processing host(s) $host_ids/$host_names");
	my @host_names = split( /;/, $host_names );
	my @host_ids   = split( /;/, $host_ids );
	if ( scalar(@host_ids) == 0 && scalar(@host_names) > 0 ) {
	  @host_ids = (undef) x scalar(@host_names);
	}
	my $hN = 0;
	for my $host_id (@host_ids) {
	  my $host_tax_id;
	  my $host_name = rm_sp( $host_names[$hN] );
	  if ( defined $host_id && $host_id =~ m/^[0-9]+$/ ) {
		$host_tax_id = $host_id;
	  }
	  if ( !defined $host_tax_id || $host_tax_id eq '' ) {
		$host_tax_id = $ncbi_host{ lc $host_name };
		if ( defined $host_tax_id ) {
		  $logger->debug(
			 "found host term in the file that is on known hosts list: "
			   . $host_tax_id . "/" . $host_name );
		}
	  }

	  if ( !defined $host_name || $host_name eq '' ) {
		$host_name = $ncbi_host_ids{$host_tax_id};
		if ( defined $host_name ) {
		  $logger->debug(
			 "found host name in the file that is on known hosts list: "
			   . $host_name );
		}
	  }

	  if ( defined $host_tax_id &&
		   defined $host_name &&
		   $host_tax_id ne '' &&
		   $host_name ne '' )
	  {

		$translation_ass->{host} =
		  { id => $host_tax_id, label => $host_name };

	  }
	  else {
		$logger->warn("Could not find host $host_id");
	  }
	  $hN++;
	} ## end for my $host_id (@host_ids)

	if ( !defined $translation_ass->{host}{id} ||
		 !defined $translation_ass->{host}{label} )
	{
	  my $msg = "Host ID/label not defined for $phibase_id";
	  $logger->warn($msg);
	  next;
	}

	my $found_phenotype;
	$logger->debug("Processing phenotype '$phenotype_name'")
	  ;
	for my $phenotype (@pheno_terms) {
	  my $ont_phenotype_name = lc( rm_sp( $phenotype->name() ) );
	  if ( $phenotype_name =~ /$ont_phenotype_name/i ) {
		$logger->debug( "Mapped '$phenotype_name' to term " .
						$phenotype->accession() );
		$found_phenotype = $phenotype;
		next;
	  }
	}

	if ( !defined $found_phenotype ) {
	  $msg = "Could not find phenotype '$phenotype_name'";
	  $logger->warn($msg);
	  next;
	}
	else {
	  $translation_ass->{phenotype} = {
									id => $found_phenotype->accession(),
									label => $found_phenotype->name() };
	}

	#deal with the conditions
	my %fix_cond_divergences = (
							'complementation' => 'gene complementation',
							'mutation'        => 'gene mutation', );
	$logger->debug("Processing condition(s) '$condition_names'");
	for my $condition_full_name (
						   split( /;/, lc( rm_sp($condition_names) ) ) )
	{
	  my @condition_short_names = split( /:/, $condition_full_name );
	  my $condition_last_name = rm_sp( $condition_short_names[-1] );
	  $condition_last_name = $fix_cond_divergences{$condition_last_name}
		if exists $fix_cond_divergences{$condition_last_name};
	  for my $ont_cond (@cond_terms) {
		my $ont_cond_name = lc( $ont_cond->name() );
		if ( $condition_last_name =~ /$ont_cond_name/ ) {
		  $logger->debug( "Mapped condition to " .
					 $ont_cond->accession() . "/" . $ont_cond->name() );

		  $translation_ass->{condition} = {id => $ont_cond->accession(),
										   label => $ont_cond->name() };
		  last;
		}
	  }
	}

	#deal with the publications
	if ( defined $literature_ids ) {
	  $logger->debug("Processing literature refs(s) '$literature_ids'");
	  for my $publication ( split( /;/, $literature_ids ) ) {
		push @{ $translation_ass->{pubmed} }, $publication;
	  }
	}
	if ( defined $dois ) {
	  $logger->debug("Processing literature ref(s) '$dois'");
	  for my $publication ( split( /;/, $dois ) ) {
		push @{ $translation_ass->{doi} }, $publication;
	  }
	}

	if ( $success == 0 ) {
	  $success = 1;
	}

	push
	  @{ $phibase_xrefs->{ $dba->species() }->{ $translation->dbID() }
		->{$phibase_id} }, $translation_ass;

	print $ok_outfile join( "\t",
							$dba->species(), $translation->stable_id(),
							$phibase_id ) .
	  "\n";

  } ## end for $dba ( @{$dbas} )
  continue {
	if ( defined $dba ) {
	  $dba->dbc()->disconnect_if_idle();
	}
  }

  if ( $success == 1 ) {
	$n++;
  }
  else {
	if ( $found == 0 ) {
	  $msg = "Gene not found in any target genomes";
	}
	print $fail_outfile "# $msg\n";
	print $fail_outfile "$line\n";
	$x++;
  }

} ## end LINE: while ( my $line = <$INP>)
$logger->info("Completed - matched $n and skipped $x xrefs");
close $fail_outfile;
close $ok_outfile;
close $INP;

if ( $opts->{write} ) {

  my $dbcs = {};

  $logger->info("Storing PHIbase xrefs");
  # now apply the changes
  while ( my ( $genome, $translations ) = each %$phibase_xrefs ) {

	my $dba    = $dbas_by_genome->{$genome};
	my $dbc    = $dba->dbc();
	my $dbname = $dbc->dbname();

	if ( !defined $dbcs->{$dbname} ) {
	  $logger->info( "Removing existing annotations from " . $dbname );
	  $dbc->sql_helper()->execute_update(
		-SQL => q/delete x.*,ox.*,ax.*,ag.*,oox.* from external_db e 
join xref x using (external_db_id) 
join object_xref ox using (xref_id) 
join associated_xref ax using (object_xref_id) 
join associated_group ag using (associated_group_id) 
join ontology_xref oox using (object_xref_id) 
where e.db_name='PHI'/ );

	  $dbc->sql_helper()
		->execute_update(
		  -SQL => q/delete from gene_attrib where attrib_type_id=358/ );

	  $dbc->sql_helper()
		->execute_update( -SQL =>
q/delete from gene_attrib where attrib_type_id=317 and value='PHI'/ );

	  $dbcs->{$dbname} = $dba->dbc();
	}

	my $dbentry_adaptor = $dba->get_DBEntryAdaptor();

	$logger->info("Storing xrefs for $genome");
	my $group = 0;
	my $tN    = 0;
	my $xN    = 0;
	while ( my ( $translation, $phis ) = each %$translations ) {
	  $tN++;
	  $logger->info(
				  "Storing xrefs for $genome translation $translation");
	  while ( my ( $phi, $asses ) = each %$phis ) {
		$xN++;
		$logger->info(
				   "Storing $phi for $genome translation $translation");
		my $phi_dbentry =
		  Bio::EnsEMBL::OntologyXref->new( -PRIMARY_ID  => $phi,
										   -DBNAME      => 'PHI',
										   -DISPLAY_ID  => $phi,
										   -DESCRIPTION => $phi,
										   -RELEASE     => 1,
										   -INFO_TYPE   => 'DIRECT' );
		for my $ass (@$asses) {
		  # first, work out literature
		  my $pub_name = 'PUBMED';
		  my $pubs     = $ass->{pubmed};
		  if ( !defined $pubs || scalar( @{$pubs} ) == 0 ) {
			if ( defined $ass->{doi} && scalar( @{ $ass->{doi} } ) > 0 )
			{
			  $pub_name = 'DOI';
			  $pubs     = $ass->{doi};
			}
			else {
			  $pubs = ['ND'];
			}
		  }
		  my $condition_db_entry =
			Bio::EnsEMBL::DBEntry->new(
								 -PRIMARY_ID => $ass->{condition}{id},
								 -DBNAME     => 'PHIE',
								 -RELEASE    => 1,
								 -DISPLAY_ID => $ass->{condition}{label}
			);
		  my $host_db_entry =
			Bio::EnsEMBL::DBEntry->new(
									  -PRIMARY_ID => $ass->{host}{id},
									  -DBNAME     => 'NCBI_TAXONOMY',
									  -RELEASE    => 1,
									  -DISPLAY_ID => $ass->{host}{label}
			);
		  my $phenotype_db_entry =
			Bio::EnsEMBL::DBEntry->new(
								 -PRIMARY_ID => $ass->{phenotype}{id},
								 -DBNAME     => 'PHIP',
								 -RELEASE    => 1,
								 -DISPLAY_ID => $ass->{phenotype}{label}
			);

		  my $rank = 0;
		  for my $pub (@$pubs) {
			$group++;
			my $pub_entry =
			  Bio::EnsEMBL::DBEntry->new(
									   -PRIMARY_ID => lc( rm_sp($pub) ),
									   -DBNAME     => $pub_name,
									   -DISPLAY_ID => lc( rm_sp($pub) ),
									   -INFO_TYPE  => 'DIRECT' );
			$phi_dbentry->add_associated_xref( $phenotype_db_entry,
							 $pub_entry, 'phenotype', $group, $rank++ );
			$phi_dbentry->add_associated_xref( $host_db_entry,
								  $pub_entry, 'host', $group, $rank++ );
			$phi_dbentry->add_associated_xref( $condition_db_entry,
									$pub_entry, 'experimental evidence',
									$group, $rank++ );

			$phi_dbentry->add_linkage_type( 'ND', $pub_entry );
		  }
		} ## end for my $ass (@$asses)
		$logger->debug(
		  "Storing " . $phi_dbentry->display_id() . " on translation " .
			$translation . " from " . $dba->species() .
			" from " . $dba->dbc()->dbname() . "/" . $dba->species_id );
		$dbentry_adaptor->store( $phi_dbentry, $translation,
								 'Translation' );
	  } ## end while ( my ( $phi, $asses...))
	} ## end while ( my ( $translation...))
	$logger->info( "Stored " . $xN .
		 " xrefs on " . $tN . " translations from " . $dba->species() );
  } ## end while ( my ( $genome, $translations...))

  # now add the colours

  # last step - apply the colours
  for my $dbc ( values %{$dbcs} ) {
	$logger->info( "Applying colours for " . $dbc->dbname() );
# assemble a list of gene colours, dealing with mixed_outcome as required
	my %gene_color;
	$dbc->sql_helper()->execute_no_return(
	  -SQL => q/select t.gene_id, x.display_label 
	from associated_xref a, object_xref o, xref x, transcript t, translation tl
	where a.object_xref_id = o.object_xref_id and condition_type = 'phenotype' 
        and tl.transcript_id=t.transcript_id
	and x.xref_id = a.xref_id and o.ensembl_id = tl.translation_id/,
	  -CALLBACK => sub {
		my @row   = @{ shift @_ };
		my $color = lc( $row[1] );
		$color =~ s/^\s*(.*)\s*$/$1/;
		$color =~ s/\ /_/g;
		if ( !exists( $gene_color{ $row[0] } ) ) {
		  $gene_color{ $row[0] } = $color;
		}
		else {
		  if ( $gene_color{ $row[0] } eq $color ) {
			return;
		  }
		  else {
			$gene_color{ $row[0] } = 'mixed_outcome';
		  }
		}
		return;
	  } );
	foreach my $gene ( keys %gene_color ) {
	  $logger->debug("Setting " . $gene . " as " . $gene_color{$gene} );
	  $dbc->sql_helper()->execute_update(
		-SQL =>
q/INSERT INTO gene_attrib (gene_id, attrib_type_id, value) VALUES ( ?, 358, ?)/,
		-PARAMS => [ $gene, $gene_color{$gene} ] );
	  $dbc->sql_helper()->execute_update(
		-SQL =>
q/INSERT INTO gene_attrib (gene_id, attrib_type_id, value) VALUES ( ?, 317, 'PHI')/,
		-PARAMS => [$gene] );
	}
  } ## end for my $dbc ( values %{...})

} ## end if ( $opts->{write} )

sub rm_sp {
  my $string = shift;
  $string =~ s/^\s*//;
  $string =~ s/\s*$//;
  return $string;
}

sub find_translation {
  my ( $dba, $acc, $locus, $gene_name ) = @_;
  my $translation;
  my $dbentry_adaptor = $dba->get_adaptor("DBEntry");
  my @transcripts_ids =
	$dbentry_adaptor->list_transcript_ids_by_extids($acc);
  my @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($acc);

  if ( scalar(@gene_ids) == 0 ) {
	@gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($locus);
  }
  if ( scalar(@gene_ids) == 0 ) {
	@gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($gene_name);
  }
  my $translation_adaptor = $dba->get_adaptor("Translation");
  my $transcript_adaptor  = $dba->get_adaptor("Transcript");
  my $gene_adaptor        = $dba->get_adaptor("Gene");
  my $transcript;
  my $gene;

  if ( scalar(@transcripts_ids) >= 1 ) {
	my $transcript_id = $transcripts_ids[0];
	$transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
	$translation =
	  $translation_adaptor->fetch_by_Transcript($transcript);
  }
  elsif ( scalar(@gene_ids) >= 1 ) {
	$gene = $gene_adaptor->fetch_by_dbID( $gene_ids[0] );
	my @transcripts =
	  @{ $transcript_adaptor->fetch_all_by_Gene($gene) };
	$transcript = $transcripts[0];
	$translation =
	  $translation_adaptor->fetch_by_Transcript($transcript);
  }
  return $translation;
} ## end sub find_translation
