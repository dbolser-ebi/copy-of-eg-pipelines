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
use Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor;
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

$logger->info("Loading ontology db");
my ($ont_dba_details) =
  @{ $cli_helper->get_dba_args_for_opts( $opts, 1, 'ont' ) };
my $ont_dba =
  new Bio::EnsEMBL::DBSQL::OntologyDBAdaptor(%$ont_dba_details);

$logger->info("Fetching ontologies");
# Get lists of term names for host, phenotype and condition

my $ont_adaptor = $ont_dba->get_adaptor("OntologyTerm");

my $ident_root_term = $ont_adaptor->fetch_by_accession("PHI:0");
my @ident_terms =
  @{ $ont_adaptor->fetch_all_by_ancestor_term($ident_root_term) };

Bio::EnsEMBL::Registry->load_registry_from_db( -USER => $opts->{user},
                         -PASS => $opts->{pass},
                         -HOST => $opts->{host},
                         -PORT => $opts->{port} );

print "Registry loaded\n";

$logger->info("Completed fetching PHI-base identifiers from ontology db");

# Get all phibase rows in the cvs file that are also in the core db

my $phibase_file = $opts->{file};
$logger->info("Reading $phibase_file");

open( my $INP, "<", $phibase_file ) or
  croak "Could not open $phibase_file for reading";

my %host_voc;
my %condition_voc;
my %phenotype_voc;

my $updated_dbcs = {};

# nested hash of phibase xrefs by species-translation_id-phiId
my $phibase_xrefs = {};
# hash of dbas by species

my $n = 0;
my $x = 0;
my $counter = 0;
LINE: while ( my $line = <$INP> ) {
  print "$line";
  next if $line =~ m/^phinum.*/;

  chomp $line;

  my $msg     = '';
  
  my @cols = split( ",", $line, -1 );

  my $phibase_id      = $cols[1];
  my $acc             = $cols[2];
  my $locus           = $cols[4];
  my $gene_name       = $cols[3];
  my $tax_id          = $cols[7];
  my $species_tax_id 	= $cols[6];
  my $species_name    = $cols[8];
  my $host_ids        = $cols[9];
  my $host_names      = $cols[10];
  my $phenotype_name  = $cols[11];
  my $condition_names = $cols[12];
  my $literature_ids  = $cols[13];
  my $dois            = $cols[14];

  my $direct_map			= $cols[16];
  my $mapping_method  = $cols[17];
  my $db_name 				= $cols[18];
  my $gene_stable_id 	= $cols[19];
  my $translation_stable_id = $cols[20]; 
 
  next if ! $translation_stable_id;
 
  #hack for botrytis
  if ( $db_name eq 'botrytis_cinerea_core_35_88_3' ) {
    $translation_stable_id =~ tr/g/p/; 
  }

  my $evidence;
  my $success;

  if ($mapping_method != 0) {
    $evidence = 'DIRECT';
  } else {
    $evidence = 'SEQUENCE_MATCH';
  }

  for my $var ( $phibase_id, $acc, $tax_id, $phenotype_name ) {
	 if ( !defined $var ) {
    print "problem with var\n" . $var ;
	  $success = 0;
	  $msg     = "Cannot parse line";
	  next LINE;
	 }
  }

  if ( $phibase_id =~ /;/ ) {
    next LINE;
  }

  print "Line is parsed\n";

  # get dbadaptor based on tax ID
  print "Processing entry: $phibase_id annotation on translation $translation_stable_id from $species_name ($species_tax_id) on db $db_name\n";

  my $dba;
  my $db_name_pref;
  my $species_id;
=cut
  if ( $db_name !~ /collection/ ) {
    print $db_name . "\n";
    $dbas{$db_name} = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -USER   => $opts->{user},
                     -PASS   => $opts->{pass},
                     -HOST   => $opts->{host},
                     -PORT   => $opts->{port},
                     -DBNAME => $db_name );
  } else {
    ( $db_name_pref , $species_id ) = split(/\//, $db_name );
    $dbas{$db_name} = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -USER   => $opts->{user},
                     -PASS   => $opts->{pass},
                     -HOST   => $opts->{host},
                     -PORT   => $opts->{port},
                     -DBNAME => $db_name_pref,
                     -SPECIES_ID => $species_id,
                     -SPECIES => $species_name,
                     -MULTISPECIES_DB => 1,
                      );
    
  } 
=cut
#Checking is translation exists in database

  my $translation_adaptor = Bio::EnsEMBL::Registry->get_adaptor($db_name, "Core", "Translation");

  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);

  if ($translation) {
    print "Found translation on database with dbID " . $translation->dbID() . "\n";
  } else {
    print "Could not find translation on database for $translation_stable_id for $db_name with species_id $species_id\n";
    exit;
  }
 
  my $translation_ass = {};

#Checking if host exist in the file

  $translation_ass->{host} =
      { id => $host_ids, label => $host_names };

  if ( !defined $translation_ass->{host}{id} ||
     !defined $translation_ass->{host}{label} )
  {
    print  "Host ID/label not defined for $phibase_id";
    next;
  } else {
    print "Found host for this annotation: " . $translation_ass->{host}{label} . "\n";
  } 

#Adding phenotype and condition to associated xref

  if ( !defined $phenotype_name || $phenotype_name eq "") {
    print  "Could not find phenotype\n";
    next;
  }
  else {
    $translation_ass->{phenotype} = {
                  id => "$phenotype_name",
                  label => "$phenotype_name" };
  }

#Adding experimental condition

  for my $condition_full_name ( split( /;/, lc( rm_sp($condition_names) ) ) ) {
    if ( ! ($condition_full_name eq "")) {
      push (@{$translation_ass->{condition}}, {id => $condition_full_name,
                                               label => $condition_full_name });
    }
  }

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
  push @{ $translation_ass->{evidence} }, $evidence;

  push @{ $phibase_xrefs->{ $db_name }->{ $translation->dbID() }->{$phibase_id} }, $translation_ass; 
 
  #$dba->remove_db_adaptor($db_name_pref);

}

close $INP;


if ( $opts->{write} ) {

  $logger->info("Storing PHIbase xrefs");
  # now apply the changes
  my $count_by_db;
  my $count;
  my $mapped_translations;
  while ( my ( $db_name, $translations ) = each %$phibase_xrefs ) {

    my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($db_name, "Core");
    my $dbc    = $dba->dbc();
    my $dbname = $dbc->dbname();

    $logger->info( "Removing existing annotations from " . $db_name );
    $dbc->sql_helper()->execute_update(
    -SQL => q/delete x.*,ox.*,ax.*,ag.*,oox.* from external_db e
    join xref x using (external_db_id)
    join object_xref ox using (xref_id)
    join associated_xref ax using (object_xref_id)
    join associated_group ag using (associated_group_id)
    join ontology_xref oox using (object_xref_id)
    where e.db_name='PHI'/ );

    $dbc->sql_helper()->execute_update( -SQL => q/delete from gene_attrib where attrib_type_id=358/ );

    $dbc->sql_helper()->execute_update( -SQL => q/delete from gene_attrib where attrib_type_id=317 and value='PHI'/ );
     
    my $dbentry_adaptor = $dba->get_DBEntryAdaptor();

    print "Storing xrefs for $db_name\n";
    my $group = 0;
    $count_by_db->{$db_name} = 0;
    while ( my ( $translation, $phis ) = each %$translations ) {
      $mapped_translations->{$translation} = 1;
      while ( my ( $phi, $asses ) = each %$phis ) { 

        print "Storing $phi for $db_name translation $translation\n";        
        $count_by_db->{$db_name} += 1;
        $count += 1;
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
          if ( defined $ass->{doi} && scalar( @{ $ass->{doi} } ) > 0 ) {
            $pub_name = 'DOI';
            $pubs     = $ass->{doi};
          } else {
            $pubs = ['ND'];
          }
        }
        my @conditions_db_entries;
        foreach my $cond ( @{$ass->{condition}} ) {
          my $condition_db_entry = Bio::EnsEMBL::DBEntry->new(
                 -PRIMARY_ID => $cond->{id},
                 -DBNAME     => 'PHIE',
                 -RELEASE    => 1,
                 -DISPLAY_ID => $cond->{label}
          );
          push (@conditions_db_entries, $condition_db_entry);
        }
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
                     -INFO_TYPE  => $ass->{evidence}[0] );
          $phi_dbentry->add_associated_xref( $phenotype_db_entry,
               $pub_entry, 'phenotype', $group, $rank++ );
          $phi_dbentry->add_associated_xref( $host_db_entry,
                  $pub_entry, 'host', $group, $rank++ );  

          foreach my $condition_db_entry ( @conditions_db_entries ) {
          $phi_dbentry->add_associated_xref( $condition_db_entry,
                  $pub_entry, 'experimental evidence',
                  $group, $rank++ );
          }
          $phi_dbentry->add_linkage_type( 'ND', $pub_entry );
        }
      } ## end for my $ass (@$asses)
      print "Storing " . $phi_dbentry->display_id() . " on translation " .
        $translation . " from " . $dba->species() .
        " from " . $db_name . "\n";
      $dbentry_adaptor->store( $phi_dbentry, $translation,
                 'Translation' );
      }  
			}

    print "Added $count_by_db->{$db_name} xrefs to $db_name\n";

    print "Applying colors to $db_name\n";
    
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
      } else {
        if ( $gene_color{ $row[0] } eq $color ) {
         return;
        } else {
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
  }
  print Dumper($count_by_db);
  print "Added a total of $count PHI-base xrefs to ";
  print scalar keys %{$mapped_translations} ;
  print " translations.\n";
}  

sub rm_sp {
  my $string = shift;
  $string =~ s/^\s*//;
  $string =~ s/\s*$//;
  return $string;
}






