# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

# Sets display_xref_ids for novel genes in the "to" database based
# on their orthologs in the "from" database. Can also project GO xrefs.
# Orthology relationships are read from a Compara database.

use Cwd qw/cwd/;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end eprof_dump);
use Bio::EnsEMBL::Utils::Exception;

use LWP::Simple;
use JSON;

my $method_link_type = "ENSEMBL_ORTHOLOGUES";

my ($from_species, @to_multi, $release);

my $host;
my $user = 'ensro';
my $pass;
my $port = 3306;
my $host1 = 'ens-staging';
my $user1 = 'ensro';
my $pass1;
my $port1 = 3306;
my $host2 = 'ens-staging2';
my $user2 = 'ensro';
my $pass2;
my $port2 = 3306;
my $compara_host = 'ens-livemirror';
my $compara_user = 'ensro';
my $compara_pass;
my $compara_port = 3306;
my $compara_dbname;

GetOptions('host=s'            => \$host,
           'user=s'            => \$user,
           'pass=s'            => \$pass,
           'port=s'            => \$port,
           'host1=s'           => \$host1,
           'user1=s'           => \$user1,
           'pass1=s'           => \$pass1,
           'port1=s'           => \$port1,
           'host2=s'           => \$host2,
           'user2=s'           => \$user2,
           'pass2=s'           => \$pass2,
           'port2=s'           => \$port2,
           'compara_host=s'    => \$compara_host,
           'compara_user=s'    => \$compara_user,
           'compara_pass=s'    => \$compara_pass,
           'compara_port=s'    => \$compara_port,
           'compara_dbname=s'  => \$compara_dbname,
	   'from=s'            => \$from_species,
	   'to=s'              => \@to_multi,
	   'release=i'         => \$release,
	   'help'              => sub { usage(); exit(0); });

$| = 1; # auto flush stdout

@to_multi = split(/,/,join(',',@to_multi));

# load from database and conf file
my $registry = "Bio::EnsEMBL::Registry";
$registry->no_version_check(1);

# Registryconf is either the registry configuration passed from the submit_projections.pl 
# script or a file name containing the same information that is passed on the command line.

my $args;
if ($host) {
  $registry->load_registry_from_multiple_dbs(
          {
              '-host'       => $host,
              '-user'       => $user,
              '-pass'       => $pass,
              '-port'       => $port,
              '-db_version' => $release,
          },
  );
} else {
  $registry->load_registry_from_multiple_dbs(
          {
              '-host'       => $host1,
              '-user'       => $user1,
              '-pass'       => $pass1,
              '-port'       => $port1,
              '-db_version' => $release,
          },
          {
              '-host'       => $host2,
              '-user'       => $user2,
              '-pass'       => $pass2,
              '-port'       => $port2,
              '-db_version' => $release,
          }
  );
}

my $compara_db;
if ($compara_dbname) {
  if($registry->get_DBAdaptor('multi', 'compara', 1)) {
    $registry->remove_DBAdaptor('multi', 'compara');
  }
  $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
              '-host'    => $compara_host,
              '-user'    => $compara_user,
              '-pass'    => $compara_pass,
              '-port'    => $compara_port,
              '-dbname'  => $compara_dbname,
              '-species' => 'multi',
              '-group'  => 'compara',
  );
}


# Get Compara adaptors - use the one specified on the command line, or the first one
# defined in the registry file if not specified

my $mlssa;
my $ha;
my $ma;
my $gdba;

if ($compara_db) {

   $mlssa = $compara_db->get_adaptor('MethodLinkSpeciesSet');
   $ha    = $compara_db->get_adaptor('Homology');
   $ma    = $compara_db->get_adaptor('GeneMember');
   $gdba  = $compara_db->get_adaptor('GenomeDB');

   warn "Can't connect to Compara database specified by $compara_dbname - check command-line arguments" if (!$mlssa || !$ha || !$ma ||!$gdba);

} else {

   $mlssa = @{$registry->get_all_adaptors(-group => "compara", -type => "MethodLinkSpeciesSet")}[0];
   $ha    = @{$registry->get_all_adaptors(-group => "compara", -type => "Homology")}[0];
   $ma    = @{$registry->get_all_adaptors(-group => "compara", -type => "GeneMember")}[0];
   $gdba  = @{$registry->get_all_adaptors(-group => "compara", -type => "GenomeDB")}[0];

   warn "Can't connect to Compara database from registry - check registry file settings" if (!$mlssa || !$ha || !$ma ||!$gdba);

}

die "Could not find MethodLinkSpeciesSet adaptor" if !$mlssa;
die "Could not find Homology adaptor" if !$ha;
die "Could not find GeneMember adaptor" if !$ma;
die "Could not find GenomeDB adaptor" if !$gdba;

my $from_ga = $registry->get_adaptor($from_species, 'core', 'Gene');

my $from_gene;
my $to_species;

# Always use production name
my $from_meta_container = $registry->get_adaptor($from_species, 'core', 'MetaContainer');
my ($from_production_species) = @{ $from_meta_container->list_value_by_key('species.production_name') };

foreach my $local_to_species (@to_multi) {

  $to_species = $local_to_species;
  my $to_ga   = $registry->get_adaptor($to_species, 'core', 'Gene');
  my $to_ta   = $registry->get_adaptor($to_species, 'core', 'Transcript');
  die("Can't get gene adaptor for $to_species - check database connection details; make sure meta table contains the correct species alias\n") if (!$to_ga);
  my $to_dbea = $registry->get_adaptor($to_species, 'core', 'DBEntry');

  # Always use production name
  my $to_meta_container = $registry->get_adaptor($to_species,'core','MetaContainer');
  my ($to_production_species) = @{ $to_meta_container->list_value_by_key('species.production_name')};
  my $file = "orthologs-" . $from_production_species . "-" . $to_production_species . ".new.tsv";
  open OUT, ">$file" or die "couldn't open file " . $file . " $!";
  my $datestring = localtime();
  print OUT "## " . $datestring . "\n";
  print OUT "## orthologs from $from_production_species to $to_production_species\n";
  print OUT "## compara db " . $ma->dbc->dbname() . "\n";

  # build Compara GenomeDB objects
  my $from_GenomeDB = $gdba->fetch_by_registry_name($from_species);
  my $to_GenomeDB = $gdba->fetch_by_registry_name($to_species);

  my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$from_GenomeDB, $to_GenomeDB]);

  my $from_species_alias = $from_GenomeDB->name();

  # get homologies from compara - comes back as a hash of arrays
  my $homologies = $ha->fetch_all_by_MethodLinkSpeciesSet($mlss);

  foreach my $homology (@{$homologies}) {

    my $from_member = $homology->get_Member_by_GenomeDB($from_GenomeDB)->[0];
    my $to_members  = $homology->get_Member_by_GenomeDB($to_GenomeDB); 

    my $from_stable_id = $from_member->stable_id();
    my $from_perc_id   = $from_member->perc_id();

#    my $members = $homology->get_all_GeneMembers();
#    my @to_stable_ids;
#    my $from_stable_id;
#    foreach my $member (@{$members}) {
#      if ($member->genome_db()->name() eq $from_species_alias) {
#        $from_stable_id = $member->stable_id();
#      }
#      else {
#        push(@to_stable_ids, $member->stable_id());
#      }
#    }

    print "Warning: can't find stable ID corresponding to 'from' species ($from_species_alias)\n" if (!$from_stable_id);

    my $from_translation = $from_member->get_Transcript->translation();
#    my $from_translation = get_canonical_translation($from_stable_id, $from_ga);
    my $from_uniprot;
    if ($from_translation) { $from_uniprot = get_uniprot($from_translation); }

#    foreach my $to_stable_id (@to_stable_ids) {
     foreach my $to_member (@$to_members) {
      my $to_stable_id = $to_member->stable_id();
      my $to_perc_id   = $to_member->perc_id();
      my $to_translation  = $to_member->get_Transcript->translation();
#      my $to_translation = get_canonical_translation($to_stable_id, $to_ga);
      if (!$to_translation || !$from_translation) {
        next;
      }
      my $to_uniprot = get_uniprot($to_translation);
      if (scalar(@$from_uniprot) == 0 && scalar(@$to_uniprot) == 0) {
        print OUT $from_production_species . "\t" . $from_stable_id . "\t" . $from_translation->stable_id . "\t" . "no_uniprot" . "\t" . $from_perc_id . "\t";
        print OUT $to_production_species . "\t" . $to_stable_id . "\t" . $to_translation->stable_id . "\t" . "no_uniprot" . "\t" . $to_perc_id . "\t" . $homology->description . "\n";
      } elsif (scalar(@$from_uniprot) == 0) {
        foreach my $to_xref (@$to_uniprot) {
          print OUT $from_production_species . "\t" . $from_stable_id . "\t" . $from_translation->stable_id . "\t" . "no_uniprot" . "\t" . $from_perc_id . "\t";
          print OUT $to_production_species . "\t" . $to_stable_id . "\t" . $to_translation->stable_id . "\t" . $to_xref . "\t" . $to_perc_id. "\t" . $homology->description . "\n";
        }
      } elsif (scalar(@$to_uniprot) == 0) {
        foreach my $from_xref (@$from_uniprot) {
          print OUT $from_production_species . "\t" . $from_stable_id . "\t" . $from_translation->stable_id . "\t" . $from_xref . "\t" . $from_perc_id . "\t";
          print OUT $to_production_species . "\t" . $to_stable_id . "\t" . $to_translation->stable_id . "\t" . "no_uniprot" . "\t" . $to_perc_id . "\t" . $homology->description . "\n";
        }
      }
      foreach my $to_xref (@$to_uniprot) {
        foreach my $from_xref (@$from_uniprot) {
          print OUT $from_production_species . "\t" . $from_stable_id . "\t" . $from_translation->stable_id . "\t" . $from_xref . "\t" . $from_perc_id . "\t";
          print OUT $to_production_species . "\t" . $to_stable_id . "\t" . $to_translation->stable_id . "\t" . $to_xref . "\t" . $to_perc_id . "\t" . $homology->description . "\n";
        }
      }
    }
  }
  close OUT;
}


# ----------------------------------------------------------------------
# Get the uniprot entries associated with the canonical translation

sub get_uniprot {
  my $translation = shift;
  my $uniprots = $translation->get_all_DBEntries('Uniprot%');
  my @uniprots;
  foreach my $uniprot (@$uniprots) {
    push @uniprots, $uniprot->primary_id();
  }
  return \@uniprots;
}


# ----------------------------------------------------------------------
# Get the translation associated with the gene's canonical transcript

sub get_canonical_translation {

  my $gene_stable_id = shift;
  my $ga = shift;
  my $gene = $ga->fetch_by_stable_id($gene_stable_id);

  my $canonical_transcript = $gene->canonical_transcript();

  if (!$canonical_transcript) {
    warn("Can't get canonical transcript for " . $gene->stable_id() . ", skipping this homology");
    return undef;
  }

  return $canonical_transcript->translation();;

}

# ----------------------------------------------------------------------



sub usage {

  print << "EOF";

  Sets display_xref_ids and/or GO terms for novel genes in the "to" database
  based on their orthologs in the "from" database. Orthology relationships
  are read from a Compara database.

 perl project_display_xrefs.pl {options}

 Options ([..] indicates optional):

  [--conf filepath]     the Bio::EnsEMBL::Registry configuration file. If none
                        is given, the one set in ENSEMBL_REGISTRY will be used
                        if defined, if not ~/.ensembl_init will be used.
                        Note only the Compara database needs to be defined here,
                        assuming the rest of the databases are on the server
                        defined by --registryconf
   

   --registryconf       There are two ways in which the registry configuration
                        information can be passed to the script. This information
		        is a hash that encodes the registry configuration parameters
			and can be passed as a string in a file or as a string on the 
			commandline. 


                        Note that a combination of the host/user and conf files
                        can be used. Databases specified in both will use the
                        conf file setting preferentially.

   --from string        The species to use as the source
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

   --to string          The target species.
                        (a Bio::EnsEMBL::Registry alias, defined in config file)
                        More than one target species can be specified by using
                        several --to arguments or a comma-separated list, e.g.
                        -from human -to dog -to opossum or --to dog,opossum

   --release            The current Ensembl release. Needed for projection_info
                        database.

  [--compara string]    A Compara database
                        (a Bio::EnsEMBL::Registry alias, defined in config file)
                        If not specified, the first compara database defined in
                        the registry file is used.


EOF

}
