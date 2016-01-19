#!/usr/bin/env/perl
use strict;
use warnings;

# Add a new species to the production database.

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use DBI;

my ($host, $port, $user, $pass, $dbname,
    $mhost, $mport, $muser, $mpass, $mdbname,
    $add_to_division);

GetOptions(
  "host=s", \$host,
  "P|port=i", \$port,
  "user=s", \$user,
  "p|pass=s", \$pass,
  "dbname=s", \$dbname,
  "mhost=s", \$mhost,
  "mP|mport=i", \$mport,
  "muser=s", \$muser,
  "mp|mpass=s", \$mpass,
  "mdbname=s", \$mdbname,
  "add_to_division", \$add_to_division,
);

die "--host required" unless $host;
die "--port required" unless $port;
die "--user required" unless $user;
die "--dbname required" unless $dbname;
die "--mhost required" unless $mhost;
die "--mport required" unless $mport;
die "--muser required" unless $muser;
die "--mpass required" unless $mpass;

$add_to_division = 1 unless defined $add_to_division;

$mdbname = "ensembl_production" unless $mdbname;

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
  -host   => $host,
  -port   => $port,
  -user   => $user,
  -pass   => $pass,
  -dbname => $dbname,
);

my $dsn = "DBI:mysql:host=$mhost;port=$mport;database=$mdbname";
my $dbh = DBI->connect($dsn, $muser, $mpass, { 'PrintError'=>1, 'RaiseError'=>1 });
my $sth = $dbh->prepare(
  'INSERT INTO species ('.
    'db_name, '.
    'common_name, '.
    'web_name, '.
    'is_current, '.
    'taxon, '.
    'scientific_name, '.
    'production_name, '.
    'url_name) '.
  'VALUES (?, ?, ?, 1, ?, ?, ?, ?);'
);

my $meta_container = $dba->get_MetaContainer();
my $prod_name = $meta_container->single_value_by_key('species.production_name');
print "Adding $prod_name to $mdbname...";

my $return = $sth->execute(
  $prod_name,
  $meta_container->single_value_by_key('species.common_name') || $meta_container->single_value_by_key('species.production_name'),
  $meta_container->single_value_by_key('species.scientific_name'),
  $meta_container->single_value_by_key('species.taxonomy_id'),
  $meta_container->single_value_by_key('species.scientific_name'),
  $prod_name,
  $meta_container->single_value_by_key('species.url')
);

if ($return) {
  print "Done\n";
} else {
  print "Failed\n";
}

if ($add_to_division) {
  my $div_sth = $dbh->prepare(
    'INSERT INTO division_species (division_id, species_id) '.
    'SELECT division.division_id, species.species_id FROM division, species '.
    'WHERE division.name = ? AND species.production_name = ?;'
  );
  
  my $div_name = $meta_container->single_value_by_key('species.division');
  if ($div_name) {
    print "Adding $prod_name to $div_name division...";
  
    my $div_return = $div_sth->execute($div_name, $prod_name);
    
    if ($div_return) {
      print "Done\n";
    } else {
      print "Failed\n";
    }
  }
}

