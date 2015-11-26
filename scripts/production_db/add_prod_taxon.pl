#!/usr/bin/env/perl
use strict;
use warnings;

# Add a new species to the production database.

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use DBI;

my ($host, $port, $user, $pass, $dbname,
    $mhost, $mport, $muser, $mpass, $mdbname);

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
);

die "--host required" unless $host;
die "--port required" unless $port;
die "--user required" unless $user;
die "--dbname required" unless $dbname;
die "--mhost required" unless $mhost;
die "--mport required" unless $mport;
die "--muser required" unless $muser;
die "--mpass required" unless $mpass;

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
print "Adding ".$meta_container->single_value_by_key('species.production_name')." to $mdbname...";

my $return = $sth->execute(
  $meta_container->single_value_by_key('species.production_name'),
  $meta_container->single_value_by_key('species.common_name') || $meta_container->single_value_by_key('species.production_name'),
  $meta_container->single_value_by_key('species.scientific_name'),
  $meta_container->single_value_by_key('species.taxonomy_id'),
  $meta_container->single_value_by_key('species.scientific_name'),
  $meta_container->single_value_by_key('species.production_name'),
  $meta_container->single_value_by_key('species.url')
);

if ($return) {
  print "Done\n";
} else {
  print "Failed\n";
}
