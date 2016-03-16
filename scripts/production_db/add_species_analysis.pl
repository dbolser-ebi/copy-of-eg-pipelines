#!/usr/bin/env/perl
use strict;
use warnings;

# Add links between species and analyses.

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($host, $port, $user, $pass, $dbname,
    $mhost, $mport, $muser, $mpass, $mdbname,
    $species, $type, @logic_names);

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
  "species:s", \$species,
  "type:s", \$type,
  "logic_name:s", \@logic_names,
);

die "--host required" unless $host;
die "--port required" unless $port;
die "--user required" unless $user;
$pass = "" unless $pass;
die "--dbname required" unless $dbname;
die "--mhost required" unless $mhost;
die "--mport required" unless $mport;
die "--muser required" unless $muser;
die "--mpass required" unless $mpass;
die "--species and --type are both required" if ( ($species && !$type) || (!$species && $type) );

$mdbname = "ensembl_production" unless $mdbname;
($species, $type) = $dbname =~ /^([^_]+_[^_]+)_([^_]+)/ unless $species && $type;
die "Could not determine species and type from $dbname" unless $species && $type;

my %logic_names = map { $_ => 1 } @logic_names;

my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
  -host   => $host,
  -port   => $port,
  -user   => $user,
  -pass   => $pass,
  -dbname => $dbname,
);

# Fetch all data from the master database (modelled on 'populate_analysis_description.pl).
my $species_id;
my %mdata;

my ($logic_name, %properties);
my $dsn = "DBI:mysql:host=$mhost;port=$mport;database=$mdbname";
my $dbh = DBI->connect($dsn, $muser, $mpass, { 'PrintError'=>1, 'RaiseError'=>1 });

my $sth = $dbh->prepare(
  'SELECT species_id FROM species WHERE db_name = ?;'
);
$sth->execute($species);
$sth->bind_columns(\($species_id));
$sth->fetch();

unless ($species_id) {
  print "Species '$species' doesn't exist in the production database; use add_prod_taxon.pl\n";
  exit;  
}

$sth = $dbh->prepare(
  'SELECT logic_name, analysis_description_id, default_web_data_id, default_displayable '.
  'FROM analysis_description '.
  'WHERE is_current = 1;'
);
$sth->execute();

$sth->bind_columns(\(
  $logic_name,
  $properties{'analysis_description_id'},
  $properties{'web_data_id'},
  $properties{'displayable'}
));

while ( $sth->fetch() ) {
  $mdata{$logic_name} = { %properties };
}

$sth = $dbh->prepare(
  'INSERT IGNORE INTO analysis_web_data '.
  '(analysis_description_id, web_data_id, displayable, species_id, db_type, created_at, modified_at) '.
  'VALUES (?, ?, ?, ?, ?, NOW(), NOW());'
);

my $analysis_adaptor = $core_dba->get_adaptor('Analysis');

foreach my $analysis ( @{ $analysis_adaptor->fetch_all() } ) {
  my $logic_name = $analysis->logic_name();
  if (scalar(@logic_names) == 0 || exists $logic_names{$logic_name}) {
    if (exists $mdata{$logic_name}) {
      print "Setting $logic_name and $species in $mdbname...";
      my $return = $sth->execute(
        $mdata{$logic_name}{'analysis_description_id'},
        $mdata{$logic_name}{'web_data_id'},
        $mdata{$logic_name}{'displayable'},
        $species_id,
        $type,
      );

      if ($return) {
        print "Done\n";
      } else {
        print "Failed\n";
      }
    } else {
      print "No data for $logic_name in production database.\n";
    }
  }
}
