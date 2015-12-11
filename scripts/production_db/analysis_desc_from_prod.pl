#!/usr/bin/env/perl
use strict;
use warnings;

# Update a new database with analysis descriptions from the production
# database.

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use DBI;

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
die "--pass required" unless $pass;
die "--dbname required" unless $dbname;
die "--mhost required" unless $mhost;
die "--mport required" unless $mport;
die "--muser required" unless $muser;
die "--species and --type are both required" if ( ($species && !$type) || (!$species && $type) );

$mdbname = "ensembl_production" unless $mdbname;
($species, $type) = $dbname =~ /^([^_]+_[^_]+)_([^_]+)/ unless $species && $type;
die "Could not determine species and type from $dbname" unless $species && $type;

my %logic_names = map { $_ => 1 } @logic_names;

my $new_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
  -host   => $host,
  -port   => $port,
  -user   => $user,
  -pass   => $pass,
  -dbname => $dbname,
);

# Fetch all data from the master database (lifted from 'populate_analysis_description.pl).
my %mdata;
{
  my ($logic_name, %properties);
  my $dsn = "DBI:mysql:host=$mhost;port=$mport;database=$mdbname";
  my $dbh = DBI->connect($dsn, $muser, $mpass, { 'PrintError'=>1, 'RaiseError'=>1 });

  # Load generic analyses
  my $sth = $dbh->prepare(
    'SELECT ad.logic_name, ad.description, ad.display_label, wd.data, 1 '.
    'FROM analysis_description ad '.
    'LEFT OUTER JOIN web_data wd ON ad.default_web_data_id = wd.web_data_id '.
    'WHERE ad.is_current = 1;'
  );
  $sth->execute();

  $sth->bind_columns(\(
    $logic_name,
    $properties{'description'},
    $properties{'display_label'},
    $properties{'web_data'},
    $properties{'displayable'}
  ));

  while ( $sth->fetch() ) {
    $mdata{$logic_name} = { %properties };
  }

  # Load species-specific analyses, overwriting the generic info
  $sth = $dbh->prepare(
    'SELECT ad.logic_name, ad.description, ad.display_label, wd.data, aw.displayable '.
    'FROM analysis_description ad, species s, analysis_web_data aw '.
    'LEFT OUTER JOIN web_data wd ON aw.web_data_id = wd.web_data_id '.
    'WHERE ad.analysis_description_id = aw.analysis_description_id '.
    'AND aw.species_id = s.species_id '.
    'AND s.db_name = ? '.
    'AND aw.db_type =? '
  );
  $sth->execute($species, $type);

  $sth->bind_columns(\(
    $logic_name,
    $properties{'description'},
    $properties{'display_label'},
    $properties{'web_data'},
    $properties{'displayable'}
  ));

  while ( $sth->fetch() ) {
    $mdata{$logic_name} = { %properties };
  }

  $dbh->disconnect();
}

my $new_analysis_adaptor = $new_dba->get_adaptor('Analysis');

foreach my $new_analysis ( @{ $new_analysis_adaptor->fetch_all() } ) {
  my $logic_name = $new_analysis->logic_name();
  if (scalar(@logic_names) == 0 || exists $logic_names{$logic_name}) {
    if (exists $mdata{$logic_name}) {
      print "Updating $logic_name from $mdbname...";
      $mdata{$logic_name}{'web_data'} = eval ($mdata{$logic_name}{'web_data'}) if defined $mdata{$logic_name}{'web_data'};
      $new_analysis->description($mdata{$logic_name}{'description'});
      $new_analysis->display_label($mdata{$logic_name}{'display_label'});
      $new_analysis->web_data($mdata{$logic_name}{'web_data'});
      $new_analysis->displayable($mdata{$logic_name}{'displayable'});
      if ($new_analysis_adaptor->update($new_analysis)) {
        print "Done\n";
      } else {
        print "Failed\n";
      }
    } else {
      print "No data for $logic_name.\n";
    }
  }
}
