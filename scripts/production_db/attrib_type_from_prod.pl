#!/usr/bin/env/perl
use strict;
use warnings;

# Update a database with attrib_type data from the production database.

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use DBI;

my ($host, $port, $user, $pass, $dbname,
    $mhost, $mport, $muser, $mpass, $mdbname,
    $no_insert, $no_update, $no_delete, $no_backup);

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
  "no_insert", \$no_insert,
  "no_update", \$no_update,
  "no_delete", \$no_delete,
  "no_backup", \$no_backup,
);

die "--host required" unless $host;
die "--port required" unless $port;
die "--user required" unless $user;
die "--pass required" unless $pass;
die "--dbname required" unless $dbname;
die "--mhost required" unless $mhost;
die "--mport required" unless $mport;
die "--muser required" unless $muser;
$mdbname = "ensembl_production" unless $mdbname;
$no_insert = 0 unless $no_insert;
$no_update = 0 unless $no_update;
$no_delete = 0 unless $no_delete;
$no_backup = 0 unless $no_backup;

my @cols = qw(attrib_type_id code name description);

my $backup_sql = 'CREATE TABLE attrib_type_bak AS SELECT * FROM attrib_type;';
my $select_sql = 'SELECT '.join(', ', @cols).' FROM attrib_type;';
my $insert_sql = 'INSERT IGNORE INTO attrib_type ('.join(', ', @cols).') VALUES ('.join(', ', (map {"?"} @cols)).');';
my $update_sql = 'UPDATE attrib_type SET '.join(', ', (map {"$_ = ?"} @cols)).' WHERE attrib_type_id = ?;';
my $delete_sql = 'DELETE FROM attrib_type WHERE attrib_type_id = ?;';

my %mdata;
my $mdsn = "DBI:mysql:host=$mhost;port=$mport;database=$mdbname";
my $mdbh = DBI->connect($mdsn, $muser, $mpass, { 'PrintError'=>1, 'RaiseError'=>1 });
my $msth = $mdbh->prepare($select_sql);
$msth->execute();
%mdata = %{$msth->fetchall_hashref('attrib_type_id')};
$mdbh->disconnect();

my %data;
my $dsn = "DBI:mysql:host=$host;port=$port;database=$dbname";
my $dbh = DBI->connect($dsn, $user, $pass, { 'PrintError'=>1, 'RaiseError'=>1 });
my $sth = $dbh->prepare($select_sql);
$sth->execute();
%data = %{$sth->fetchall_hashref('attrib_type_id')};

# Backup existing table
if (!$no_backup) {
  $sth = $dbh->prepare($backup_sql);
  $sth->execute();
}

# Add new data
if (!$no_insert) {
  $sth = $dbh->prepare($insert_sql);
  foreach my $attrib_type_id (keys %mdata) {
    if (!exists $data{$attrib_type_id}) {
      $sth->execute(map {$mdata{$attrib_type_id}{$_}} @cols);
      print "Inserted data for attrib_type_id $attrib_type_id (".$mdata{$attrib_type_id}{'code'}.")\n";
    }
  }
}

# Update data
if (!$no_update) {
  $sth = $dbh->prepare($update_sql);
  foreach my $attrib_type_id (keys %mdata) {
    if (exists $data{$attrib_type_id}) {
      if (join('', (map {$mdata{$attrib_type_id}{$_} || ''} @cols)) ne 
          join('', (map {$data{$attrib_type_id}{$_} || ''} @cols))
      ) {
        $sth->execute((map {$mdata{$attrib_type_id}{$_}} @cols), $attrib_type_id);
        print "Updated data for attrib_type_id $attrib_type_id (".$mdata{$attrib_type_id}{'code'}.")\n";
      }
    }
  }
}

# Delete data
if (!$no_delete) {
  $sth = $dbh->prepare($delete_sql);
  foreach my $attrib_type_id (keys %data) {
    if (!exists $mdata{$attrib_type_id}) {
      $sth->execute($attrib_type_id);
      print "Deleted data for attrib_type_id $attrib_type_id (".$data{$attrib_type_id}{'code'}.")\n";
    }
  }
}
