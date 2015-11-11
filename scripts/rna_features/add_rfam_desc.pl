#!/usr/bin/env/perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Path::Tiny qw(path);

# This script adds descriptions to the Rfam.cm file, so that they
# appear in the results when you use the file.
# The cm_file should be downloaded and unzipped from: 
#  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
# The family_file should be downloaded and unzipped from:
#  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
#
# Once you've got a new Rfam.cm file, save it (as ensgen) to:
#  /nfs/panda/ensemblgenomes/external/Rfam/<rfam_release>/Rfam.cm
# and generate the indices with:
#  cd /nfs/panda/ensemblgenomes/external/Rfam/<rfam_release>
#  /nfs/panda/ensemblgenomes/external/bin/cmpress Rfam.cm

my ($cm_file, $family_file);

GetOptions(
  "cm_file=s", \$cm_file,
  "family_file=s", \$family_file,
);

die '-cm_file is required and must exist' unless $cm_file && -e $cm_file;
die '-family_file is required and must exist' unless $family_file && -e $family_file;

my $family_path = path($family_file);
my $families = $family_path->slurp;
my %families = $families =~ /^([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)/gm;

my $cm_path = path($cm_file);
my $cm = $cm_path->slurp;
$cm =~ s/^(ACC\s+)(\S+)(\nSTATES\s+)/$1$2\nDESC     $families{$2}$3/gm;
$cm =~ s/^(ACC\s+)(\S+)(\nLENG\s+)/$1$2\nDESC  $families{$2}$3/gm;
$cm_path->spew($cm);
