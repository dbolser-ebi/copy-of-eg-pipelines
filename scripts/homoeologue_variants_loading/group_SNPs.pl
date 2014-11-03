#!/usr/bin/env perl

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


#
# Generate a tab-delimited file, specifying, for each group of variant loci, in which component genome it is present
# perl group_SNPs.pl -variants wheat10.txt > groups.txt
#
# e.g.
# 127219	A	D
# means that the loci group is present in 'A' and 'D', but not in 'B'
#

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

my $filename = undef;

GetOptions( "variants=s",    \$filename,
    );

my $groups_href = {};

my $in_fh = FileHandle->new();
$in_fh->open("<$filename") or die "can't open file with variants, $filename!\n";

print STDERR "Parsing $filename...\n";

while (<$in_fh>) {
    my $line = $_;
    chomp $line;
    my @cols = split (',',$line);

    my $group_id = $cols[0];
    my $ref_component_genome = $cols[1];
    my $alt_component_genome = $cols[6];

    my $group_href = {};
    if (defined $groups_href->{$group_id}) {
	$group_href = $groups_href->{$group_id};
    }
    else {
	$groups_href->{$group_id} = $group_href;
    }
    
    $group_href->{$ref_component_genome} = 1;
    $group_href->{$alt_component_genome} = 1;

}
$in_fh->close();

print STDERR "Parsing file done, reporting groups now...\n";

# Let's report now whether a SNP locus is part of two components or the whole three

foreach my $group_id (keys (%$groups_href)) {
    my $group_href = $groups_href->{$group_id};
    my @component_genomes_aref = keys (%$group_href);

    @component_genomes_aref = sort (@component_genomes_aref);

    print "$group_id\t" . join ("\t", @component_genomes_aref) . "\n";
}
