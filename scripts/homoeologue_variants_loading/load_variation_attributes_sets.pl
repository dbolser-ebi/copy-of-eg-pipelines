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
# perl load_variation_attributes_sets.pl -synonyms wheat10.synonyms -host [db host] -port [db port] -user [db user] -pass [db pass] -dbname [db name]
#
# Load variation_set and variation_attrib tables from a file like wheat10.synonyms
# wheat10.synonyms e.g.
# 564321	A	EPlTAEV01493686	A	D	EPlTAEV01493687	G
# 544357	A	EPlTAEV01445614	T	D	EPlTAEV01445616	G	B	EPlTAEV01445615	G
# 81834	A	EPlTAEV02101942	G	B	EPlTAEV02101943	A
#

use warnings;
use strict;
use FileHandle;
use DBI;
use Getopt::Long;

my ($filename, $dbname, $host, $port, $user, $pass);

GetOptions( "host=s",       \$host,
            "user=s",       \$user,
            "pass=s",       \$pass,
            "port=i",       \$port,
            "dbname=s",     \$dbname,
	    "synonyms=s",   \$filename,
          );

my $verbose = 0;

if (! defined $filename) {
  die "specify a synonyms file name\n";
}
    
my $dsn = "DBI:mysql:database=$dbname;host=$host;port=$port";
my $dbh = DBI->connect("$dsn",$user,$pass, {RaiseError => 1});
 
print STDERR "Loading mapping from $filename...\n";

# Key is the component genome
# Value is the stable id

my $attributes_href = {};

# Parse the group mapping file

my $synonyms_fh = new FileHandle;
$synonyms_fh->open("<$filename") or die "can't open synonyms file, $filename!\n";
while (<$synonyms_fh>) {
  my $line = $_;
  chomp $line;
  my @cols = split("\t",$line);

  my $attributes_per_locus_href = {};
  
  for (my $i=1; $i<@cols; $i+=3) {
      
      # Key is the component genome
      # Value is a hash ref with the stable id and the ref. allele

      $attributes_per_locus_href->{$cols[$i]} = {
          'stable_id' => $cols[$i+1],
          'allele' => $cols[$i+2],
      }
  }

  # Key is the group_id
  # Value is the mapping data for a given locus
  $attributes_href->{$cols[0]} = $attributes_per_locus_href;
  
}
$synonyms_fh->close();

print STDERR "Parsing mapping done\n";

# Todo: Make sure these IDs are right

my $A_attrib_type_id = 420;
my $B_attrib_type_id = 421;
my $D_attrib_type_id = 422;

my $vs_AB = 4;
my $vs_AD = 5;
my $vs_BD = 6;
my $vs_ABD = 7;
my $vs_ABD_D_differ = 8;
my $vs_ABD_B_differ = 9;
my $vs_ABD_A_differ = 10;

print STDERR "Loading variation_attrib and variation_set_variation mapping...\n";

my $get_variation_id_sql = "SELECT variation_id FROM variation WHERE name = ?";
my $get_variation_id_sth = $dbh->prepare($get_variation_id_sql);

my $insert_va_sql = "INSERT INTO variation_attrib SELECT variation_id, ?,? FROM variation WHERE name = ?";
my $insert_va_sth = $dbh->prepare($insert_va_sql);

my $insert_vsv_sql = "INSERT INTO variation_set_variation VALUES (?,?)";
my $insert_vsv_sth = $dbh->prepare($insert_vsv_sql);

foreach my $group_id (keys (%$attributes_href)) {
    
    my $attributes_per_locus_href = $attributes_href->{$group_id};

    my @component_genomes = keys (%$attributes_per_locus_href);
    
    # Associate variations and variation_attribs

    foreach my $component_genome (@component_genomes) {
	my $attribute_href = $attributes_per_locus_href->{$component_genome};
	my $stable_id = $attribute_href->{'stable_id'};
	
	foreach my $xref_component_genome (keys (%$attributes_per_locus_href)) {
	    if ($xref_component_genome ne $component_genome) {

		my $attrib_type_id = undef;
		if ($xref_component_genome eq "A") {
		    $attrib_type_id = $A_attrib_type_id;
		}
		elsif ($xref_component_genome eq "B") {
		    $attrib_type_id = $B_attrib_type_id;
		}
		elsif ($xref_component_genome eq "D") {
		    $attrib_type_id = $D_attrib_type_id;
		}
		else {
		    warn("Houston, we have a problem!\n");
		}

		my $xref_attribute_href = $attributes_per_locus_href->{$xref_component_genome};
		my $xref_stable_id = $xref_attribute_href->{'stable_id'};
		
		if ($verbose) {
		    warn("loading attrib_type_id, xref_stable_id, $attrib_type_id, $xref_stable_id, for stable_id, $stable_id, from component, $component_genome\n");
		}
		
		$insert_va_sth->execute($attrib_type_id, $xref_stable_id, $stable_id);
	    }
	}
	
    }

    # Associate variations and variation_sets

    my $vs_id = undef;

    if (@component_genomes == 3) {
        # also need to check the component alleles, so let's process the allele info for that later on
        
        # Get which variation set it falls in
        
	if ($verbose) {
	    warn("locus present in the 3 components, finding out which variation_set it falls in\n");
	}
	
        $vs_id = get_vs_id ($attributes_per_locus_href);
             
    }
    elsif (($component_genomes[0] =~ /A|B/) && ($component_genomes[1] =~ /A|B/)) {
	$vs_id = $vs_AB;
    }
    elsif (($component_genomes[0] =~ /A|D/) && ($component_genomes[1] =~ /A|D/)) {
	$vs_id = $vs_AD;
    }
    elsif (($component_genomes[0] =~ /B|D/) && ($component_genomes[1] =~ /B|D/)) {
	$vs_id = $vs_BD;
    }
    else {
	warn("what am I doing here !!!\n");
    }
    
    if (! defined $vs_id) {
	
	# Skip it
        # This is in the case where the variant alleles are identical on the forward strand
	# don't know what to do here yet
	
	next;
    }

    foreach my $component_genome (@component_genomes) {
	
	my $attribute_href = $attributes_per_locus_href->{$component_genome};
	my $stable_id = $attribute_href->{'stable_id'};
	
	$get_variation_id_sth->execute($stable_id);
	my ($variation_id) = $get_variation_id_sth->fetchrow_array();

	if ($verbose) {
	    print STDERR "inserting vvs row with variation_id, variation_set_id, $variation_id, $vs_id\n";
	}

	$insert_vsv_sth->execute($variation_id,$vs_id);
    }

}

$get_variation_id_sth->finish();

$insert_va_sth->finish();
$insert_vsv_sth->finish();

# End

sub get_vs_id {
    my ($attributes_per_locus_href) = @_;
    
    my $vs_id = undef;
    
    my $A_allele = $attributes_per_locus_href->{'A'}->{'allele'};
    my $B_allele = $attributes_per_locus_href->{'B'}->{'allele'};
    my $D_allele = $attributes_per_locus_href->{'D'}->{'allele'};
    
    if ($A_allele eq $B_allele) {
        if ($A_allele eq $D_allele) {
	    
            # Can not happen, means it is not a variant !!
	    # Well probably it can happen when the alt_allele is on the reverse strand

	    warn ("Allele A, B and D, $A_allele, $B_allele and $D_allele\n");
	    warn ("That's not a variant then!\n");
	    
        }
        else {
            $vs_id = $vs_ABD_D_differ;
        }
    }
    else {
        if ($A_allele eq $D_allele) {
            $vs_id = $vs_ABD_B_differ;
        }
        elsif ($B_allele eq $D_allele) {
	    $vs_id = $vs_ABD_A_differ;
	}
	else {
            $vs_id = $vs_ABD;
        }
    }
    
    return $vs_id;
    
}
