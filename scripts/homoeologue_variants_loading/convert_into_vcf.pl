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
# e.g. perl convert_into_vcf.pl -variants wheat10.txt -groups groups.txt
#
# Generate 4 vcf files, and a stable_id mapping file:
# * Three vcf files for each pair combination
# * A vcf file where a SNP is shared by the three component genomes (which is at the moment in one of the 3 vcf files above)
# * A variant stable ids mapping file
#

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

my ($variants_text_file, $groups_text_file);

my $variants_text_fh = new FileHandle;
my $groups_text_fh = new FileHandle;

GetOptions( 
	    "variants=s",   \$variants_text_file,
            "groups=s",     \$groups_text_file,
          );

my $groups_href = {};

print STDERR "Processing groups file...\n";

$groups_text_fh->open("<$groups_text_file") or die "can't open file, $groups_text_file!\n";
while (<$groups_text_fh>) {
    my $line = $_;
    chomp $line;
    my @cols = split ("\t",$line);

    my $group_id = shift @cols;

    $groups_href->{$group_id} = \@cols;
    
}
$groups_text_fh->close();

my $index = 1;
my $prefix = "EPlTAEV";

my $group_index_href = {};

# Keys are 'AB', 'AD', 'BD' and 'ABD'

my $variants_per_genome_set_href = {};

print STDERR "Processing text variants file...\n";

$variants_text_fh->open("<$variants_text_file") or die "can't open file, $variants_text_file!\n";
while (<$variants_text_fh>) {
    my $line = $_;
    chomp $line;
    my @cols = split (',',$line);
    
    my $start = $cols[3];
    my $end = $cols[4];
    if ($start > $end) {
	$start = $end;
    }

    my $group_index = $cols[0];

    if (! defined $group_index) {
	die "group index not set in line, $line\n";
    }
    elsif ($group_index eq "") {
        warn ("group index eq ''!\n");
    }

    my $ref_allele = $cols[10] || $cols[8];
    my $alt_allele = $cols[11] || $cols[9];
   
    #
    # '-' is wrong actually, this is against vcf specs
    # Fix by Paul so this should not happen anymore anyway
    #
    if ($ref_allele eq '-') {
	$ref_allele = '.';
    }
    if (!defined $alt_allele || $alt_allele eq '-') {
	$alt_allele = '.';
    }
    #
    
    #print STDERR "ref_allele, alt_allele, $ref_allele, $alt_allele\n";

    my $ref_component_genome = $cols[1];

    # Set a stable_id
    
    my $stable_id = sprintf "%s%08d", $prefix , $index+1;   
    
    # Add group_index info to group_index_href hash so we can map SNP loci between genome components

    my $group_aref = $groups_href->{$group_index};

    if (! defined $group_index_href->{$group_index}) {

	# new group

	my $mapping_href = {
	    $ref_component_genome => {
		'stable_id' => $stable_id,
		'ref_allele' => $ref_allele,
		'alt_allele' => $alt_allele,
		'start' => $start,
	    }
	};
	$group_index_href->{$group_index} = $mapping_href;
    }
    else {

	# existing group

	my $mapping_href = $group_index_href->{$group_index};
	if (!defined $mapping_href->{$ref_component_genome}) {

	    # new component genome

	    $mapping_href->{$ref_component_genome} = {
		'stable_id' => $stable_id,
		'ref_allele' => $ref_allele,
		'alt_allele' => $alt_allele,
		'start' => $start,
	    };
	}
	else {

	    # existing component genome, so filter out and merge the alt allele

	    if ($start != $mapping_href->{$ref_component_genome}->{'start'}) {
		# Bug in the produced data
		# Somehow we have sometimes more than one variants loci (in another component genome) for a given homoeologue
		# should not happen as the genes are one to one homoeologues only
		# so filter the second (arbitral decision)

		next;
		
	    }
	    else {
		# We have already assigned a stable id, so let's use it then
		
		# This is in the case a variant is present is the three component genomes
		# Add the new alt allele to the list of alt_alleles
		
		$stable_id = $mapping_href->{$ref_component_genome}->{'stable_id'};
		
		if (@$group_aref == 3) {
		    
		    # WHEAT ABD variant, so we expect duplication in there

		    if ($alt_allele ne $mapping_href->{$ref_component_genome}->{'alt_allele'}) {

			# Todo: in these cases, somehow, we still need to remove the first occurence then
			# which has been already printed out, so we'll have to change a lot of the logic in there

			$alt_allele .= "," . $mapping_href->{$ref_component_genome}->{'alt_allele'};
		    }
		    else {
			# then we don't need to print it as in it an exact duplication
			
			next;
		    }
		}
		else {
		    # non WHEAT ABD, so we don't expect any duplication
		    # so skip the second occurrence then
		    next;
		}
	    }
	}
    }
    
    my $group_key = undef;

    if (@$group_aref == 3) {
	$group_key = "ABD";

    }
    elsif ($group_aref->[0] eq "A" && $group_aref->[1] eq "B") {
	$group_key = "AB";
    }
    elsif ($group_aref->[0] eq "A" && $group_aref->[1] eq "D") {
	$group_key = "AD";
    }
    else {
	$group_key = "BD";
    }

    if (!defined $variants_per_genome_set_href->{$group_key}) {
	$variants_per_genome_set_href->{$group_key} = {};
    }
    my $variants_per_group_href = $variants_per_genome_set_href->{$group_key};
    if (! defined $variants_per_group_href->{$group_index}) {
	$variants_per_group_href->{$group_index} = {};
    }
    my $variants_per_genome_href = $variants_per_group_href->{$group_index};
    if (! defined $variants_per_genome_href->{$ref_component_genome}) {
	$variants_per_genome_href->{$ref_component_genome} = [];
    }
    else {
	# get rid of the first element as it holds the wrong alt_allele info
	# it should contains one element
	$variants_per_genome_href->{$ref_component_genome} = [];
    }
    # an array so we keep them in order!
    my $variants_aref = $variants_per_genome_href->{$ref_component_genome};
    my $variant_href = {
	'seq_name' => $cols[2],
	'start' => $start,
	'variant_name' => $stable_id,
	'ref_allele' => $ref_allele,
	'alt_allele' => $alt_allele,
    };

    push(@$variants_aref, $variant_href);
    
    # Let's output later as for wheat_ABD we need to keep one instance of the variant loci (instead of the two)
    #print $out_fh "$cols[2]\t$start\t$stable_id\t$ref_allele\t$alt_allele\t.\t.\t.\n";
    
    $index++;
}

$variants_text_fh->close();

print STDERR "Processing text variants file done\n";

print STDERR "Mapping variants between component genomes...\n";

# For each group, report the list of variant stable ids

my $syn_fh = new FileHandle;
$syn_fh->open (">wheat10.synonyms") or die "can't open synonym output file!\n";

foreach my $group_index (keys (%$group_index_href)) {
    my $stable_ids_href = $group_index_href->{$group_index};
    my @component_genomes_aref = keys ($stable_ids_href);
    
    if (@component_genomes_aref > 1) {

	print $syn_fh "$group_index";
	
	foreach my $component_genome (@component_genomes_aref) {
	    print $syn_fh "\t$component_genome\t" . $stable_ids_href->{$component_genome}->{'stable_id'} . "\t" . $stable_ids_href->{$component_genome}->{'ref_allele'};
	}
	print $syn_fh "\n";
    }
    else {
	warn("This situation can not possibly happen!!\n");
    }

}

$syn_fh->close();

print STDERR "Mapping variants done\n";

print STDERR "now let's print out the vcf files...\n";

my $wheat_AB_fh = new FileHandle;
$wheat_AB_fh->open(">wheat_AB.vcf") or die "can't open wheat_AB vcf file for writing!\n";
my $wheat_AD_fh = new FileHandle;
$wheat_AD_fh->open(">wheat_AD.vcf") or die "can't open wheat_AD vcf file for writing!\n";
my $wheat_BD_fh = new FileHandle;
$wheat_BD_fh->open(">wheat_BD.vcf") or die "can't open wheat_BD vcf file for writing!\n";
my $wheat_ABD_fh = new FileHandle;
$wheat_ABD_fh->open(">wheat_ABD.vcf") or die "can't open wheat_ABD vcf file for writing!\n";

# Add the header lines to the VCF files

# Let's go through the variants now

foreach my $genome_set_key (keys (%$variants_per_genome_set_href)) {
    my $out_fh = new FileHandle;
    if ($genome_set_key eq "AB") {
	$out_fh->open(">wheat_AB.vcf") or die "can't open wheat_AB vcf file for writing!\n";
    }
    elsif ($genome_set_key eq "AD") {
	$out_fh->open(">wheat_AD.vcf") or die "can't open wheat_AD vcf file for writing!\n";
    }
    elsif ($genome_set_key eq "BD") {
	$out_fh->open(">wheat_BD.vcf") or die "can't open wheat_BD vcf file for writing!\n";
    }
    elsif ($genome_set_key eq "ABD") {
	$out_fh->open(">wheat_ABD.vcf") or die "can't open wheat_ABD vcf file for writing!\n";
    }
    else {
	die "genome_set, $genome_set_key, is unknown!\n";
    }
    
    print $out_fh "##fileformat=VCFv4.2\n##fileDate=20140728\n##reference=ftp://ftp.ensemblgenomes.org/pub/plants/release-23/fasta/triticum_aestivum/dna/\n";
    print $out_fh '##INFO=<ID=NA,Number=1,Type=String,Description="Unknown">' . "\n";
    print $out_fh '##FILTER=<ID=NA,Description="Unknown">' . "\n";
    print $out_fh '##FORMAT=<ID=NA,Number=1,Type=String,Description="Unknown">'. "\n";
    print $out_fh "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

    my $variants_per_group_href = $variants_per_genome_set_href->{$genome_set_key};
    foreach my $group_id_key (keys (%$variants_per_group_href)) {
	my $variants_per_genome_href = $variants_per_group_href->{$group_id_key};
	foreach my $genome_component_name (keys (%$variants_per_genome_href)) {
	    my $variants_aref = $variants_per_genome_href->{$genome_component_name};
	    my $variant_href = $variants_aref->[0];

	    print $out_fh $variant_href->{'seq_name'} . "\t" . $variant_href->{'start'} . "\t" . $variant_href->{'variant_name'} . "\t" . $variant_href->{'ref_allele'} . "\t" . $variant_href->{'alt_allele'} . "\t.\t.\t.\n";
	}
    } 

    $out_fh->close();
}

