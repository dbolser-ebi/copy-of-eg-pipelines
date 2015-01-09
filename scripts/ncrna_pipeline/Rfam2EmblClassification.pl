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

use warnings;
use strict;
use Getopt::Long;

use vars qw($opt_full);

my $_debug = 0;

GetOptions("full=s");

if ((!defined $opt_full) || (! -f $opt_full)) {
    if (defined $opt_full) {
	print STDERR "Can't find Rfam.full input file, $opt_full!\n";
    }
    else {
	print STDERR "No Rfam.full input file name specified\n";
    }
    print STDERR "Usage, perl Rfam2EmblClassification.pl -full <path_to_Rfam.full>\n";
    exit 1;
}

my $cisreg_models    = [];
my $sRNA_models      = [];
my $snoRNA_models    = [];
my $antisense_models = [];
my $rfam_models      = [];

my $rfam_acc;
my $rfam_name;
my $rfam_desc;
my $rfam_type;

my $rfam_model_name_href = {};
my $rfam_model_desc_href = {};

open FULL, "<$opt_full" or die "can't find Rfam.full input file, $opt_full!\n";
while (<FULL>) {
    my $line = $_;
    
    if ($line =~ /^#=GF AC\s+(\w+)/) {

	$rfam_acc = $1;
	
	if ($_debug) {
	    # print STDERR "parsed Rfam model, $rfam_acc\n";
	}
    }

    if ($line =~ /^#=GF ID\s+(.+)/) {
	$rfam_name = $1;
	chomp $rfam_name;
	
	if ($_debug) {
	    # print STDERR "parsed Rfam id, $rfam_name\n";
	}

	$rfam_model_name_href->{$rfam_acc} = $rfam_name;

    }

    if ($line =~ /^#=GF DE\s+(.+)/) {
	$rfam_desc = $1;
	chomp $rfam_desc;
	$rfam_model_desc_href->{$rfam_acc} = $rfam_desc;

	if ($_debug) {
	    # print STDERR "parsed Rfam desc, $rfam_desc\n";
	}
    }
    
    if ($line =~ /^#=GF TP\s+(.+)/) {
	$rfam_type = $1;
	chomp $rfam_type;
	
	if ($_debug) {
	    # print STDERR "parsed Rfam type, $rfam_type\n\n";
	}
	
	if ($rfam_type =~ /cis-reg/i) {
	    if (! (($rfam_name =~ /stem.loop/) || ($rfam_desc =~ /stem.loop/))) {
		
		if ($_debug) {
		    # print STDERR "$rfam_acc / $rfam_name is a cis-reg feature\n";
		}
				
		push (@$cisreg_models, $rfam_acc);

	    }
	    else {
		if ($_debug) {
		    print STDERR "$rfam_acc / $rfam_name is a stem-loop feature\n";
		}
	    }
	}
	elsif ($rfam_type =~ /sRNA/) {
	    push (@$sRNA_models, $rfam_acc);
	}
	elsif ($rfam_type =~ /antisense/i) {
	    push (@$antisense_models, $rfam_acc);
	}
	elsif ($rfam_type =~ /snorna/i) {
	    push (@$snoRNA_models, $rfam_acc);
	}
    }

}
close FULL;

if ($_debug) {
    print STDERR "found " . @$cisreg_models . " cis-reg models\n";
    print STDERR "found " . @$snoRNA_models . " snoRNA models\n";
    print STDERR "found " . @$sRNA_models . " sRNA models\n";
    print STDERR join (',', @$sRNA_models) . "\n";
}

foreach my $rfam_acc (keys %$rfam_model_name_href) {
    
    my $rfam_name = $rfam_model_name_href->{$rfam_acc};
    
    if (! defined $rfam_name) {
	die "Error, can't parse Rfam name for Rfam_acc, $rfam_acc\n";
    }
    
    # Default is ncRNA
    my $embl_feature_type = "ncRNA";
    my $embl_feature_class;
    my $biotype;
    my $kingdom = "";

    if ($rfam_name =~ /rrna/i) {
	$embl_feature_type = "rRNA";
	$biotype = "rRNA";

	$embl_feature_class = $rfam_name;
	$embl_feature_class =~ s/_rRNA/ ribosomal RNA/;
    }
    elsif ($rfam_name =~ /^tRNA$/i) {
	$embl_feature_type = "tRNA";
	$embl_feature_class = "tRNA";
	$biotype = "tRNA";
    }
    elsif ($rfam_name =~ /^tmRNA$/i) {
	$embl_feature_type = "tmRNA";
	$embl_feature_class = "tmRNA";
	$biotype = "tmRNA";
	$kingdom = "prokaryota";
    }
    elsif ($rfam_name =~ /^6S$/i) {
	$embl_feature_type = "misc_RNA";
	$embl_feature_class = "undefined";
	$biotype = "misc_RNA";
    }
    elsif ($rfam_name =~ /^Histone3|s2m|G-CSF_SLDE|AMV_RNA1_SL|Pospi_RY|HIV_GSL3|HCV_.*SL.*$/) {
	$embl_feature_type = "stem_loop";
	$embl_feature_class = "undefined";
	$biotype = "undefined";
    }
    elsif (($rfam_name =~ /ribozyme/) || ($rfam_name =~ /^hairpin$/i)) {
	$embl_feature_type = "ncRNA";
	$embl_feature_class = "ribozyme";
	$biotype = "ribozyme";
	$kingdom = "eukaryota";
    }
    elsif (is_in_array ($rfam_acc, $cisreg_models)) {
	$embl_feature_type = "misc_structure";
	$embl_feature_class = "undefined";
	$biotype = "undefined";
    }
    elsif (is_in_array ($rfam_acc, $antisense_models)) {
	$embl_feature_type = "ncRNA";
	$embl_feature_class = "antisense";
	$biotype = "antisense";
	# also archaea ???
	$kingdom = "bacteria";
    }
    elsif (is_in_array ($rfam_acc, $sRNA_models)) {
	$embl_feature_type = "misc_RNA";
	$embl_feature_class = "undefined";
	$biotype = "misc_RNA";
    }
    elsif (is_in_array ($rfam_acc, $snoRNA_models)) {
	$embl_feature_type = "ncRNA";
	$embl_feature_class = "snoRNA";
	$biotype = "snoRNA";
	$kingdom = "eukaryota";
    }

    if (($embl_feature_type eq "ncRNA") && (! defined $embl_feature_class)) {
	($embl_feature_class, $biotype, $kingdom) = define_ncRNA_class ($rfam_name);
	# print STDERR "returned kingdom, $kingdom\n";
    }

    if ($embl_feature_type eq "ncRNA") {
	if (! defined $embl_feature_class) {
	    print STDERR "could not figure out the EMBL ncRNA class ";
	    print STDERR "for Rfam model, $rfam_acc, with description, $rfam_name\n";
	    exit 1;
	}
    }

    # Finally, if ncRNA is not defined,
    # report it as a misc_RNA

    if (($embl_feature_type eq "ncRNA") && ($embl_feature_class eq "undefined")) {
	$embl_feature_type = "misc_RNA";
	# for ensembl purposes, all undefined ncRNAs get a class 'misc_RNA' which will be their ensembl biotype.
	$embl_feature_class = "misc_RNA";
	$biotype = "misc_RNA";
    }
    elsif ($embl_feature_type eq "misc_RNA") {
	# for ensembl purposes, all misc_RNAs get a class 'misc_RNA' which will be their ensembl biotype.
	$embl_feature_class = "misc_RNA";
	$biotype = "misc_RNA";
    }
    
    my $rfam_desc = $rfam_model_desc_href->{$rfam_acc};
    
    if (!defined $biotype) {
	print STDERR "not defined - $rfam_acc\t$rfam_name\t$rfam_desc\t$embl_feature_type\t$embl_feature_class\t$biotype\t$kingdom\n";
    }
    
    print STDOUT "$rfam_acc\t$rfam_name\t$rfam_desc\t$embl_feature_type\t$embl_feature_class\t$biotype\t$kingdom\n";

}

#
## The end
#


sub define_ncRNA_class {
    my ($rfam_descr) = @_;
    my $embl_feature_class;
    my $biotype;
    my $kingdom = "";

    if ($rfam_descr =~ /sno/i) {
	$embl_feature_class = "snoRNA";
	$biotype = "snoRNA";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /^U3|U8|U54|U98|SCARNA\d+/) {

	# Was guide RNA in the Rfam classification,
	# check that the type line has moved to snoRNA from guide as they said it will
	# Rfam 8.1, yes moved done
	
	$embl_feature_class = "snoRNA";
	$biotype = "snoRNA";
	$kingdom = "eukaryota";
    }
    elsif (($rfam_descr =~ /^U\d+$/) || ($rfam_descr =~ /^U\d+atac$/) || ($rfam_descr =~ /^U1_yeast$/)) {
	$embl_feature_class = "snRNA";
	$biotype = "snRNA";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /^tRNA$/) {
	# do not process it as it is not a ncRNA
	# no class
    }
    elsif ($rfam_descr =~ /vault/i) {
	$embl_feature_class = "vault_RNA";
	$biotype = "vault_RNA";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /hammerhead/i) {
	$embl_feature_class = "hammerhead_ribozyme";
	$biotype = "ribozyme";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /RNaseP/) {
	$embl_feature_class = "RNase_P_RNA";
	$biotype = "P_RNA";
	if ($rfam_descr =~ /nuc/i) {
	    $kingdom = "eukaryota";
	}
	elsif ($rfam_descr =~ /arch/i) {
	    $kingdom = "archaea";
	}
	elsif ($rfam_descr =~ /bact/i) {
	    $kingdom = "bacterial";
	}
	else {
	    die "can't define kingdom associated with rfam descr, $rfam_descr!\n";
	}
    }
    elsif ($rfam_descr =~ /^RNase_MRP$/) {
	$embl_feature_class = "RNase_MRP_RNA";
	#$biotype = "MRP_RNA";
	$biotype = "RNase_MRP_RNA";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /SRP/) {
	$embl_feature_class = "SRP_RNA";
	$biotype = "SRP_RNA";
	if ($rfam_descr =~ /dictyostelium|euk|plant|fungi|metazoa|protozoa/i) {
	    $kingdom = "eukaryota";
	}
	elsif ($rfam_descr =~ /bact/i) {
	    $kingdom = "bacterial";
	}
	elsif ($rfam_descr =~ /archaea/i) {
	    $kingdom = "archaea";
	}
	else {
	    die "can't define kingdom associated with rfam descr, $rfam_descr!\n";
	}
    }
    elsif ($rfam_descr =~ /^Y$/) {
	$embl_feature_class = "Y_RNA";
	$biotype = "Y_RNA";
	# what kingdom ?
    }
    elsif ($rfam_descr =~ /telomerase/i) {
	$embl_feature_class = "telomerase_RNA";
	$biotype = "telomerase_RNA";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /^let-7|mir-?\w+|lin-4$/i) {
	$embl_feature_class = "miRNA";
	$biotype = "miRNA";
	$kingdom = "eukaryota";
    }
    elsif ($rfam_descr =~ /intron_gp/i) {
	# Todo: Should probably be mapped to intron - it is not a gene, so should not be a ncRNA
	$embl_feature_class = "autocatalytically_spliced_intron";
	$biotype = "undefined";
    }
    # elsif ($rfam_descr =~ /^6S|DsrA|CsrB|Spot_42|GcvB|tmRNA|SECIS|Histone3$/) {
    else {
	# maybe ncRNA as a default, with an undefined class for the moment
	# or misc_RNA ???
	$embl_feature_class = "undefined";
	$biotype = "misc_RNA";
    }
    
    # print STDERR "returning kingdom, $kingdom\n";

    return ($embl_feature_class, $biotype, $kingdom);
}

sub is_in_array {
    my ($rfam_acc, $a_ref) = @_;
    
    foreach my $model (@$a_ref) {
	if ($rfam_acc eq $model) {
	    return 1;
	}
    }
    
    return 0;
    
}
