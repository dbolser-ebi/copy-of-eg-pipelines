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


# Add RFAM xrefs to rRNAs predicted by RNAmmer,
# Required by the ncRNA tree pipeline from Ensembl Compara

use warnings;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBEntry;

my ($dbhost, $dbport, $dbuser, $dbpass, $species);

# Get the options

GetOptions(
	   'dbhost=s'    => \$dbhost,
	   'species=s'   => \$species,
	   'dbuser=s'    => \$dbuser,
	   'dbpass=s'    => \$dbpass,
	   'dbport=s'    => \$dbport,
    );


# Load the Registry

Bio::EnsEMBL::Registry->load_registry_from_db (
    -host => $dbhost,
    -user => $dbuser,
    -pass => $dbpass,
    -port => $dbport,
    );

if (!defined $dbhost || !defined $dbport || !defined $dbuser || !defined $dbpass || !defined $species) {
    die "missing argument!\n";
}

# They don't have 28S in RFAM ?

my $mapping_href = {
    '5S_rRNA'   => 'RF00001',
    '5.8S_rRNA' => 'RF00002',
    '18S_rRNA'  => 'RF01960',
};
my $rfams_href = {
    'RF00001' => {'name' => '5S_rRNA', 'desc' => '5S ribosomal RNA'},
    'RF00002' => {'name' => '5.8S_rRNA', 'desc' => '5.8S ribosomal RNA'},
    'RF01960' => {'name' => 'SSU_rRNA_eukarya', 'desc' => 'Eukaryotic small subunit ribosomal RNA'},
};

my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor("$species", "core", "gene");
my $dbentry_adaptor = Bio::EnsEMBL::Registry->get_adaptor("$species", "core", "DBEntry");

if (!defined $gene_adaptor) {
    die "can't get gene adaptor for species, $species, for server, $dbhost!\n";
}

my $genes_aref = $gene_adaptor->fetch_all_by_biotype('rRNA');

print STDERR "Processing " . @$genes_aref. " rRNA genes\n";

foreach my $gene (@$genes_aref) {
    
    if (defined $gene->external_name()) {
	
	# Get the gene name if there is one
	
	my $display_name = $gene->external_name();
	
	# Skip 28S
	# 'PK-G12rRNA' from RFAM, so already there
	
	if ($display_name eq "28S_rRNA" || $display_name eq "PK-G12rRNA") {
	    next;
	}
	
	# Map to Rfam
	
	if (!defined $mapping_href->{$display_name}) {
	    die "can't get rfam_acc mapping for gene, '" . $gene->stable_id() . "' with display_name, '$display_name'!\n";
	}
	my $rfam_acc = $mapping_href->{$display_name};
	my $rfam_entry_href = $rfams_href->{$rfam_acc};
	
	# Create a DBEntry and store it, attach to the current gene object
	
	my $db_entry = Bio::EnsEMBL::DBEntry->new (
	    -primary_id => $rfam_acc,
	    -version => 1,
	    -dbname  => 'RFAM',
	    -display_id => $rfam_entry_href->{'name'},
	    -description => $rfam_entry_href->{'desc'},
	    -info_type => 'DIRECT',
	    );
	
	# print STDERR "Attaching dbEntry, '$rfam_acc', to gene with internal_id, " . $gene->dbID() . "\n";
	$gene->add_DBEntry($db_entry);
	
	$dbentry_adaptor->store($db_entry,$gene->dbID(),'Gene',1);
    }
}

