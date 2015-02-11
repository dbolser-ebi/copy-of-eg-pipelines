#!/usr/bin/env perl
=head1 LICENSE

  Copyright (c) 1999-2014 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 DESCRIPTION

This script is used to add UniProt cross-references to one or more Ensembl cores, using UniParc annotation

=head1 EXAMPLE

 perl -I modules scripts/xrefs_pipeline/load_uniprot.pl -host mysql-eg-devel-1 -port 4126
  -user ensrw -pass scr1b3d1 -dbname  schizosaccharomyces_pombe_core_22_75_2 -uniprothost whisky.ebi.ac.uk
   -uniprotport 1531 -uniprotuser spselect -uniprotpass spselect -uniprotdbname SWPREAD
    -uniprotdriver Oracle -uniparchost ora-vm-004.ebi.ac.uk -uniparcport 1551
     -uniparcuser uniparc_read -uniparcpass uniparc -uniparcdbname UAPRO -uniparcdriver Oracle

=head1 USAGE

  --user=user                      username for the core database server

  --pass=pass                      password for core database server

  --host=host                      release core server 

  --port=port                      port for release database server 
  
  --dbname=dbname                  name of core database
  
  --rfamdriver=dbname           driver to use for rfam database

  --rfamuser=user               username for the rfam database

  --rfampass=pass               password for rfam database

  --rfamhost=host               server where the rfam database is stored

  --rfamport=port               port for rfam database
  
  --rfamdbname=dbname           name/SID of rfam database to process
  
  --verbose                        Increase logging level to debug

=head1 AUTHOR

dstaines

=cut

use warnings;
use strict;

use Bio::EnsEMBL::Utils::CliHelper;
use Log::Log4perl qw/:easy/;
use Pod::Usage;
use Bio::EnsEMBL::EGPipeline::Rfam::RfamLoader;
use Data::Dumper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();
# get the basic options for connecting to a database server
my $optsd = [@{$cli_helper->get_dba_opts()}, @{$cli_helper->get_dba_opts('rfam')}];
push(@{$optsd}, "verbose");

my $opts = $cli_helper->process_args($optsd, \&pod2usage);

if ($opts->{verbose}) {
  Log::Log4perl->easy_init($DEBUG);
} else {
  Log::Log4perl->easy_init($INFO);
}

my $logger = get_logger();

$logger->info("Connecting to Rfam database");
my ($rfam_dba) = @{$cli_helper->get_dbas_for_opts($opts, 1, 'rfam')};

my $loader = Bio::EnsEMBL::EGPipeline::Rfam::RfamLoader->new(
    -RFAM_DBA => $rfam_dba
);

$logger->info("Connecting to core database(s)");
for my $core_dba_details (@{$cli_helper->get_dba_args_for_opts($opts)}) {
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$core_dba_details});
  $logger->info("Processing " . $dba->species());
  $loader->add_genes($dba);
}
