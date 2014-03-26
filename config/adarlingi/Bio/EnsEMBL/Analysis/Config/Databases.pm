=head1 LICENSE

Copyright [1999-2014] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 NAME

    Bio::EnsEMBL::Analysis::Config::Databases

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::Databases;

=head1 DESCRIPTION

    Databases.pm is the main configuration file which holds the different
    parameters (usernames, hosts, passwords, ports, database-names) to
    connect to different databases used in the Ensembl-Analysis pipeline.

    It imports and sets a number of standard global variables into the
    calling package. Without arguments all the standard variables are set,
    and with a list, only those variables whose names are provided are
    set. The module will die if a variable which doesn't appear in its
    C<%Config> hash is asked to be set.

    A common way to get an DBAdaptor in a module is:

        print "Loading database : ".  $$DATABASES{REFERENCE_DB}{"-dbname"} . "\n";

        my $ref_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $$DATABASES{REFERENCE_DB} } );

    OR if you write a RunnableDB:

        use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
        @ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild );

        my $genes_db = $self->get_dbadaptor("GENEBUILD_DB");

    OR for normal scripts :

        use Bio::EnsEMBL::Analysis::Tools::Utilities qw( get_db_adaptor_by_string );
        get_db_adaptor_by_string("GENEBUILD_DB");

    The variables can also be references to arrays or hashes.

    Edit C<%Config> to add or alter variables.

    All the variables are in capitals, so that they resemble
    environment variables.

    Databases is a pure ripoff of humConf written by James Gilbert.
    humConf is based upon ideas from the standard perl Env environment
    module.

=head1 LICENSE

    Copyright (c) 1999-2009 The European Bioinformatics Institute and
    Genome Research Limited.  All rights reserved.

    This software is distributed under a modified Apache license.
    For license details, please see
      http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

    Please email comments or questions to the public Ensembl
    developers list at <ensembl-dev@ebi.ac.uk>.

    Questions may also be sent to the Ensembl help desk at
    <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Analysis::Config::Databases;

use strict;
use vars qw(%Config);

%Config = (

  DATABASES => {

    # The REFERENCE_DB (formely known as GB_DB) holds sequence,
    # repeats and features from raw computes (e.g. ab-inito
    # predictions, dna- or protein alignments )

    REFERENCE_DB => { 
                      -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                      -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                      -dbport => '4126',
                      -dbuser => 'ensrw',
                      -dbpass => 'scr1b3d1',
    },

    # The GENEWISE_DB holds genes made by FPC_TargettedGenewise or
    # FPC_BlastMiniGenewise (TGE_gw or similarity_genewise - genes ) (
    # formerly GB_GW_DB )

    GENEWISE_DB => { 
                     -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                     -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                     -dbport => '4126',
                     -dbuser => 'ensrw',
                     -dbpass => 'scr1b3d1',
    },

    # The EXONERATE_DB ( formerly GB_cDNA ) holds alignments to cDNA's
    # or EST's gene-structtures made by exonerate (Exonerate2Genes.pm)

    EXONERATE_DB => {
                      -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                      -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                      -dbport => '4126',
                      -dbuser => 'ensrw',
                      -dbpass => 'scr1b3d1',
    },

    # The BLESSED_DB (formerly GB_BLESSED) holds the 'blessed'
    # gene-set ( if there is one )

    BLESSED_DB => { 
                    -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                    -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                    -dbport => '4126',
                    -dbuser => 'ensrw',
                    -dbpass => 'scr1b3d1',
    },

    # The UTR_DB (formerly GB_COMB) holds genes made by the
    # UTR-addtion-run Combine_Genewises_and_E2Gs.pm writes to UTR_DB

    UTR_DB => { 
                -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                -dbport => '4126',
                -dbuser => 'ensrw',
                -dbpass => 'scr1b3d1',
    },

    # GENEBUILD_DB (formerly GB_FINALDB) is the Database where
    # GeneBuilder.pm writes its results to this database and
    # The Pseudogene-code READS from this database

    GENEBUILD_DB => { 
                      -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                      -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                      -dbport => '4126',
                      -dbuser => 'ensrw',
                      -dbpass => 'scr1b3d1',
    },

    # PSEUDO_DB holds the pseudo-genes

    PSEUDO_DB => { 
                   -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                   -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                   -dbport => '4126',
                   -dbuser => 'ensrw',
                   -dbpass => 'scr1b3d1',
    },

    # COALESCER_DB is the DB where TranscriptCoalescer writes its
    # results to

    COALESCER_DB => { 
                      -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                      -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                      -dbport => '4126',
                      -dbuser => 'ensrw',
                      -dbpass => 'scr1b3d1',
    },

    # OrthologueEvalutor, FindMissingGenes, etc. write it's results into
    # ORTHOLOGUE_DB

    ORTHOLOGUE_DB => { 
                       -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                       -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                       -dbport => '4126',
                       -dbuser => 'ensrw',
                       -dbpass => 'scr1b3d1',
    },

    # This database is use mainly with human an mouse for merging the
    # havana get set with the ensembl genebuild

    HAVANA_DB => { 
                   -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                   -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                   -dbport => '4126',
                   -dbuser => 'ensrw',
                   -dbpass => 'scr1b3d1',
    },

    # This database is use mainly with human variation genewise alignments

    VARIATION_DB => {
                      -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                      -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                      -dbport => '4126',
                      -dbuser => 'ensrw',
                      -dbpass => 'scr1b3d1',
    },

    # Database which stores the ditags...

    DITAG_DB => {
                  -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                  -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                  -dbport => '4126',
                  -dbuser => 'ensrw',
                  -dbpass => 'scr1b3d1',
    },

    #############
    ### Databases for RunnableDB/CopyGenes.pm
    #############

    COPY_SOURCE_DB => { 
                        -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                        -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                        -dbport => '4126',
                        -dbuser => 'ensrw',
                        -dbpass => 'scr1b3d1',
    },

    COPY_TARGET_DB => {
                        -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                        -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                        -dbport => '4126',
                        -dbuser => 'ensrw',
                        -dbpass => 'scr1b3d1',
    },

    KILL_LIST_DB => { 
                      -dbname => 'mhinsley_anopheles_darlingi_core_20_73_1',
                      -dbhost => 'mysql-eg-devel-1.ebi.ac.uk',
                      -dbport => '4126',
                      -dbuser => 'ensrw',
                      -dbpass => 'scr1b3d1',
    },
  },
  #
  # database with dna + repeat feature and repeat consensus table 
  #

  DNA_DBNAME => "REFERENCE_DB",

  #
  # Arrayref. to distributed dbs ( see POD on top of file )   
  #

  DISTRIBUTED_DBS => {
                       DITAG_DB_DIST => [ "DB_COPY_1","DB_COPY_2"],
                     },
  #
  # MAIN_REFERENCE_DB should point to the database which is used as reference db. The MAIN_REFERENCE_DB
  # key is used in some cloud-based modules to automatically create new databaes. The 'vital tables' array
  # specifies which tables are dumped by some modules.
  #

  MAIN_REFERENCE_DB => "REFERENCE_DB",
  VITAL_TABLES => ["analysis",
                   "analysis_description",
                   "assembly",
                   "assembly_exception",
                   "attrib_type",
                   "coord_system",
                   "external_db",
                   "meta",
                   "meta_coord",
                   "seq_region",
                   "seq_region_attrib",
                  ],

);

sub import {
  my ($callpack) = caller(0);    # Name of the calling package
  my $pack = shift;              # Need to move package off @_

  # Get list of variables supplied, or else all
  my @vars = @_ ? @_ : keys(%Config);
  return unless @vars;

  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join( ' ', map { '$' . $_ } @vars ) . ")";
  die $@ if $@;

  foreach (@vars) {
    if ( defined $Config{$_} ) {
      no strict 'refs';
      # Exporter does a similar job to the following
      # statement, but for function names, not
      # scalar variables:
      *{"${callpack}::$_"} = \$Config{$_};
    } else {
      die "Error: Config: $_ not known\n";
    }
  }
} ## end sub import

1;
