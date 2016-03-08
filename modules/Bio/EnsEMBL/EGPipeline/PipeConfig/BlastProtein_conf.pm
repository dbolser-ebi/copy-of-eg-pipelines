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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::PipeConfig::BlastProtein_conf

=head1 DESCRIPTION

Align amino acids sequences against a set of core databases.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::BlastProtein_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');
use File::Spec::Functions qw(catdir);

sub default_options {
  my $self = shift @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'blast_protein_'.$self->o('ensembl_release'),
    
    species      => [],
    antispecies  => [],
    division     => [],
    run_all      => 0,
    meta_filters => {},
    
    # Parameters for dumping and splitting Fasta query files.
    max_seq_length          => 10000000,
    max_seq_length_per_file => $self->o('max_seq_length'),
    max_seqs_per_file       => undef,
    max_files_per_directory => 100,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),
    
    max_hive_capacity => 100,
    
    # Source for proteomes can be one of: 'file', 'database', 'uniprot'
    proteome_source   => 'file',
    logic_name_prefix => $self->o('proteome_source'),
    
    # This parameter should only be specified if proteome_source = 'file'
    db_fasta_file => undef,
    
    # Dump the proteomes of other species and align them.
    source_species      => [],
    source_antispecies  => [],
    source_division     => [],
    source_run_all      => 0,
    source_meta_filters => {},
    
    # Taxonomic level is one of 'fungi', 'invertebrates', 'plants';
    # the UniProt source can be 'sprot' or 'trembl'.
    uniprot_ftp_uri  => 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions',
    taxonomic_level  => 'invertebrates',
    uniprot_source   => 'sprot',
    taxonomic_levels => [$self->o('taxonomic_level')],
    uniprot_sources  => [$self->o('uniprot_source')],
    uniprot_dir      => catdir($self->o('pipeline_dir'), 'uniprot'),
    
    # Default is to use ncbi-blast; wu-blast is possible, for backwards
    # compatability, but it is now unsupported and orders of magnitude
    # slower than ncbi-blast. Setting BLASTMAT is irrelevant for ncbi-blast,
    # but the Ensembl code requires it to be set as an environment variable,
    # pointing to a directory that exists.
    blast_type       => 'ncbi',
    blast_dir        => '/nfs/panda/ensemblgenomes/external/ncbi-blast-2+',
    makeblastdb_exe  => catdir($self->o('blast_dir'), 'bin/makeblastdb'),
    blastp_exe       => catdir($self->o('blast_dir'), 'bin/blastp'),
    blastx_exe       => catdir($self->o('blast_dir'), 'bin/blastx'),
    blast_matrix     => undef,
    blast_threads    => 3,
    blast_parameters => '-word_size 3 -num_alignments 100000 -num_descriptions 100000 -lcase_masking -seg yes -num_threads '.$self->o('blast_threads'),
    
    # For wu-blast, set the following  parameters instead of the above.
    # blast_type       => 'wu',
    # blast_dir        => '/nfs/panda/ensemblgenomes/external/wublast',
    # makeblastdb_exe  => catdir($self->o('blast_dir'), 'wu-formatdb'),
    # blastp_exe       => catdir($self->o('blast_dir'), 'blastp'),
    # blastx_exe       => catdir($self->o('blast_dir'), 'blastx'),
    # blast_matrix     => catdir($self->o('blast_dir'), 'matrix'),
    # blast_parameters => '-W 3 -B 100000 -V 100000 -hspmax=0 -lcmask -wordmask=seg',
    
    # For parsing the output.
    output_regex     => '^\s*(\w+)',
    pvalue_threshold => 0.01,
    filter_prune     => 1,
    filter_min_score => 200,
    
    # By default, do both blastp and blastx.
    # Specify blast_db if you don't want it in the directory alongside db_file.
    blast_db => undef,
    blastp   => 1,
    blastx   => 1,
    
    # For blastx results, filter so that only best hit per seq_region gets saved.
    unique_hits => 1,
    
    # Generate a GFF file for loading into, e.g., WebApollo
    create_gff    => 0,
    gt_exe        => '/nfs/panda/ensemblgenomes/external/genometools/bin/gt',
    gff3_tidy     => $self->o('gt_exe').' gff3 -tidy -sort -retainids',
    gff3_validate => $self->o('gt_exe').' gff3validator',
    
    analyses =>
    [
      {
        'logic_name'    => 'blastp',
        'db'            => $self->o('proteome_source'),
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
      },
      
      {                 
        'logic_name'    => 'blastx',
        'db'            => $self->o('proteome_source'),
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'otherfeatures',
      },
      
    ],

    # Remove existing *_align_features; if => 0 then existing analyses
    # and their features will remain, with the logic_name suffixed by '_bkp'.
    delete_existing => 1,

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 0,

  };
}

sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  
  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );
  
  return $options;
}

sub hive_meta_table {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  
  return [
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('pipeline_dir'),
    $self->db_cmd("CREATE TABLE split_genome (species varchar(100) NOT NULL, split_file varchar(255) NOT NULL)"),
    $self->db_cmd("CREATE TABLE split_proteome (species varchar(100) NOT NULL, split_file varchar(255) NOT NULL)"),
  ];
}

sub pipeline_analyses {
  my $self = shift @_;
  
  return [
    {
      -logic_name      => 'InitialisePipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -input_ids       => [{}],
      -flow_into       => {
                            '1->A' => ['SpeciesFactory'],
                            'A->1' => ['BlastProtein'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count => 1,
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            variation_flow  => 0,
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['TargetDatabase'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'TargetDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::TargetDatabase',
      -max_retry_count => 1,
      -parameters      => {
                            blastp => $self->o('blastp'),
                            blastx => $self->o('blastx'),
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2->A' => ['BackupCoreDatabase'],
                            'A->2' => ['DumpProteome'],
                            '3->B' => ['CheckOFDatabase'],
                            'B->3' => ['DumpGenome'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'BackupCoreDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -can_be_empty      => 1,
      -max_retry_count => 1,
      -parameters      => {
                            output_file => catdir($self->o('pipeline_dir'), '#species#', 'core_bkp.sql.gz'),
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'CheckOFDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::CheckOFDatabase',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {},
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BackupOFDatabase'],
                            '3' => ['CreateOFDatabase'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'BackupOFDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DatabaseDumper',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            db_type     => 'otherfeatures',
                            output_file => catdir($self->o('pipeline_dir'), '#species#', 'of_bkp.sql.gz'),
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'CreateOFDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::CreateOFDatabase',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {},
      -rc_name         => 'normal',
      -flow_into       => ['BackupOFDatabase'],
    },

    {
      -logic_name        => 'DumpProteome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpProteome',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              proteome_dir => catdir($self->o('pipeline_dir'), '#species#', 'proteome'),
                              use_dbID     => 1,
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SplitProteome'],
    },

    {
      -logic_name        => 'SplitProteome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              fasta_file              => '#proteome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => [ ':////split_proteome' ],
                            }
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              genome_dir => catdir($self->o('pipeline_dir'), '#species#', 'genome'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SplitGenome'],
    },

    {
      -logic_name        => 'SplitGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => [ ':////split_genome' ],
                            }
    },

    {
      -logic_name      => 'BlastProtein',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastProtein',
      -max_retry_count => 0,
      -parameters      => {
                            proteome_source => $self->o('proteome_source'),
                            db_fasta_file   => $self->o('db_fasta_file'),
                          },
      -flow_into       => {
                            '2' => ['FetchFile'],
                            '3' => ['SourceSpeciesFactory'],
                            '4' => ['FetchUniprot'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'FetchFile',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -can_be_empty    => 1,
      -max_retry_count => 0,
      -flow_into       => ['CreateBlastDB'],
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'SourceSpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            species         => $self->o('source_species'),
                            antispecies     => $self->o('source_antispecies'),
                            division        => $self->o('source_division'),
                            run_all         => $self->o('source_run_all'),
                            meta_filters    => $self->o('source_meta_filters'),
                            chromosome_flow => 0,
                            variation_flow  => 0,
                            species_varname => 'source_species',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['FetchDatabase'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name        => 'FetchDatabase',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpProteome',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              species      => '#source_species#',
                              proteome_dir => catdir($self->o('pipeline_dir'), 'proteomes'),
                              file_varname => 'db_fasta_file',
                            },
      -rc_name           => 'normal',
      -flow_into         => ['CreateBlastDB'],
    },

    {
      -logic_name      => 'FetchUniprot',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot',
      -can_be_empty    => 1,
      -max_retry_count => 2,
      -parameters      => {
                            ftp_uri          => $self->o('uniprot_ftp_uri'),
                            taxonomic_levels => $self->o('taxonomic_levels'),
                            uniprot_sources  => $self->o('uniprot_sources'),
                            out_dir          => $self->o('uniprot_dir'),
                          },
      -rc_name         => 'normal',
      -flow_into       => ['CreateBlastDB'],
    },

    {
      -logic_name      => 'CreateBlastDB',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::CreateBlastDB',
      -max_retry_count => 2,
      -parameters      => {
                            makeblastdb_exe   => $self->o('makeblastdb_exe'),
                            blast_db_type     => 'prot',
                            proteome_source   => $self->o('proteome_source'),
                            logic_name_prefix => $self->o('logic_name_prefix'),
                          },
      -rc_name         => 'normal',
      -flow_into       => ['TargetSpeciesFactory'],
    },

    {
      -logic_name      => 'TargetSpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count => 1,
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            variation_flow  => 0,
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['AnalysisFactory'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'AnalysisFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::AnalysisFactory',
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            analyses => $self->o('analyses'),
                            blastp   => $self->o('blastp'),
                            blastx   => $self->o('blastx'),
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2->A' => ['AnalysisSetupCore'],
                            'A->4' => ['FetchProteomeFiles'],
                            '3->B' => ['AnalysisSetupOF'],
                            'B->5' => ['FetchGenomeFiles'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'AnalysisSetupCore',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -can_be_empty    => 1,
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            db_backup_required => 1,
                            db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'core_bkp.sql.gz'),
                            delete_existing    => $self->o('delete_existing'),
                            production_lookup  => $self->o('production_lookup'),
                            production_db      => $self->o('production_db'),
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'AnalysisSetupOF',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup',
      -can_be_empty    => 1,
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            db_type            => 'otherfeatures',
                            db_backup_required => 1,
                            db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'of_bkp.sql.gz'),
                            delete_existing    => $self->o('delete_existing'),
                            production_lookup  => $self->o('production_lookup'),
                            production_db      => $self->o('production_db'),
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'FetchProteomeFiles',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchSplitFiles',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            seq_type => 'proteome',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2->A' => ['BlastPFactory'],
                            'A->1' => ['UniqueBlastPHits'],
                          },
    },

    {
      -logic_name      => 'FetchGenomeFiles',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchSplitFiles',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            seq_type => 'genome',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2->A' => ['BlastXFactory'],
                            'A->1' => ['UniqueBlastXHits'],
                          },
    },

    {
      -logic_name      => 'BlastPFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastFactory',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            max_seq_length => $self->o('max_seq_length'),
                            queryfile      => '#split_file#',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BlastP'],
                          },
    },

    {
      -logic_name      => 'BlastXFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastFactory',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            max_seq_length => $self->o('max_seq_length'),
                            queryfile      => '#split_file#',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BlastX'],
                          },
    },

    {
      -logic_name      => 'BlastP',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            logic_name       => '#logic_name_prefix#_blastp',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'pep',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '8GB_threads',
      -flow_into       => {
                            '-1' => ['BlastP_HighMem'],
                          },
    },

    {
      -logic_name      => 'BlastP_HighMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            logic_name       => '#logic_name_prefix#_blastp',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'pep',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '16GB_threads',
      -flow_into       => {
                            '-1' => ['BlastP_HigherMem'],
                          },
    },

    {
      -logic_name      => 'BlastP_HigherMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            logic_name       => '#logic_name_prefix#_blastp',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'pep',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                          },
      -rc_name         => '32GB_threads',
    },

    {
      -logic_name      => 'BlastX',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'otherfeatures',
                            logic_name       => '#logic_name_prefix#_blastx',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'dna',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '8GB_threads',
      -flow_into       => {
                            '-1' => ['BlastX_HighMem'],
                          },
    },

    {
      -logic_name      => 'BlastX_HighMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'otherfeatures',
                            logic_name       => '#logic_name_prefix#_blastx',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'dna',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '16GB_threads',
      -flow_into       => {
                            '-1' => ['BlastX_HigherMem'],
                          },
    },

    {
      -logic_name      => 'BlastX_HigherMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'otherfeatures',
                            logic_name       => '#logic_name_prefix#_blastx',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'dna',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                          },
      -rc_name         => '32GB_threads',
    },

    {
      -logic_name      => 'UniqueBlastPHits',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::UniqueHits',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            db_type     => 'core',
                            unique_hits => $self->o('unique_hits'),
                            logic_name  => '#logic_name_prefix#_blastx',
                          },
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'UniqueBlastXHits',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::UniqueHits',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            db_type     => 'otherfeatures',
                            unique_hits => $self->o('unique_hits'),
                            logic_name  => '#logic_name_prefix#_blastx',
                            create_gff  => $self->o('create_gff'),
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['GFF3Dump'],
                          },
    },

    {
      -logic_name        => 'GFF3Dump',
      -module            => 'Bio::EnsEMBL::EGPipeline::FileDump::GFF3Dumper',
      -analysis_capacity => 10,
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              db_type           => 'otherfeatures',
                              feature_type      => ['ProteinAlignFeature'],
                              include_scaffold  => 0,
                              remove_separators => 1,
                              results_dir       => $self->o('pipeline_dir'),
                              out_file_stem     => '#logic_name_prefix#_blastx.gff3',
                            },
      -rc_name           => 'normal',
      -flow_into         => ['GFF3Validate'],
    },

    {
      -logic_name        => 'GFF3Validate',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -analysis_capacity => 10,
      -can_be_empty      => 1,
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              cmd => $self->o('gff3_tidy').' #out_file# > #out_file#.sorted; '.
                                      'mv #out_file#.sorted #out_file#; '.
                                      $self->o('gff3_validate').' #out_file#',
                            },
      -rc_name           => 'normal',
    }

  ];
}

sub resource_classes {
  my ($self) = @_;
  
  my $blast_threads = $self->o('blast_threads');

  return {
    %{$self->SUPER::resource_classes},
    '8GB_threads' => {'LSF' => '-q production-rh6 -n ' . ($blast_threads + 1) . ' -R "span[hosts=1]" -M 8000 -R "rusage[mem=8000,tmp=8000]"'},
    '16GB_threads' => {'LSF' => '-q production-rh6 -n ' . ($blast_threads + 1) . ' -R "span[hosts=1]" -M 16000 -R "rusage[mem=16000,tmp=16000]"'},
    '32GB_threads' => {'LSF' => '-q production-rh6 -n ' . ($blast_threads + 1) . ' -R "span[hosts=1]" -M 32000 -R "rusage[mem=32000,tmp=32000]"'},
  }
}

1;