=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::PipeConfig::ShortReadAlignment_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.3;
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');
use File::Spec::Functions qw(catdir);

# To-do list:
# Work out how to stop STAR from filtering short reads in EST-mode.
# For STAR report its statistics rather than the less meaningful BAM stats.
# Allow runs to be specified instead of studies

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    pipeline_name => 'short_read_alignment_'.$self->o('ensembl_release'),

    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},

    # This pipeline can align data from one or more files, or direct
    # from ENA; it _could_ use data from both sources, but you're liable
    # to get in a muddle if you do that, so it's not recommended.
    seq_file         => [],
    seq_file_pair    => [],
    run              => [],
    study            => [],
    merge_level      => 'run',
    tax_id_restrict  => 1,
    
    # RNA-seq options
    ini_type      => 'rnaseq_align',
    bigwig        => 0,
    vcf           => 0,
    use_csi       => 0,
    clean_up      => 1,
    
    # Parameters for dumping and splitting Fasta DNA query files.
    max_seq_length_per_file => 30000000,
    max_seqs_per_file       => undef,
    max_files_per_directory => 50,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),

    # Parameters for repeatmasking the genome files.
    repeat_masking     => 'soft',
    repeat_logic_names => [],
    min_slice_length   => 0,

    # Aligner options.
    aligner    => 'star',
    threads    => 4,
    data_type  => 'rnaseq',
    read_type  => 'default',

    logic_name => $self->o('data_type').'_'.$self->o('aligner'),

    # Some of the aligners have newer versions, but it's not a given that
    # these will be better than the version we've used up till now. So the
    # latter is the default, but you can experiment with the latest versions
    # by commenting/uncommenting below.
    bowtie2_dir  => '/nfs/panda/ensemblgenomes/external/bowtie2-2.2.6',
    
    bwa_dir      => '/nfs/panda/ensemblgenomes/external/bwa',
    #bwa_dir      => '/nfs/panda/ensemblgenomes/external/bwa0.7.12_x64-Linux',
    
    gsnap_dir    => '/nfs/panda/ensemblgenomes/external/gmap-gsnap/bin',
    #gsnap_dir    => '/nfs/panda/ensemblgenomes/external/gmap-gsnap-2015-11-20/bin',
    
    star_dir     => '/nfs/panda/ensemblgenomes/external/STAR',
    #star_dir     => '/nfs/panda/ensemblgenomes/external/STAR_2.4.2a.Linux_x86_64',

    # Different aligners have different memory requirements; unless explicitly
    # over-ridden, use defaults, which should work on a genome that isn't too
    # fragmented, of size < 1Gb. (Values here are MB.)
    index_memory      => undef,
    index_memory_high => undef,
    align_memory      => undef,
    align_memory_high => undef,
    
    index_memory_default => {
      'bowtie2' =>  8000,
      'bwa'     => 16000,
      'gsnap'   => 16000,
      'star'    => 32000,
    },
    index_memory_high_default => {
      'bowtie2' => 16000,
      'bwa'     => 32000,
      'gsnap'   => 32000,
      'star'    => 64000,
    },
    align_memory_default => {
      'bowtie2' =>  8000,
      'bwa'     => 32000,
      'gsnap'   => 32000,
      'star'    => 16000,
    },
    align_memory_high_default => {
      'bowtie2' => 16000,
      'bwa'     => 64000,
      'gsnap'   => 64000,
      'star'    => 32000,
    },

    samtools_dir  => '/nfs/panda/ensemblgenomes/external/samtools',
    bedtools_dir  => '/nfs/panda/ensemblgenomes/external/bedtools/bin',
    ucscutils_dir => '/nfs/panda/ensemblgenomes/external/ucsc_utils',

    # Remove existing alignments; if => 0 then existing analyses
    # and their features will remain, with the logic_name suffixed by '_bkp'.
    delete_existing => 1,

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 1,
  };
}

# Force an automatic loading of the registry in all workers.
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;

  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );
  
  return $options;
}

# Ensures that species output parameter gets propagated implicitly.
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
    'mkdir -p '.catdir($self->o('pipeline_dir'), $self->o('aligner')),
  ];
}

sub pipeline_analyses {
  my ($self) = @_;
  
  # The analyses are defined within a function, to allow inheriting conf
  # files to easily modify the core functionality of this pipeline.
  my $alignment_analyses = $self->alignment_analyses();
  $self->modify_analyses($alignment_analyses);
  
  return $alignment_analyses;
}

sub aligner_parameters {
  my ($self, $aligner, $data_type, $read_type) = @_;
  
  my %aligner_classes =
  (
    'bowtie2' => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::Bowtie2Aligner',
    'bwa'     => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::BwaAligner',
    'gsnap'   => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::GsnapAligner',
    'star'    => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::StarAligner',
  );
  my $aligner_class = $aligner_classes{$aligner};
  
  my %aligner_dirs =
  (
    'bowtie2' => $self->o('bowtie2_dir'),
    'bwa'     => $self->o('bwa_dir'),
    'gsnap'   => $self->o('gsnap_dir'),
    'star'    => $self->o('star_dir'),
  );
  my $aligner_dir = $aligner_dirs{$aligner};
  
  $read_type = 'long_reads' if $data_type !~ /rna_?seq/i;
  
  return ($aligner_class, $aligner_dir, $read_type);
}

sub alignment_analyses {
  my ($self) = @_;
  
  my ($aligner_class, $aligner_dir, $read_type) =
    $self->aligner_parameters(
      $self->o('aligner'),
      $self->o('data_type'),
      $self->o('read_type')
    );
  
  my $dir = catdir($self->o('pipeline_dir'), $self->o('aligner'));
  
  return
  [
    {
      -logic_name        => 'SpeciesFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              species         => $self->o('species'),
                              antispecies     => $self->o('antispecies'),
                              division        => $self->o('division'),
                              run_all         => $self->o('run_all'),
                              meta_filters    => $self->o('meta_filters'),
                              chromosome_flow => 0,
                              variation_flow  => 0,
                            },
      -input_ids         => [ {} ],
      -flow_into         => {
                              '2' => ['DumpGenome'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -analysis_capacity => 5,
      -batch_size        => 2,
      -max_retry_count   => 1,
      -parameters        => {
                              genome_dir         => catdir($dir, '#species#'),
                              repeat_masking     => $self->o('repeat_masking'),
                              repeat_logic_names => $self->o('repeat_logic_names'),
                              min_slice_length   => $self->o('min_slice_length'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['BigWigFlowControl'],
    },

    {
      -logic_name        => 'BigWigFlowControl',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FlowControl',
      -max_retry_count   => 1,
      -parameters        => {
                              control_value => $self->o('bigwig'),
                              control_flow  => { '0' => '1', '1' => '2', },
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '1' => ['IndexGenome'],
                              '2' => ['SequenceLengths'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'SequenceLengths',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceLengths',
      -analysis_capacity => 5,
      -batch_size        => 2,
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              fasta_file  => '#genome_file#',
                              length_file => '#genome_file#'.'.lengths.txt',
                            },
      -rc_name           => 'normal',
      -flow_into         => ['IndexGenome'],
    },

    {
      -logic_name        => 'IndexGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::IndexGenome',
      -analysis_capacity => 5,
      -batch_size        => 2,
      -max_retry_count   => 1,
      -parameters        => {
                              aligner_class => $aligner_class,
                              aligner_dir   => $aligner_dir,
                              samtools_dir  => $self->o('samtools_dir'),
                              threads       => $self->o('threads'),
                              memory_mode   => 'default',
                              escape_branch => -1,
                            },
      -rc_name           => 'index_default',
      -flow_into         => {
                              '-1' => ['IndexGenome_HighMem'],
                               '1' => ['SequenceFactory'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'IndexGenome_HighMem',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::IndexGenome',
      -analysis_capacity => 5,
      -batch_size        => 2,
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              aligner_class => $aligner_class,
                              aligner_dir   => $aligner_dir,
                              samtools_dir  => $self->o('samtools_dir'),
                              threads       => $self->o('threads'),
                              memory_mode   => 'himem',
                            },
      -rc_name           => 'index_himem',
      -flow_into         => {
                              '1' => ['SequenceFactory'],
                            },
    },

    {
      -logic_name        => 'SequenceFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              seq_file         => $self->o('seq_file'),
                              seq_file_pair    => $self->o('seq_file_pair'),
                              run              => $self->o('run'),
                              study            => $self->o('study'),
                              merge_level      => $self->o('merge_level'),
                              data_type        => $self->o('data_type'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '3' => ['SplitSeqFile'],
                              '4' => ['PairSeqFile'],
                              '5' => ['SRASeqFile'],
                            },
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name        => 'SplitSeqFile',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SplitSeqFile',
      -analysis_capacity => 5,
      -batch_size        => 4,
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                              out_dir                 => catdir($dir, '#species#', 'seqs'),
                              delete_existing_files   => 0,
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2->A' => ['AlignSequence'],
                              'A->1' => ['MergeBam'],
                            },
    },

    {
      -logic_name        => 'PairSeqFile',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -can_be_empty      => 1,
      -parameters        => {},
      -rc_name           => 'normal',
      -flow_into         => {
                              '1->A' => ['AlignSequence'],
                              'A->1' => ['MergeBam'],
                            },
    },

    {
      -logic_name        => 'SRASeqFile',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SRASeqFile',
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              work_directory  => catdir($dir, '#species#'),
                              merge_level     => $self->o('merge_level'),
                              tax_id_restrict => $self->o('tax_id_restrict'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2->A' => ['AlignSequence'],
                              'A->1' => ['MergeBam'],
                            },
    },

    {
      -logic_name        => 'AlignSequence',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::AlignSequence',
      -analysis_capacity => 25,
      -max_retry_count   => 1,
      -parameters        => {
                              aligner_class  => $aligner_class,
                              aligner_dir    => $aligner_dir,
                              samtools_dir   => $self->o('samtools_dir'),
                              threads        => $self->o('threads'),
                              read_type      => $read_type,
                              clean_up       => $self->o('clean_up'),
                              escape_branch  => -1,
                            },
      -rc_name           => 'align_default',
      -flow_into         => {
                              '-1' => ['AlignSequence_HighMem'],
                               '1' => {
                                        ':////accu?merge_ids={bam_file}' => {'merge_ids' => '#merge_id#'},
                                      },
                            },
    },

    {
      -logic_name        => 'AlignSequence_HighMem',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::AlignSequence',
      -analysis_capacity => 25,
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              aligner_class  => $aligner_class,
                              aligner_dir    => $aligner_dir,
                              samtools_dir   => $self->o('samtools_dir'),
                              threads        => $self->o('threads'),
                              read_type      => $read_type,
                              clean_up       => $self->o('clean_up'),
                            },
      -rc_name           => 'align_himem',
      -flow_into         => {
                               '1' => {
                                        ':////accu?merge_ids={bam_file}' => {'merge_ids' => '#merge_id#'},
                                      },
                            },
    },

    {
      -logic_name        => 'MergeBam',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::MergeBam',
      -max_retry_count   => 1,
      -parameters        => {
                              merge_ids      => '#merge_ids#',
                              work_directory => catdir($dir, '#species#'),
                              samtools_dir   => $self->o('samtools_dir'),
                              vcf            => $self->o('vcf'),
                              use_csi        => $self->o('use_csi'),
                              clean_up       => $self->o('clean_up'),
                              bigwig         => $self->o('bigwig'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '3' => ['WriteIniFile'],
                              '4' => ['CreateBigWig'],
                            },
    },

    {
      -logic_name        => 'CreateBigWig',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::CreateBigWig',
      -can_be_empty      => 1,
      -max_retry_count   => 1,
      -parameters        => {
                              bedtools_dir  => $self->o('bedtools_dir'),
                              ucscutils_dir => $self->o('ucscutils_dir'),
                              length_file   => '#genome_file#'.'.lengths.txt',
                              clean_up      => $self->o('clean_up'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['WriteIniFile'],
    },

    {
      -logic_name        => 'WriteIniFile',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteIniFile',
      -max_retry_count   => 1,
      -parameters        => {
                              work_directory => catdir($dir, '#species#'),
                              run            => $self->o('run'),
                              study          => $self->o('study'),
                              merge_level    => $self->o('merge_level'),
                              ini_type       => $self->o('ini_type'),
                              bigwig         => $self->o('bigwig'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['EmailReport'],
    },

    {
      -logic_name        => 'EmailReport',
      -module            => 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::EmailReport',
      -max_retry_count   => 1,
      -parameters        => {
                              email        => $self->o('email'),
                              subject      => 'Short Read Alignment pipeline: Report for #species#',
                              samtools_dir => $self->o('samtools_dir'),
                            },
      -rc_name           => 'normal',
    },
  ];
}

sub modify_analyses {
  my ($self, $analyses) = @_;
}

sub resource_classes {
  my ($self) = @_;
  
  my $threads                   = $self->o('threads');
  my $aligner                   = $self->o('aligner');
  my $index_memory_default      = $self->o('index_memory_default');
  my $index_memory_high_default = $self->o('index_memory_high_default');
  my $align_memory_default      = $self->o('align_memory_default');
  my $align_memory_high_default = $self->o('align_memory_high_default');
  
  my $index_mem   = $self->o('index_memory')      || $$index_memory_default{$aligner};
  my $index_himem = $self->o('index_memory_high') || $$index_memory_high_default{$aligner};
  my $align_mem   = $self->o('align_memory')      || $$align_memory_default{$aligner};
  my $align_himem = $self->o('align_memory_high') || $$align_memory_high_default{$aligner};
  
  return {
    %{$self->SUPER::resource_classes},
    'index_default' => {'LSF' => '-q production-rh6 -n '. ($threads + 1) .' -M '.$index_mem.' -R "rusage[mem='.$index_mem.',tmp=16000] span[hosts=1]"'},
    'index_himem'   => {'LSF' => '-q production-rh6 -n '. ($threads + 1) .' -M '.$index_himem.' -R "rusage[mem='.$index_himem.',tmp=16000] span[hosts=1]"'},
    'align_default' => {'LSF' => '-q production-rh6 -n '. ($threads + 1) .' -M '.$align_mem.' -R "rusage[mem='.$align_mem.',tmp=16000] span[hosts=1]"'},
    'align_himem'   => {'LSF' => '-q production-rh6 -n '. ($threads + 1) .' -M '.$align_himem.' -R "rusage[mem='.$align_himem.',tmp=16000] span[hosts=1]"'},
  }
}

1;
