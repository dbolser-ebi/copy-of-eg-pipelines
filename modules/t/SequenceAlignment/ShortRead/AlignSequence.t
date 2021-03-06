# Copyright [2016] EMBL-European Bioinformatics Institute
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

# IMPORTANT: You need to run this test with at
# least 8GB of memory and 8GB of /tmp space.

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);
use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::AlignSequence;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $test_files_dir     = catdir($FindBin::Bin, '../../test-files/ShortReadAlignment');
my $genome_file        = catdir($test_files_dir, 'agam_genome.fa');
my $single_read_file   = catdir($test_files_dir, 'agam_single_read.fastq');
my $paired_read_1_file = catdir($test_files_dir, 'agam_paired_read_1.fastq');
my $paired_read_2_file = catdir($test_files_dir, 'agam_paired_read_2.fastq');

my $tmp_dir = "/tmp/$ENV{'USER'}/alignsequence";

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::AlignSequence';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(max_intron_length);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj  = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::AlignSequence->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('threads'),  4,         'param_defaults method: threads');
is($obj->param('run_mode'), 'default', 'param_defaults method: run_mode');
is($obj->param('gtf_file'), undef,     'param_defaults method: gtf_file');
is($obj->param('clean_up'), 1,         'param_defaults method: clean_up');

# Set some params that would otherwise come via the pipeline.
$obj->param('species', $species);
$obj->param('escape_branch', undef);
$obj->param('max_intron', 1);

# Test basic functionality for each aligner.
my %aligners =
(
  bowtie2 => {
               class         => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::Bowtie2Aligner',
               dir           => undef,
               version       => '2.3.2',
               single_total  => 25000,
               single_mapped => 1625,
               paired_total  => 50000,
               paired_mapped => 3746,
             },
  bwa     => {
               class         => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::BwaAligner',
               dir           => undef,
               version       => '0.7.16a-r1181',
               single_total  => 25000,
               single_mapped => 2137,
               paired_total  => 50000,
               paired_mapped => 3705,
             },
  gsnap   => {
               class         => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::GsnapAligner',
               dir           => undef,
               version       => '2017-01-14',
               single_total  => 25967,
               single_mapped => 2382,
               paired_total  => 50030,
               paired_mapped => 3921,
             },
  hisat2  => {
               class         => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::HISAT2Aligner',
               dir           => undef,
               version       => '2.0.5',
               single_total  => 25906,
               single_mapped => 2250,
               paired_total  => 50029,
               paired_mapped => 3734,
             },
  star    => {
               class         => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::StarAligner',
               dir           => '/nfs/panda/ensemblgenomes/external/STAR',
               version       => '2.3.1z',
               single_total  => 8630,
               single_mapped => 8630,
               paired_total  => 3688,
               paired_mapped => 3688,
             },
  tophat2 => {
               class         => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::TopHat2Aligner',
               dir           => undef,
               version       => '2.1.1',
               single_total  => 3158,
               single_mapped => 3158,
               paired_total  => 3626,
               paired_mapped => 3626,
             },
);

$obj->param('samtools_dir', undef);

foreach my $aligner (sort keys %aligners) {
  my $class         = $aligners{$aligner}{class};
  my $dir           = $aligners{$aligner}{dir};
  my $version       = $aligners{$aligner}{version};
  my $single_total  = $aligners{$aligner}{single_total};
  my $single_mapped = $aligners{$aligner}{single_mapped};
  my $paired_total  = $aligners{$aligner}{paired_total};
  my $paired_mapped = $aligners{$aligner}{paired_mapped};
  my $intron_length = $aligner =~ /^(hisat2|star|tophat2)$/ ? 25000 : undef;
  
  $obj->param('aligner_class', $class);
  $obj->param('aligner_dir',   $dir);
  
  $obj->fetch_input();
  
  my $aligner_object = $obj->param('aligner_object');
  isa_ok($aligner_object, $class, 'fetch_input method');
  is($aligner_object->version,             $version,       "fetch_input method: $aligner version correct");
  is($aligner_object->{max_intron_length}, $intron_length, "fetch_input method: $aligner max_intron_length correct");
  is($aligner_object->{threads},           4,              "fetch_input method: $aligner threads correct");
  
  # Have per-aligner directories for indexes and results
  my $results_dir = catdir($tmp_dir, $aligner);
  make_path($results_dir);

  my $genome_link        = catdir($results_dir, 'agam_genome.fa');
  my $single_read_link   = catdir($results_dir, 'agam_single_read.fastq');
  my $paired_read_1_link = catdir($results_dir, 'agam_paired_read_1.fastq');
  my $paired_read_2_link = catdir($results_dir, 'agam_paired_read_2.fastq');
  
  symlink $genome_file, $genome_link;
  symlink $single_read_file, $single_read_link;
  symlink $paired_read_1_file, $paired_read_1_link;
  symlink $paired_read_2_file, $paired_read_2_link;
  
  $obj->param('genome_file', $genome_link);
  
  unless ($aligner_object->index_exists($genome_link)) {
    $aligner_object->index_file($genome_link);
  }
  
  # Single read file
  $obj->param('seq_file_1', $single_read_link);
  $obj->param('seq_file_2', undef);
  
  $obj->run();
  
  is(-e $obj->param('sam_file'), 1, "run method: $aligner sam_file exists (single read)");
  
  # More than just checking that a file is produced,
  # check the results by converting to BAM, then pull out stats.
  $aligner_object->sam_to_bam($obj->param('sam_file'));
  (my $single_bam_file = $obj->param('sam_file')) =~ s/\.sam$/\.bam/;
  
  my $single_stats = $aligner_object->get_bam_stats($single_bam_file);
  my ($total_single_reads)  = $single_stats =~ /^(\d+).*in total/m;
  my ($mapped_single_reads) = $single_stats =~ /^(\d+).*mapped/m;
  
  is($total_single_reads,  $single_total,  "run method: $aligner number of reads (single read)");
  is($mapped_single_reads, $single_mapped, "run method: $aligner mapped reads (single read)");

  # Paired read files
  $obj->param('seq_file_1', $paired_read_1_link);
  $obj->param('seq_file_2', $paired_read_2_link);
  
  $obj->run();
  
  is(-e $obj->param('sam_file'), 1, "run method: $aligner sam_file exists (paired read)");

  # More than just checking that a file is produced,
  # check the results by converting to BAM, then pull out stats.
  $aligner_object->sam_to_bam($obj->param('sam_file'));
  (my $paired_bam_file = $obj->param('sam_file')) =~ s/\.sam$/\.bam/;
  
  my $paired_stats = $aligner_object->get_bam_stats($paired_bam_file);
  my ($total_paired_reads)  = $paired_stats =~ /^(\d+).*in total/m;
  my ($mapped_paired_reads) = $paired_stats =~ /^(\d+).*mapped/m;
  
  is($total_paired_reads,  $paired_total,  "run method: $aligner number of reads (paired read)");
  is($mapped_paired_reads, $paired_mapped, "run method: $aligner mapped reads (paired read)");
  
  remove_tree($results_dir);

  # STAR puts files in the current dir, and there is no way
  # to redirect them elsewhere...
  if ($aligner eq 'star') {
    unlink 'Log.final.out';
    unlink 'Log.out';
    unlink 'Log.progress.out';
    unlink 'Log.std.out';
    unlink 'SJ.out.tab';
  }
}

done_testing();
