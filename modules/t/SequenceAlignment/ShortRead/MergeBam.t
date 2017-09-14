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

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use File::Compare qw(compare);
use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);
use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::MergeBam;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $test_files_dir  = catdir($FindBin::Bin, '../../test-files/ShortReadAlignment');
my $genome_file     = catdir($test_files_dir, 'agam_genome.fa');
my $bam_1_file      = catdir($test_files_dir, 'agam_single_read.sorted.bam');
my $bam_2_file      = catdir($test_files_dir, 'agam_paired_read.sorted.bam');
my $merged_bam_file = catdir($test_files_dir, 'agam_merged.sorted.bam');
my $bam_index_file  = catdir($test_files_dir, 'agam_merged.sorted.bai');
my $bam_cindex_file = catdir($test_files_dir, 'agam_merged.sorted.csi');
my $vcf_file        = catdir($test_files_dir, 'agam_merged.sorted.vcf');

my $samtools_dir = undef;
my $bcftools_dir = undef;

my $tmp_dir = "/tmp/$ENV{'USER'}/mergebam";

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::MergeBam';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(merge_bam bam_files_from_db);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::MergeBam->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('samtools_dir'), $samtools_dir, 'param_defaults method: samtools_dir');
is($obj->param('bcftools_dir'), $bcftools_dir, 'param_defaults method: bcftools_dir');
is($obj->param('threads'),      1,             'param_defaults method: threads');
is($obj->param('clean_up'),     1,             'param_defaults method: clean_up');
is($obj->param('use_csi'),      0,             'param_defaults method: use_csi');
is($obj->param('vcf'),          0,             'param_defaults method: vcf');

# Set some params that would otherwise come via the pipeline.
$obj->param('species', $species);
$obj->param('results_dir', $tmp_dir);
$obj->param('fasta_file', $genome_file);

$obj->fetch_input();

my $bam_1_link = catdir($tmp_dir, 'agam_single_read.sorted.bam');
my $bam_2_link = catdir($tmp_dir, 'agam_paired_read.sorted.bam');

my $assembly = $obj->param('assembly');
is($assembly, 'AgamP4', 'fetch_input method: assembly');

# Run with standard parameters
{
make_path($tmp_dir);

symlink $bam_1_file, $bam_1_link;
symlink $bam_2_file, $bam_2_link;

my $bam_files = [$bam_1_link, $bam_2_link];
my $tmp_merged_bam_file = catdir($tmp_dir, "test_merge_$assembly.bam");
my $tmp_bam_index_file  = catdir($tmp_dir, "test_merge_$assembly.bam.bai");

my $cmds = $obj->merge_bam($bam_files, $tmp_merged_bam_file);

ok(-s $tmp_merged_bam_file > 0, 'merge_bam method: bam file exists');
ok(-s $tmp_bam_index_file  > 0, 'merge_bam method: bam index exists');
ok(! -e $bam_1_link, 'merge_bam method: source bam deleted');
ok(! -e $bam_2_link, 'merge_bam method: source bam deleted');

remove_tree($tmp_dir);
}

# Leave component bam files alone; alternative index type; create vcf file.
{
make_path($tmp_dir);

$obj->param('clean_up', 0);
$obj->param('use_csi', 1);
$obj->param('vcf', 1);

symlink $bam_1_file, $bam_1_link;
symlink $bam_2_file, $bam_2_link;

my $bam_files = [$bam_1_link, $bam_2_link];
my $tmp_merged_bam_file = catdir($tmp_dir, "test_merge_$assembly.bam");
my $tmp_bam_cindex_file = catdir($tmp_dir, "test_merge_$assembly.bam.csi");
my $tmp_vcf_file        = catdir($tmp_dir, "test_merge_$assembly.vcf");

my $cmds = $obj->merge_bam($bam_files, $tmp_merged_bam_file);

ok(-s $tmp_merged_bam_file > 0, 'merge_bam method: bam file exists');
ok(-s $tmp_bam_cindex_file > 0, 'merge_bam method: bam index exists');
ok(-s $tmp_vcf_file        > 0, 'merge_bam method: vcf file exists');
is(compare($vcf_file, $tmp_vcf_file), 0, 'merge_bam method: vcf contents correct');
ok(-e $bam_1_link, 'merge_bam method: source bam not deleted');
ok(-e $bam_2_link, 'merge_bam method: source bam not deleted');

remove_tree($tmp_dir);
}

done_testing();
