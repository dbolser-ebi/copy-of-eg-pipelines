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
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteIniFile;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $test_files_dir = catdir($FindBin::Bin, '../../test-files');
my $bam_file       = catdir($test_files_dir, 'agam_single_read.sorted.bam');
my $bigwig_file    = catdir($test_files_dir, 'agam_single_read.sorted.bigwig');
my $ini_1_file     = catdir($test_files_dir, 'agam_single_read.1.ini');
my $ini_2_file     = catdir($test_files_dir, 'agam_single_read.2.ini');
my $ini_3_file     = catdir($test_files_dir, 'agam_single_read.3.ini');

my $tmp_dir        = "/tmp/$ENV{'USER'}/writeinifile";
my $tmp_ini_1_file = catdir($tmp_dir, 'test_1.ini');
my $tmp_ini_2_file = catdir($tmp_dir, 'test_2.ini');
my $tmp_ini_3_file = catdir($tmp_dir, 'test_3.ini');

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteIniFile';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(sra_desc ena_link);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteIniFile->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('bigwig'),   0,              'param_defaults method: bigwig');
is($obj->param('ini_type'), 'rnaseq_align', 'param_defaults method: ini_type');

# Set some params that would otherwise come via the pipeline.
$obj->param('species', $species);
$obj->param('results_dir', $tmp_dir);

# bam file, no ENA links
{
make_path($tmp_dir);

$obj->param('merged_bam_file', $bam_file);
$obj->param('merge_level',     'taxon');
$obj->param('merge_id',        'test_1');
$obj->param('merge_label',     'Anopheles gambiae');
$obj->param('run_ids',         '');

$obj->run();

is($obj->param('ini_file'), $tmp_ini_1_file, "run method: ini filename");

is(-e $tmp_ini_1_file, 1, "run method: ini file exists");

is(compare($ini_1_file, $tmp_ini_1_file), 0, "run method: file contents correct");

remove_tree($tmp_dir);
}

# bigwig file, ENA link
{
make_path($tmp_dir);

$obj->param('bigwig',      1);
$obj->param('bw_file',     $bigwig_file);
$obj->param('merge_level', 'run');
$obj->param('merge_id',    'test_2');
$obj->param('merge_label', 'SRR1927173');
$obj->param('run_ids',     'SRR1927173');

$obj->run();

is($obj->param('ini_file'), $tmp_ini_2_file, "run method: ini filename");

is(-e $tmp_ini_2_file, 1, "run method: ini file exists");

is(compare($ini_2_file, $tmp_ini_2_file), 0, "run method: file contents correct");

remove_tree($tmp_dir);
}

# bigwig file, multiple ENA links
{
make_path($tmp_dir);

$obj->param('bigwig',      1);
$obj->param('bw_file',     $bigwig_file);
$obj->param('merge_level', 'study');
$obj->param('merge_id',    'test_3');
$obj->param('merge_label', 'SRP056482');
$obj->param('run_ids',     'SRR1927173,SRR1927174');

$obj->run();

is($obj->param('ini_file'), $tmp_ini_3_file, "run method: ini filename");

is(-e $tmp_ini_3_file, 1, "run method: ini file exists");

is(compare($ini_3_file, $tmp_ini_3_file), 0, "run method: file contents correct");

remove_tree($tmp_dir);
}

done_testing();
