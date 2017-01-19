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
use File::Spec::Functions qw(catdir);
use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::CreateBigWig;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $test_files_dir  = catdir($FindBin::Bin, '../../test-files/ShortReadAlignment');
my $length_file     = catdir($test_files_dir, 'agam_length.txt');
my $bam_file        = catdir($test_files_dir, 'agam_single_read.sorted.bam');
my $bigwig_file     = catdir($test_files_dir, 'agam_single_read.sorted.bigwig');
my $tmp_bigwig_file = catdir($test_files_dir, 'agam_single_read.sorted.bw');

my $bedtools_dir  = '/nfs/panda/ensemblgenomes/external/bedtools/bin';
my $ucscutils_dir = '/nfs/panda/ensemblgenomes/external/ucsc_utils';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::CreateBigWig';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
can_ok($module_name, @hive_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::CreateBigWig->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('bedtools_dir'),  $bedtools_dir,  'param_defaults method: bedtools_dir');
is($obj->param('ucscutils_dir'), $ucscutils_dir, 'param_defaults method: ucscutils_dir');
is($obj->param('clean_up'),      1,              'param_defaults method: clean_up');

# Set some params that would otherwise come via the pipeline.
$obj->param('length_file', $length_file);
$obj->param('merged_bam_file', $bam_file);

$obj->run();

is($obj->param('bw_file'), $tmp_bigwig_file, "run method: bigwig filename");

is(-e $tmp_bigwig_file, 1, "run method: bigwig file exists");

is(compare($bigwig_file, $tmp_bigwig_file), 0, "run method: file contents correct");

my $cmds =
  "$bedtools_dir/bedtools genomecov ".
  " -g $length_file ".
  " -ibam $bam_file ".
  " -bg ".
  " -split ".
  " > $bam_file.wig ; ".
  "$ucscutils_dir/wigToBigWig ".
  " $bam_file.wig ".
  " $length_file ".
  " $tmp_bigwig_file ";
is($obj->param('cmds'), $cmds, "run method: bigwig command line");

unlink($tmp_bigwig_file);

done_testing();
