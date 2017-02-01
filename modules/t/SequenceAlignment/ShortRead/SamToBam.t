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
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SamToBam;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $test_files_dir = catdir($FindBin::Bin, '../../test-files/ShortReadAlignment');
my $test_sam_file  = catdir($test_files_dir, 'agam_single_read.sam');

my $tmp_dir = "/tmp/$ENV{'USER'}/samtobam";

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SamToBam';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
can_ok($module_name, @hive_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj  = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SamToBam->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('threads'),  4,    'param_defaults method: threads');
is($obj->param('memory'),   8000, 'param_defaults method: memory');
is($obj->param('clean_up'), 1,    'param_defaults method: clean_up');

$obj->param('samtools_dir', '/nfs/panda/ensemblgenomes/external/samtools');

$obj->fetch_input();

my $aligner_object = $obj->param('aligner_object');
my $class = 'Bio::EnsEMBL::EGPipeline::Common::Aligner';
isa_ok($aligner_object, $class, 'fetch_input method');

make_path($tmp_dir);

my $test_sam_link = catdir($tmp_dir, 'agam_single_read.sam');

symlink $test_sam_file, $test_sam_link;

# Single read file
$obj->param('sam_file', $test_sam_link);
  
$obj->run();
  
is(-e $obj->param('bam_file'), 1, "run method: bam_file exists");

remove_tree($tmp_dir);

done_testing();
