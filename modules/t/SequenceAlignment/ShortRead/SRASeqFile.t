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

use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);
use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SRASeqFile;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $single_read_run = 'SRR1927173';
my $paired_read_run = 'SRR2094641';

my $tmp_dir = "/tmp/$ENV{'USER'}/sraseqfile";

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SRASeqFile';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(retrieve_files);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj  = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SRASeqFile->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);



{
make_path($tmp_dir);

$obj->param('work_dir', $tmp_dir);
$obj->param('run_id', $single_read_run);

$obj->run();

my @output = @{$obj->param('output')};

is(scalar(@output), 1, "run method: one set of output parameters ($single_read_run)");

my $output      = $output[0];
my $read_file_1 = $$output{'seq_file_1'};
my $read_file_2 = $$output{'seq_file_2'};
my $sam_file    = $$output{'sam_file'};

is(-e $read_file_1, 1, "run method: single read file exists ($single_read_run)");
ok(!defined($read_file_2), "run method: single read file does not exist ($single_read_run)");

is(-s $read_file_1, 1975556812, "run method: single read file size ($single_read_run)");

is($sam_file, catdir($tmp_dir, "$single_read_run.sam"), "run method: sam file name ($single_read_run)");

remove_tree($tmp_dir);
}

{
make_path($tmp_dir);

$obj->param('work_dir', $tmp_dir);
$obj->param('run_id', $paired_read_run);

$obj->run();

my @output = @{$obj->param('output')};

is(scalar(@output), 1, "run method: one set of output parameters ($paired_read_run)");

my $output      = $output[0];
my $read_file_1 = $$output{'seq_file_1'};
my $read_file_2 = $$output{'seq_file_2'};
my $sam_file    = $$output{'sam_file'};

is(-e $read_file_1, 1, "run method: paired read file exists ($paired_read_run)");
is(-e $read_file_2, 1, "run method: paired read file exists ($paired_read_run)");

is(-s $read_file_1, 3305191542, "run method: paired read file size ($paired_read_run)");
is(-s $read_file_2, 3274992152, "run method: paired read file size ($paired_read_run)");

is($sam_file, catdir($tmp_dir, "$paired_read_run.sam"), "run method: sam file name ($paired_read_run)");

remove_tree($tmp_dir);
}


done_testing();
