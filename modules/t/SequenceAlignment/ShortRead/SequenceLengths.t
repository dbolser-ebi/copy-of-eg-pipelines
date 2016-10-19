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
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceLengths;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $test_files_dir = catdir($FindBin::Bin, '../../test-files');
my $genome_file    = catdir($test_files_dir, 'agam_genome.fa');
my $length_file   = catdir($test_files_dir, 'agam_length.txt');

my $tmp_dir = "/tmp/$ENV{'USER'}/sequencelength";
my $tmp_length_file = catdir($tmp_dir, 'agam_length.txt');

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceLengths';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
can_ok($module_name, @hive_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj = Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceLengths->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

$obj->param('fasta_file', $genome_file);
$obj->param('length_file', $tmp_length_file);

make_path($tmp_dir);

$obj->run();

is(-e $tmp_length_file, 1, "run method: length file exists");

is(compare($length_file, $tmp_length_file), 0, "run method: lengths correct");

remove_tree($tmp_dir);

done_testing();
