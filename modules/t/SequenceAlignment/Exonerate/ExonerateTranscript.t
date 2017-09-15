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

use Cwd;
use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);
use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::ExonerateTranscript;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../../';

my $cwd = cwd();

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $aa = $dba->get_adaptor('Analysis');

my $tmp_dir = "/tmp/$ENV{'USER'}/exoneratetranscript";
make_path($tmp_dir);

my $test_files_dir = catdir($FindBin::Bin, '../../test-files/Exonerate');
my $genome_file    = catdir($tmp_dir, 'agam_genome.fa');
my $cdna_file      = catdir($tmp_dir, 'agam_cdna.fa');

symlink(catdir($test_files_dir, 'agam_genome.fa'), $genome_file);
symlink(catdir($test_files_dir, 'agam_cdna.fa'), $cdna_file);

my $exonerate_exe    = 'exonerate';
my $exonerate_params = '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --bestn 10 --score 1000 --maxintron 25000 ';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::ExonerateTranscript';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(fetch_runnable filter_output results_by_index make_genes);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj  = Bio::EnsEMBL::EGPipeline::SequenceAlignment::Exonerate::ExonerateTranscript->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('coverage'),       90, 'param_defaults method: coverage');
is($obj->param('percent_id'),     97, 'param_defaults method: percent_id');
is($obj->param('best_in_genome'), 1,  'param_defaults method: best_in_genome');

# Set some params that would otherwise come via the pipeline.
$obj->param('workdir',     $tmp_dir);
$obj->param('species',     $species);
$obj->param('queryfile',   $genome_file);
$obj->param('seq_file',    $cdna_file);
$obj->param('server_file', '');
$obj->param('seq_type',    'dna');
$obj->param('biotype',     'cdna');
$obj->param('logic_name',  'cdna_e2g');

{

# create analysis in core db.
my $analysis = Bio::EnsEMBL::Analysis->new(
  -logic_name   => $obj->param('logic_name'),
  -program_file => $exonerate_exe,
  -parameters   => $exonerate_params,
  -module       => $module_name,
);
my $analysis_id = $aa->store($analysis);

$obj->fetch_input();
$obj->run();

is_rows(33, $dba, 'dna_align_feature', 'WHERE analysis_id = ?', [$analysis_id]);
is_rows( 9, $dba, 'gene',              'WHERE analysis_id = ?', [$analysis_id]);
is_rows( 9, $dba, 'transcript',        'WHERE analysis_id = ?', [$analysis_id]);

}

chdir($cwd);
remove_tree($tmp_dir);

done_testing();
