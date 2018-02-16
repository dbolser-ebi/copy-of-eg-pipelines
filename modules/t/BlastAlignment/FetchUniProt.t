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

use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot;

my $tmp_dir = "/tmp/$ENV{'USER'}/fetchuniprot";
make_path($tmp_dir);

my $module_name    = 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(convert_to_fasta);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj  = Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchUniprot->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('ebi_path'),     '/ebi/ftp/pub/databases/uniprot/current_release/knowledgebase', 'param_defaults method: ebi_path');
is($obj->param('ftp_uri'),      'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase', 'param_defaults method: ftp_uri');
is($obj->param('uniprot_db'),   'sprot', 'param_defaults method: uniprot_db');
is($obj->param('file_varname'), 'uniprot_fasta_file', 'param_defaults method: file_varname');

# Set some params that would otherwise come via the pipeline.
$obj->param('out_dir',         $tmp_dir);

{
$obj->fetch_input();
is($obj->param('uniprot_file'), 'uniprot_sprot.fasta.gz', 'fetch_input method: uniprot_file');
is($obj->param('output_file'),  catdir($tmp_dir, 'uniprot_sprot.fasta'), 'fetch_input method: output_file');

$obj->run();
ok(-e $obj->param('output_file'),  'run method: output_file exists');

my $output = path($obj->param('output_file'));
my @lines  = $output->lines( { chomp => 1, count => 1 } );
my ($id) = $lines[0] =~ /^>(\S+)/;
ok($id !~ /\|/, 'run method: sequence ID format');
ok($id !~ /^sp/, 'run method: sequence ID format');
}

{
$obj->param('taxonomic_level', 'invertebrates');

$obj->fetch_input();
is($obj->param('uniprot_file'), 'uniprot_sprot_invertebrates.dat.gz', 'fetch_input method: uniprot_file');
is($obj->param('output_file'),  catdir($tmp_dir, 'sprot_invertebrates.fa'), 'fetch_input method: output_file');

$obj->run();
ok(-e $obj->param('output_file'),  'run method: output_file exists');

my $output = path($obj->param('output_file'));
my @lines  = $output->lines( { chomp => 1, count => 1 } );
my ($id) = $lines[0] =~ /^>(\S+)/;
ok($id !~ /\|/, 'run method: sequence ID format');
ok($id !~ /^sp/, 'run method: sequence ID format');
}

remove_tree($tmp_dir);

done_testing();
