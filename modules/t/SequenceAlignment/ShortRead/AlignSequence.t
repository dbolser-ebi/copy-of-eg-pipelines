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
my %aligner_classes =
(
  bowtie2 => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::Bowtie2Aligner',
  bwa     => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::BwaAligner',
  gsnap   => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::GsnapAligner',
  hisat2  => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::HISAT2Aligner',
  star    => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::StarAligner',
  tophat2 => 'Bio::EnsEMBL::EGPipeline::Common::Aligner::TopHat2Aligner',
);

my %aligner_dirs =
(
  bowtie2 => '/nfs/panda/ensemblgenomes/external/bowtie2',
  bwa     => '/nfs/panda/ensemblgenomes/external/bwa',
  gsnap   => '/nfs/panda/ensemblgenomes/external/gmap-gsnap/bin',
  hisat2  => '/nfs/panda/ensemblgenomes/external/hisat2',
  star    => '/nfs/panda/ensemblgenomes/external/STAR',
  tophat2 => '/nfs/panda/ensemblgenomes/external/tophat2',
);

my %aligner_versions =
(
  bowtie2 => '2.2.6',
  bwa     => '0.6.1-r104',
  gsnap   => '2012-11-09',
  hisat2  => '2.0.4',
  star    => 'unknown',
  tophat2 => '2.1.1',
);

$obj->param(samtools_dir => '/nfs/panda/ensemblgenomes/external/samtools');

foreach my $aligner (sort keys %aligner_classes) {
  my $class   = $aligner_classes{$aligner};
  my $dir     = $aligner_dirs{$aligner};
  my $version = $aligner_versions{$aligner};
  my $intron_length = $aligner =~ /^(hisat2|star|tophat2)$/ ? 23 : undef;
  
  $obj->param('aligner_class', $class);
  $obj->param('aligner_dir',   $dir);
  
  $obj->fetch_input();
  
  my $aligner_object = $obj->param('aligner_object');
  isa_ok($aligner_object, $class, 'fetch_input method');
  is($aligner_object->version,             $version,       "fetch_input method: $aligner version correct");
  is($aligner_object->{max_intron_length}, $intron_length, "fetch_input method: $aligner max_intron_length correct");
  is($aligner_object->{threads},           4,              "fetch_input method: $aligner threads correct");
}

done_testing();
