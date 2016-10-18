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

use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::LoadGFF3::DeleteGenes;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $module_name  = 'Bio::EnsEMBL::EGPipeline::LoadGFF3::DeleteGenes';
my @hive_methods = qw(param_defaults fetch_input run write_output);
can_ok($module_name, @hive_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $dg_obj  = Bio::EnsEMBL::EGPipeline::LoadGFF3::DeleteGenes->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$dg_obj->input_job($job_obj);

# These are the modules param_defaults.
$dg_obj->param('db_type', 'core');

# Species is a mandatory parameter.
$dg_obj->param('species', $species);

my $logic_name = 'vectorbase_anopheles_gambiae';

my $aa = $dba->get_adaptor('Analysis');
my $analysis = $aa->fetch_by_logic_name($logic_name);

is_rows(77, $dba, 'gene', 'where analysis_id = ? ', [$analysis->dbID]);

# Deleting a gene affects a _lot_ of tables...
$testdb->hide($dbtype, 
  qw(dna_align_feature exon exon_transcript external_synonym gene gene_attrib 
     identity_xref meta_coord object_xref ontology_xref protein_feature 
     supporting_feature transcript transcript_attrib 
     transcript_supporting_feature translation translation_attrib xref));

$dg_obj->param('logic_name', $logic_name);
$dg_obj->run($dba);
is_rows(0, $dba, 'gene', 'where analysis_id = ? ', [$analysis->dbID]);

$testdb->restore($dbtype, 
  qw(dna_align_feature exon exon_transcript external_synonym gene gene_attrib 
     identity_xref meta_coord object_xref ontology_xref protein_feature 
     supporting_feature transcript transcript_attrib 
     transcript_supporting_feature translation translation_attrib xref));

done_testing();
