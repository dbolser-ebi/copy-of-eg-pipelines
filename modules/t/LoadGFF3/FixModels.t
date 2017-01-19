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
use Bio::EnsEMBL::EGPipeline::LoadGFF3::FixModels;
use Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);
my $ta     = $dba->get_adaptor('Transcript');

my $gff3_file    = $FindBin::Bin.'/../test-files/LoadGFF3/agam_seqedits.gff3';
my $fasta_file   = $FindBin::Bin.'/../test-files/LoadGFF3/agam.fa';
my $protein_file = $FindBin::Bin.'/../test-files/LoadGFF3/agam_protein.fa';

my $gene_source = 'Ensembl';
my $logic_name  = 'gff3_test';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::LoadGFF3::FixModels';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(
  fix_models shift_translation shift_gene update_gene_position
  load_fasta set_protein_coding update_translation_start);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

# Need to load some genes which need seq edits
my $lg_obj  = Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3->new;
my $lg_job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$lg_obj->input_job($lg_job_obj);

$lg_obj->param('db_type',        'core');
$lg_obj->param('species',        $species);
$lg_obj->param('gff3_file',      $gff3_file);
$lg_obj->param('fasta_file',     $fasta_file);
$lg_obj->param('logic_name',     $logic_name);
$lg_obj->param('gene_source',    $gene_source);
$lg_obj->param('gene_types',     ['gene']);
$lg_obj->param('mrna_types',     ['mRNA']);
$lg_obj->param('exon_types',     ['exon']);
$lg_obj->param('cds_types',      ['CDS']);
$lg_obj->param('utr_types',      ['five_prime_UTR', 'three_prime_UTR']);
$lg_obj->param('ignore_types',   ['region', 'chromosome']);
$lg_obj->param('types_complete', 1);
$lg_obj->param('nontranslating', 'nontranslating_CDS');
$lg_obj->param('polypeptides',   0);
$lg_obj->param('prediction',     0);
$lg_obj->param('use_name_field', 'no');
$lg_obj->param('stable_ids',     {});

$lg_obj->run();

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $fm_obj = Bio::EnsEMBL::EGPipeline::LoadGFF3::FixModels->new;
my $fm_job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$fm_obj->input_job($fm_job_obj);

# Set and check default parameters.
my $param_defaults = $fm_obj->param_defaults();
$fm_obj->input_job->param_init($param_defaults);
is($fm_obj->param('db_type'), 'core', 'param_defaults method: db_type');

# Mandatory parameters.
$fm_obj->param('species',            $species);
$fm_obj->param('logic_name',         $logic_name);
$fm_obj->param('protein_fasta_file', $protein_file);

# run method
{
my $transcript_1 = $ta->fetch_by_stable_id('test_protein_shift_start-RA');

my $aa_seq_1_before = '*KKV*HF*PVSLSASKSLALKQKFENFLTLVMLIHRITCDIMHYICLIMIIQTDYDLH*KNVPCQICSFR*RYFLPKNHTFVPRCS*ALIFC*HQRCV*NLCNRLMTIPRVICLF*SRRNTYSITKVVAGT*TKFFIILS*RKNQCTII*PL*SCAIQFVKR*WQMVKVLLLVYGLK*SCAITKFTWANGLSIILNKILPFVGLMFQS*RGRNVARSFRKIK*SYQNLIAE*LKHNYVLRTRKIVP*LSFVNLVRQDHYL*PSVTQFTW*GCLQFILMIVMCKLRCLIRFHHFWIGLKQSYGPILNLC';
my $aa_seq_1_after  = 'KKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPL';

is($transcript_1->translate->seq(), $aa_seq_1_before, 'run method (before): amino acid sequence as expected');
is($transcript_1->translation->start(), 95, 'run method (before): sequence start correct');

$fm_obj->run();

# Need to re-fetch from the database.
$transcript_1 = $ta->fetch_by_stable_id('test_protein_shift_start-RA');

is($transcript_1->translate->seq(), $aa_seq_1_after, 'run method (after): amino acid sequence as expected');

my $transcript_1_seq_edits = $transcript_1->get_all_SeqEdits();
is(scalar(@$transcript_1_seq_edits), 0, 'run method: no trancript-level seq edits');

my $translation_1_seq_edits = $transcript_1->translation->get_all_SeqEdits();
is(scalar(@$translation_1_seq_edits), 0, 'run method: no translation-level seq edits');

is($transcript_1->translation->start, 97, 'run method (after): sequence start shifted correctly');
is($transcript_1->biotype, 'protein_coding', 'run method (after): transcript biotype correct');
}

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

done_testing();
