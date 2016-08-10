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

use FindBin;
use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);
my $ga     = $dba->get_adaptor('Gene');
my $pta    = $dba->get_adaptor('PredictionTranscript');

my $gff3_file        = $FindBin::Bin.'/../test-files/agam.gff3';
my $ncrna_gff3_file  = $FindBin::Bin.'/../test-files/agam_ncrna.gff3';
my $pseudo_gff3_file = $FindBin::Bin.'/../test-files/agam_pseudo.gff3';
my $transl_gff3_file = $FindBin::Bin.'/../test-files/agam_transl.gff3';
my $fasta_file       = $FindBin::Bin.'/../test-files/agam.fa';

my $gene_source = 'Ensembl';
my $logic_name  = 'gff3_test';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(
  add_exons add_predicted_transcript add_pseudogenic_transcript add_transcript 
  add_transcripts add_translation check_db check_seq_ids correct_exon_overlaps 
  get_cds get_cds_id get_polypeptide get_utr infer_exons infer_transcript 
  infer_translation load_db load_genes new_exon new_gene new_predicted_exon 
  new_predicted_transcript new_transcript new_translation set_exon_phase 
  set_nontranslating_gene set_nontranslating_transcript translation_coordinates 
);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $lg_obj  = Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$lg_obj->input_job($job_obj);

# Analysis against which genes will be stored.
my $aa       = $dba->get_adaptor('Analysis');
my $analysis = $aa->fetch_by_logic_name($logic_name);

# These are the modules param_defaults.
$lg_obj->param('db_type', 'core');

# Mandatory parameters.
$lg_obj->param('species', $species);
$lg_obj->param('gff3_file', $gff3_file);
$lg_obj->param('fasta_file', $fasta_file);
$lg_obj->param('logic_name', $logic_name);
$lg_obj->param('gene_source', $gene_source);
$lg_obj->param('nontranslating', 'nontranslating_CDS');
$lg_obj->param('polypeptides', 0);
$lg_obj->param('prediction', 0);
$lg_obj->param('use_name_field', 'no');
$lg_obj->param('stable_ids', {});

# There's a lot to test here, so this file is broken into chapters.
# Chapter 1: basic functionality tested with a single simple gene,
#            nothing, quirky, standard GFF3.
# Chapter 2: loading names instead of IDs, or as xrefs.
# Chapter 3: generating predicted transcripts rather than genes.
# Chapter 4: ncRNA genes (including tests for genes on -ve strand).
# Chapter 5: pseudogenes, which can be specified in many formats.
# Chapter 6: translation scenarios, e.g. no exons, odd phases,
#            using 'polypeptide' rows from the GFF3 file.

################################################################################
# CHAPTER 1: BASIC FUNCTIONALITY
################################################################################

note('##############################################################################');
note('Starting Chapter 1: basic functionality');
note('##############################################################################');

# To test basic functionality with 'toy' GFF3 files,
# start out by allowing only the most basic GFF3 terms.
$lg_obj->param('gene_types',   ['gene']);
$lg_obj->param('mrna_types',   ['mRNA']);
$lg_obj->param('exon_types',   ['exon']);
$lg_obj->param('cds_types',    ['CDS']);
$lg_obj->param('utr_types',    ['five_prime_UTR', 'three_prime_UTR']);
$lg_obj->param('ignore_types', []);
$lg_obj->param('types_complete', 1);

# load_db method
my $db = $lg_obj->load_db($gff3_file, $fasta_file);
isa_ok($db, 'Bio::DB::SeqFeature::Store', 'load_db method: db');
is($db->get_features_by_type('gene'), 1, 'load_db method: one gene');
is($db->get_features_by_type('mRNA'), 1, 'load_db method: one transcript');
is($db->get_features_by_type('exon'), 5, 'load_db method: five exons');
is($db->get_features_by_type('CDS'), 1, 'load_db method: five CDS');
is($db->get_features_by_type('five_prime_UTR'), 1, 'load_db method: one 5\' UTR');
is($db->get_features_by_type('three_prime_UTR'), 1, 'load_db method: one 3\' UTR');

# check_db method
throws_ok { $lg_obj->check_db($db) } qr/Unrecognised types in GFF3 file/, 'check_db method: throw on unrecognised types';

$lg_obj->param('types_complete', 0);
lives_ok { $lg_obj->check_db($db) } 'check_db method: unrecognised types are fine if types_complete = 0';

$lg_obj->param('ignore_types', ['region', 'chromosome']);
$lg_obj->param('types_complete', 1);
lives_ok { $lg_obj->check_db($db) } 'check_db method: ignore types in explicit ignore_types list';

# fetch_slices method
my %slices = $lg_obj->fetch_slices($dba);
is(scalar(keys %slices), 3, 'fetch_slice method: hash with three keys returned');
ok(exists($slices{'2L'}), 'fetch_slice method: slice for chr 2L');
ok(exists($slices{'CM000356.1'}), 'fetch_slice method: slice for chr 2L synonym (CM000356.1)');
ok(exists($slices{'Mt'}), 'fetch_slice method: slice for chr Mt');
isa_ok($slices{'2L'}, 'Bio::EnsEMBL::Slice', 'fetch_slice method: hash value');

# check_seq_ids method
my @gff_genes = $db->get_features_by_type('gene');
my $gff_gene  = $gff_genes[0];

lives_ok { $lg_obj->check_seq_ids(\%slices, \@gff_genes) } 'check_seq_ids method: fine if sequence name matches';

$gff_gene->seq_id('CM000356.1');
lives_ok { $lg_obj->check_seq_ids(\%slices, \@gff_genes) } 'check_seq_ids method: fine if sequence synonym matches';

$gff_gene->seq_id('chr2L');
throws_ok { $lg_obj->check_seq_ids(\%slices, \@gff_genes) } qr/Unrecognised sequences in GFF3 file/, 'check_seq_ids method: throw on unrecognised sequences';

$gff_gene->seq_id('2L');

# new_gene method
my $slice = $slices{$gff_gene->seq_id};
my $gene  = $lg_obj->new_gene($gff_gene, $slice, $analysis);
isa_ok($gene, 'Bio::EnsEMBL::Gene', 'new_gene method: gene');
is($gene->stable_id, 'AGAP004700_test', 'new_gene method: stable_id correct');
is($gene->biotype, 'protein_coding', 'new_gene method: biotype correct');
is($gene->source, $gene_source, 'new_gene method: source correct');
is($gene->seq_region_name, '2L', 'new_gene method: seq_region_name correct');
is($gene->start, 2013634, 'new_gene method: start correct');
is($gene->end, 2015127, 'new_gene method: end correct');
is($gene->strand, 1, 'new_gene method: strand correct');
is($gene->analysis->logic_name, $logic_name, 'new_gene method: analysis correct');

# new_transcript method
my @gff_transcripts = $gff_gene->get_SeqFeatures('mRNA');
my $gff_transcript  = $gff_transcripts[0];

my $transcript = $lg_obj->new_transcript($gff_transcript, $gene);
isa_ok($transcript, 'Bio::EnsEMBL::Transcript', 'new_transcript method: transcript');
is($transcript->stable_id, 'AGAP004700_test-RA', 'new_transcript method: stable_id correct');
is($transcript->biotype, 'protein_coding', 'new_transcript method: biotype correct');
is($transcript->source, $gene_source, 'new_transcript method: source correct');
is($transcript->seq_region_name, '2L', 'new_transcript method: seq_region_name correct');
is($transcript->start, 2013634, 'new_transcript method: start correct');
is($transcript->end, 2015127, 'new_transcript method: end correct');
is($transcript->strand, 1, 'new_transcript method: strand correct');
is($transcript->analysis->logic_name, $logic_name, 'new_transcript method: analysis correct');

# new_exon method
my @gff_exons = sort sort_coding $gff_transcript->get_SeqFeatures('exon');
my $gff_exon  = $gff_exons[0];
my $exon_id   = 'AGAP004700_test-RA-E1';

my $exon = $lg_obj->new_exon($gff_exon, $transcript, $exon_id);
isa_ok($exon, 'Bio::EnsEMBL::Exon', 'new_exon method: exon');
is($exon->stable_id, $exon_id, 'new_exon method: stable_id correct');
is($exon->seq_region_name, '2L', 'new_exon method: seq_region_name correct');
is($exon->start, 2013634, 'new_exon method: start correct');
is($exon->end, 2013907, 'new_exon method: end correct');
is($exon->strand, 1, 'new_exon method: strand correct');
is($exon->phase, -1, 'new_exon method: phase initially set to -1');
is($exon->end_phase, -1, 'new_exon method: end phase initially set to -1');

# add_exons method
$lg_obj->add_exons(\@gff_exons, $transcript, $lg_obj->param('prediction'));
my @exons = @{ $transcript->get_all_Exons };
is(scalar(@exons), 5, 'add_exons method: added 5 exons to transcript');
is($exon->rank($transcript), 1, 'add_exons method: rank correct');

# get_cds method
my @gff_cds = $lg_obj->get_cds($gff_transcript);
is(scalar(@gff_cds), 5, 'get_cds method: split single CDS object into five segments');

# get_cds_id method
my ($cds_id, undef) = $lg_obj->get_cds_id($gff_transcript);
is($cds_id, 'AGAP004700_test-PA', 'get_cds_id method: id correct');

# get_utr method
my @gff_utr = $lg_obj->get_utr($gff_transcript);
my $five_prime_utr = $gff_utr[0];
my $three_prime_utr = $gff_utr[-1];
is(scalar(@gff_utr), 2, 'get_utr method: two UTRs');

# infer_translation method
my ($translation_id, undef, $genomic_start, $genomic_end) =
  $lg_obj->infer_translation($gff_transcript, $transcript);
is($translation_id, 'AGAP004700_test-PA', 'infer_translation method: stable_id correct');
is($genomic_start, 2013727, 'infer_translation method: start of translation correct');
is($genomic_end, 2014922, 'infer_translation method: end of translation correct');

# translation_coordinates method
my ($start_exon, $end_exon, $seq_start, $seq_end) =
  $lg_obj->translation_coordinates($transcript, $genomic_start, $genomic_end);
is($start_exon->start, 2013634, 'translation_coordinates method: start exon correct');
is($end_exon->end, 2015127, 'translation_coordinates method: end exon correct');
is($seq_start, 94, 'translation_coordinates method: translation start correct');
is($seq_end, 341, 'translation_coordinates method: translation end correct');
is($five_prime_utr->length, $seq_start - 1, 'translation_coordinates method: translation start correct');
is($three_prime_utr->length, $gff_exons[-1]->length - $seq_end, 'translation_coordinates method: translation end correct');

# new_translation method
my $translation = $lg_obj->new_translation($translation_id, $start_exon, $end_exon, $seq_start, $seq_end);
isa_ok($translation, 'Bio::EnsEMBL::Translation', 'new_translation method: translation');
is($translation->stable_id, 'AGAP004700_test-PA', 'new_translation method: stable_id correct');
is($translation->start_Exon->stable_id, 'AGAP004700_test-RA-E1', 'new_translation method: start exon correct');
is($translation->end_Exon->stable_id, 'AGAP004700_test-RA-E5', 'new_translation method: end exon correct');
is($translation->start, 94, 'new_translation method: seq_start correct');
is($translation->end, 341, 'new_translation method: seq_end correct');

$transcript->translation($translation);

# set_exon_phase method
$lg_obj->set_exon_phase($transcript);
my $exons = $transcript->get_all_Exons;
is($$exons[0]->phase, -1, 'set_exon_phase method: exon 1 phase correct');
is($$exons[0]->end_phase, 1, 'set_exon_phase method: exon 1 end phase correct');
is($$exons[1]->phase, 1, 'set_exon_phase method: exon 2 phase correct');
is($$exons[1]->end_phase, 0, 'set_exon_phase method: exon 2 end phase correct');
is($$exons[2]->phase, 0, 'set_exon_phase method: exon 3 phase correct');
is($$exons[2]->end_phase, 2, 'set_exon_phase method: exon 3 end phase correct');
is($$exons[3]->phase, 2, 'set_exon_phase method: exon 4 phase correct');
is($$exons[3]->end_phase, 1, 'set_exon_phase method: exon 4 end phase correct');
is($$exons[4]->phase, 1, 'set_exon_phase method: exon 5 phase correct');
is($$exons[4]->end_phase, -1, 'set_exon_phase method: exon 5 end phase correct');

# load_genes method
my $cdna_seq   = 'TTCGTCTGAGGTATATTATAGGAGATTCCTAAATCTAAATTCAGTTAGAGTGAATAGTTTAGAGTGAAGCGAATCAGTCGCACAAGTGTTAACATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAGAAAGTAGTGCATACAAACTGGTTAGGTTTGATGTGCAACATTCGTACTAAGTTGATAGTGATTTTGTAAATTTCACTAAGATTTCAACGGTTCCCTGTAATATAAAAGTATGAAATTATAGGCATTACTATATTATCAAATCATTGATTAACAAAATCGATAATTTAATATAAAATATAAACCGTTAAATAAAAGATGGTTGTTT';
my $coding_seq = 'ATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAG';
my $aa_seq     = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPL';

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->load_genes($db);

my $db_gene = $ga->fetch_by_stable_id('AGAP004700_test');
isa_ok($db_gene, 'Bio::EnsEMBL::Gene', 'load_genes method: gene');

my $db_transcripts = $db_gene->get_all_Transcripts();
is(scalar(@{$db_transcripts}), 1, 'load_genes method: one transcript stored');

my $db_transcript = $$db_transcripts[0];
is($db_transcript->stable_id, 'AGAP004700_test-RA', 'load_genes method: transcript stable_id correct');

my $db_exons = $db_transcript->get_all_Exons();
is(scalar(@{$db_exons}), 5, 'load_genes method: five exons stored');

my $db_cdna_seq = $db_transcript->seq->seq();
is($db_cdna_seq, $cdna_seq, 'load_genes method: correct cDNA sequence stored');

my $db_coding_seq = $db_transcript->translateable_seq();
is($db_coding_seq, $coding_seq, 'load_genes method: correct coding sequence stored');

my $db_aa_seq = $db_transcript->translate->seq();
is($db_aa_seq, $aa_seq, 'load_genes method: correct amino acid sequence stored');

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

################################################################################
# CHAPTER 2: NAME ATTRIBUTE
################################################################################

note('##############################################################################');
note('Starting Chapter 2: Name attribute');
note('##############################################################################');

# load_genes method, using names as stable IDs
{
$lg_obj->param('use_name_field', 'stable_id');

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->load_genes($db);

my $db_gene = $ga->fetch_by_stable_id('test_gene');
isa_ok($db_gene, 'Bio::EnsEMBL::Gene', 'load_genes method (Name as stable_id): gene');

my $db_transcripts = $db_gene->get_all_Transcripts();
is(scalar(@{$db_transcripts}), 1, 'load_genes method (Name as stable_id): one transcript stored');

my $db_transcript = $$db_transcripts[0];
is($db_transcript->stable_id, 'test_transcript', 'load_genes method (Name as stable_id): transcript stable_id correct');

my $db_translation = $db_transcript->translation;
is($db_translation->stable_id, 'test_translation', 'load_genes method (Name as stable_id): translation stable_id correct');

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));
}

# load_genes method, using names as xrefs
{
$lg_obj->param('use_name_field', 'xref');
$lg_obj->param('xref_gene_external_db', 'RefSeq_gene_name');
$lg_obj->param('xref_transcript_external_db', 'RefSeq_mRNA');
$lg_obj->param('xref_translation_external_db', 'RefSeq_peptide');

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->load_genes($db);

my $db_gene = $ga->fetch_by_stable_id('AGAP004700_test');
isa_ok($db_gene, 'Bio::EnsEMBL::Gene', 'load_genes method (Name as xref): gene');

my $db_transcripts = $db_gene->get_all_Transcripts();
is(scalar(@{$db_transcripts}), 1, 'load_genes method (Name as xref): one transcript stored');

my $db_transcript = $$db_transcripts[0];
is($db_transcript->stable_id, 'AGAP004700_test-RA', 'load_genes method (Name as xref): transcript stable_id correct');

my $db_translation = $db_transcript->translation;
is($db_translation->stable_id, 'AGAP004700_test-PA', 'load_genes method (Name as stable_id): translation stable_id correct');

my $gene_xrefs = $db_gene->get_all_DBEntries('RefSeq_gene_name');
is(scalar(@{$gene_xrefs}), 1, 'load_genes method (Name as xref): one gene xref stored');

my $gene_xref = $$gene_xrefs[0];
is($gene_xref->primary_id, 'test_gene', 'load_genes method (Name as xref): gene xref primary_id correct');
is($gene_xref->display_id, 'test_gene', 'load_genes method (Name as xref): gene xref display_id correct');
is($gene_xref->analysis->logic_name, $logic_name, 'load_genes method (Name as xref): gene xref analysis correct');

my $transcript_xrefs = $db_transcript->get_all_DBEntries('RefSeq_mRNA');
is(scalar(@{$transcript_xrefs}), 1, 'load_genes method (Name as xref): one gene xref stored');

my $transcript_xref = $$transcript_xrefs[0];
is($transcript_xref->primary_id, 'test_transcript', 'load_genes method (Name as xref): transcript xref primary_id correct');
is($transcript_xref->display_id, 'test_transcript', 'load_genes method (Name as xref): transcript xref display_id correct');
is($transcript_xref->analysis->logic_name, $logic_name, 'load_genes method (Name as xref): transcript xref analysis correct');

my $translation_xrefs = $db_translation->get_all_DBEntries('RefSeq_peptide');
is(scalar(@{$translation_xrefs}), 1, 'load_genes method (Name as xref): one gene xref stored');

my $translation_xref = $$translation_xrefs[0];
is($translation_xref->primary_id, 'test_translation', 'load_genes method (Name as xref): translation xref primary_id correct');
is($translation_xref->display_id, 'test_translation', 'load_genes method (Name as xref): translation xref display_id correct');
is($translation_xref->analysis->logic_name, $logic_name, 'load_genes method (Name as xref): translation xref analysis correct');

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));
}

# run method, using non-unique names as stable IDs
{
$lg_obj->param('gff3_file', $ncrna_gff3_file);
$lg_obj->param('gene_types', ['gene', 'miRNA_gene']);
$lg_obj->param('mrna_types', ['mRNA', 'miRNA', 'tRNA']);

$lg_obj->param('use_name_field', 'stable_id');

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->run();

my $gene_1 = $ga->fetch_by_stable_id('mir_test');
isa_ok($gene_1, 'Bio::EnsEMBL::Gene', 'run method: gene');

my $transcripts_1 = $gene_1->get_all_Transcripts();
is(scalar(@{$transcripts_1}), 1, 'run method: one transcript stored');

my $transcript_1 = $$transcripts_1[0];
is($transcript_1->stable_id, 'AGAP013620_test_1-RA', 'run method: transcript stable_id correct');

my $gene_2 = $ga->fetch_by_stable_id('mir_test_1');
isa_ok($gene_2, 'Bio::EnsEMBL::Gene', 'run method: gene');

my $transcripts_2 = $gene_2->get_all_Transcripts();
is(scalar(@{$transcripts_2}), 1, 'run method: one transcript stored');

my $transcript_2 = $$transcripts_2[0];
is($transcript_2->stable_id, 'AGAP013620_test_2-RA', 'run method: transcript stable_id correct');

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->param('gff3_file', $gff3_file);
$lg_obj->param('gene_types', ['gene']);
$lg_obj->param('mrna_types', ['mRNA']);
}

$lg_obj->param('use_name_field', 'no');

################################################################################
# CHAPTER 3: PREDICTED TRANSCRIPTS
################################################################################

note('##############################################################################');
note('Starting Chapter 3: predicted transcripts');
note('##############################################################################');

$lg_obj->param('prediction', 1);

# new_predicted_transcript method
my $ptranscript = $lg_obj->new_predicted_transcript($gff_transcript, $gene);
isa_ok($ptranscript, 'Bio::EnsEMBL::PredictionTranscript', 'new_predicted_transcript method: transcript');
is($ptranscript->display_label, 'AGAP004700_test-RA', 'new_predicted_transcript method: display_label correct');
is($ptranscript->seq_region_name, '2L', 'new_predicted_transcript method: seq_region_name correct');
is($ptranscript->start, 2013634, 'new_predicted_transcript method: start correct');
is($ptranscript->end, 2015127, 'new_predicted_transcript method: end correct');
is($ptranscript->strand, 1, 'new_predicted_transcript method: strand correct');
is($ptranscript->analysis->logic_name, $logic_name, 'new_predicted_transcript method: analysis correct');

# new_predicted_exon method
my $pexon = $lg_obj->new_predicted_exon($gff_exon, $ptranscript);
isa_ok($pexon, 'Bio::EnsEMBL::PredictionExon', 'new_predicted_exon method: exon');
is($pexon->seq_region_name, '2L', 'new_predicted_exon method: seq_region_name correct');
is($pexon->start, 2013634, 'new_predicted_exon method: start correct');
is($pexon->end, 2013907, 'new_predicted_exon method: end correct');
is($pexon->strand, 1, 'new_predicted_exon method: strand correct');
is($pexon->phase, -1, 'new_predicted_exon method: phase correct');
is($pexon->score, 200, 'new_predicted_exon method: score correct');
is($pexon->p_value, 0.01, 'new_predicted_exon method: p_value correct');

# add_exons method
$lg_obj->add_exons(\@gff_exons, $ptranscript, $lg_obj->param('prediction'));
my @pexons = @{ $ptranscript->get_all_Exons };
is(scalar(@pexons), 5, 'add_exons method: added 5 exons to transcript');
is($pexon->rank($ptranscript), 1, 'add_exons method: rank correct');

# load_genes method
$testdb->hide($dbtype, qw(prediction_exon prediction_transcript));

$lg_obj->load_genes($db);

my $db_ptranscripts = $pta->fetch_all_by_logic_name($logic_name);
is(scalar(@{$db_ptranscripts}), 1, 'load_genes method: one prediction transcript stored');

my $db_ptranscript = $$db_ptranscripts[0];
isa_ok($db_ptranscript, 'Bio::EnsEMBL::PredictionTranscript', 'load_genes method: transcript');
is($db_ptranscript->display_label, 'AGAP004700_test-RA', 'load_genes method: transcript stable_id correct');

my $db_pexons = $db_ptranscript->get_all_Exons();
is(scalar(@{$db_pexons}), 5, 'load_genes method: five exons stored');

my $db_cdna_seq_p = $db_ptranscript->seq->seq();
is($db_cdna_seq_p, $cdna_seq, 'load_genes method: correct cDNA sequence stored');

$testdb->restore($dbtype, qw(prediction_exon prediction_transcript));

$lg_obj->param('prediction', 0);

################################################################################
# CHAPTER 4: NCRNA GENES
################################################################################

note('##############################################################################');
note('Starting Chapter 4: ncRNA genes');
note('##############################################################################');

$lg_obj->param('gff3_file', $ncrna_gff3_file);
$lg_obj->param('gene_types', ['gene', 'miRNA_gene']);
$lg_obj->param('mrna_types', ['mRNA', 'miRNA', 'tRNA']);

# run method
my $mirna_seq = 'ATAGGATTTGGCTTCACTGCTAAGATTCATTGAACTATTTTTTAAAAATAATTTTCTACATGATAATTTGTA';
my $trna_seq  = 'TATATTTTAGTGTATGATGCACAAAAGATTTTGATTCTTTTAGAAACAGTTTAATTCTGTTAAGTATAA';

$testdb->hide($dbtype, qw(exon exon_transcript gene meta_coord transcript));

$lg_obj->run();

my $rna_genes = $ga->fetch_all_by_logic_name($logic_name);
is(scalar(@{$rna_genes}), 3, 'run method: three genes stored');

my $rna_gene_1 = $ga->fetch_by_stable_id('AGAP013620_test_1');
isa_ok($rna_gene_1, 'Bio::EnsEMBL::Gene', 'run method: gene');
is($rna_gene_1->biotype, 'miRNA', 'run method: gene biotype correct');

my $rna_transcripts_1 = $rna_gene_1->get_all_Transcripts();
is(scalar(@{$rna_transcripts_1}), 1, 'run method: one transcript stored');

my $rna_transcript_1 = $$rna_transcripts_1[0];
is($rna_transcript_1->stable_id, 'AGAP013620_test_1-RA', 'run method: transcript stable_id correct');
is($rna_transcript_1->biotype, 'miRNA', 'run method: transcript biotype correct');

my $rna_exons_1 = $rna_transcript_1->get_all_Exons();
is(scalar(@{$rna_exons_1}), 1, 'run method: one exon stored');

my $rna_transcript_1_seq = $rna_transcript_1->seq->seq();
is($rna_transcript_1_seq, $mirna_seq, 'run method: correct DNA sequence stored');

my $rna_gene_2 = $ga->fetch_by_stable_id('AGAP013620_test_2');
isa_ok($rna_gene_2, 'Bio::EnsEMBL::Gene', 'run method: gene');
is($rna_gene_2->biotype, 'miRNA', 'run method: gene biotype correct');

my $rna_transcripts_2 = $rna_gene_2->get_all_Transcripts();
is(scalar(@{$rna_transcripts_2}), 1, 'run method: one transcript stored');

my $rna_transcript_2 = $$rna_transcripts_2[0];
is($rna_transcript_2->stable_id, 'AGAP013620_test_2-RA', 'run method: transcript stable_id correct');
is($rna_transcript_2->biotype, 'miRNA', 'run method: transcript biotype correct');

my $rna_exons_2 = $rna_transcript_2->get_all_Exons();
is(scalar(@{$rna_exons_2}), 1, 'run method: one exon stored');

my $rna_transcript_2_seq = $rna_transcript_2->seq->seq();
is($rna_transcript_2_seq, $mirna_seq, 'run method: correct DNA sequence stored');

my $rna_gene_3 = $ga->fetch_by_stable_id('AGAP028358_test');
isa_ok($rna_gene_3, 'Bio::EnsEMBL::Gene', 'run method: gene');
is($rna_gene_3->biotype, 'tRNA', 'run method: gene biotype correct');

my $rna_transcripts_3 = $rna_gene_3->get_all_Transcripts();
is(scalar(@{$rna_transcripts_3}), 1, 'run method: one transcript stored');

my $rna_transcript_3 = $$rna_transcripts_3[0];
is($rna_transcript_3->stable_id, 'AGAP028358_test-RA', 'run method: transcript stable_id correct');
is($rna_transcript_3->biotype, 'tRNA', 'run method: transcript biotype correct');

my $rna_exons_3 = $rna_transcript_3->get_all_Exons();
is(scalar(@{$rna_exons_3}), 1, 'run method: one exon stored');

my $rna_transcript_3_seq = $rna_transcript_3->seq->seq();
is($rna_transcript_3_seq, $trna_seq, 'run method: correct DNA sequence stored');

$testdb->restore($dbtype, qw(exon exon_transcript gene meta_coord transcript));

$lg_obj->param('gene_types', ['gene']);
$lg_obj->param('mrna_types', ['mRNA']);
$lg_obj->param('gff3_file', $gff3_file);

################################################################################
# CHAPTER 5: PSEUDOGENES
################################################################################

note('##############################################################################');
note('Starting Chapter 5: pseudogenes');
note('##############################################################################');

$lg_obj->param('gff3_file', $pseudo_gff3_file);
$lg_obj->param('gene_types', ['gene', 'pseudogene']);
$lg_obj->param('mrna_types', ['transcript', 'pseudogenic_transcript']);
$lg_obj->param('exon_types', ['exon', 'pseudogenic_exon']);

my $pseudo_db = $lg_obj->load_db($pseudo_gff3_file, $fasta_file);
$lg_obj->check_db($pseudo_db);

my @pseudo_gff_genes = (
  $pseudo_db->get_features_by_type('gene'),
  $pseudo_db->get_features_by_type('pseudogene'),
);

foreach my $pseudo_gff_gene (@pseudo_gff_genes) {
  my $slice = $slices{$pseudo_gff_gene->seq_id};
  my $pseudo_gene = $lg_obj->new_gene($pseudo_gff_gene, $slice, $analysis);
  
  # add_pseudogenic_transcript method
  if ($pseudo_gene->stable_id eq 'pseudo_test_4') {
    my $pseudo_transcript = $lg_obj->add_pseudogenic_transcript($pseudo_db, $pseudo_gff_gene, $pseudo_gene);
    isa_ok($pseudo_transcript, 'Bio::EnsEMBL::Transcript', 'add_pseudogenic_transcript method: transcript');
    is($pseudo_gene->biotype, 'pseudogene', 'add_pseudogenic_transcript method: gene biotype correct');
    is($pseudo_transcript->biotype, 'pseudogene', 'add_pseudogenic_transcript method: transcript biotype correct');
    is($pseudo_transcript->seq_region_name, '2L', 'add_pseudogenic_transcript method: transcript seq_region_name correct');
    is($pseudo_transcript->start, 2439935, 'add_pseudogenic_transcript method: transcript start correct');
    is($pseudo_transcript->end, 2440006, 'add_pseudogenic_transcript method: transcript end correct');
    is($pseudo_transcript->strand, 1, 'add_pseudogenic_transcript method: transcript strand correct');
  }
  
  # infer_exons method
  if ($pseudo_gene->stable_id eq 'pseudo_test_8') {
    my $pseudo_gff_transcript = $lg_obj->infer_transcript($pseudo_db, $pseudo_gff_gene);
    my $pseudo_gff_exons = $lg_obj->infer_exons($pseudo_db, $pseudo_gff_transcript);
    
    my @pseudo_gff_exons = $pseudo_gff_transcript->get_SeqFeatures('exon');
    is(scalar(@pseudo_gff_exons), 1, 'infer_exons method: one gff_exon added to gff_transcript');
    
    my $pseudo_gff_exon = $$pseudo_gff_exons[0];
    is($pseudo_gff_exon->seq_id, '2L', 'infer_exons method: exon seq_id correct');
    is($pseudo_gff_exon->start, 2439935, 'infer_exons method: exon start correct');
    is($pseudo_gff_exon->end, 2440006, 'infer_exons method: exon end correct');
    is($pseudo_gff_exon->strand, 1, 'infer_exons method: exon strand correct');
  }
  
  # set_nontranslating_gene method
  if ($pseudo_gene->stable_id eq 'pseudo_test_10') {
    my @pseudo_gff_transcripts = $pseudo_gff_gene->get_SeqFeatures('transcript');
    my $pseudo_gff_transcript = $pseudo_gff_transcripts[0];
    
    my $pseudo_transcript = $lg_obj->add_transcript($pseudo_db, $pseudo_gff_transcript, $pseudo_gene);
    $pseudo_gene->add_Transcript($pseudo_transcript);
    
    $lg_obj->add_translation($pseudo_db, $pseudo_gff_transcript, $pseudo_gene, $pseudo_transcript);
    
    $lg_obj->set_nontranslating_gene($pseudo_gene);
    is($pseudo_gene->biotype, 'nontranslating_CDS', 'set_nontranslating_gene method: gene biotype correct');
    is($pseudo_transcript->biotype, 'nontranslating_CDS', 'set_nontranslating_gene method: transcript biotype correct');
    
    $lg_obj->param('nontranslating', 'pseudogene');
    $lg_obj->set_nontranslating_gene($pseudo_gene);
    is($pseudo_gene->biotype, 'pseudogene', 'set_nontranslating_gene method: gene biotype correct');
    is($pseudo_transcript->biotype, 'pseudogene', 'set_nontranslating_gene method: transcript biotype correct');
  }
}

# load_genes method
my $pseudo_seq = 'ATAGGATTTGGCTTCACTGCTAAGATTCATTGAACTATTTTTTAAAAATAATTTTCTACATGATAATTTGTA';

$testdb->hide($dbtype, qw(exon exon_transcript gene meta_coord transcript));

$lg_obj->load_genes($pseudo_db);

my $pseudo_genes = $ga->fetch_all_by_logic_name($logic_name);
is(scalar(@{$pseudo_genes}), 10, 'load_genes method: ten pseudogenes stored');

foreach my $pseudo_gene (@{$pseudo_genes}) {
  is($pseudo_gene->biotype, 'pseudogene', 'load_genes method: gene biotype correct for '.$pseudo_gene->stable_id);
  my $pseudo_transcripts = $pseudo_gene->get_all_Transcripts();
  is(scalar(@{$pseudo_transcripts}), 1, 'load_genes method: one transcript stored');
  
  my $pseudo_transcript = $$pseudo_transcripts[0];
  is($pseudo_transcript->biotype, 'pseudogene', 'load_genes method: transcript biotype correct for '.$pseudo_transcript->stable_id);
}

$testdb->restore($dbtype, qw(exon exon_transcript gene meta_coord transcript));

$lg_obj->param('gene_types', ['gene']);
$lg_obj->param('mrna_types', ['mRNA']);
$lg_obj->param('exon_types', ['exon']);
$lg_obj->param('gff3_file', $gff3_file);
$lg_obj->param('nontranslating', 'nontranslating_CDS');

################################################################################
# CHAPTER 6: TRANSLATION SCENARIOS
################################################################################

note('##############################################################################');
note('Starting Chapter 6: translation scenarios');
note('##############################################################################');

$lg_obj->param('gff3_file', $transl_gff3_file);
$lg_obj->param('ignore_types', ['region', 'chromosome', 'polypeptide', 'protein']);
$lg_obj->param('polypeptides', 1);

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->run();

my $test_genes = $ga->fetch_all_by_logic_name($logic_name);
is(scalar(@{$test_genes}), 11, 'run method: eleven genes stored');

foreach my $test_gene (@{$test_genes}) {
  note('gene:'.$test_gene->stable_id);
  if ($test_gene->stable_id eq 'test_negative_strand') {
    my $cdna   = 'CTGGCATTGAAATGCTTCATTTCAATATAAAAAGAACATTTGTTTAGTGTATATTATTCAAACTCGTTTCACCAAATAGAGTTTTTGTATATATATACAGCAAATTTCTTTTGAGCGGTACTCAGTACAAAACTAACTTTGTTTCCTTGTACCCGCAGATCGTTACGTCCGACGGTTACCTATTCAACTAAAACAAACGGATCGAACGAACGGACTTAGGAAACATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAATCAAAGTCATTCATAGATAAAAAAAATTCTACTTTTAACTTTTTTCC';
    my $coding = 'ATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAA';
    my $aa     = 'MRLTFTLILLGICIPLQAKPQLLGGGGGLLGTGVGGSTGLLGTGVLSSGTTSGLLGTGLGATTGVLGSGILGTGTTSGLLGTGLGATTGVLGTGLLGTGTTSGLLGTGLGATSGVLGTGLLGTGTTSGLLGTGILGTGTTSGLLGTGIGATSGLLGTGILSGGATTTTAAVTTTTAAATTTAVPTTTTAVPTTTTAVPTTTTAAATTTVTATTTTTTTTIAPATTTNAPTTTTIAPATTTNAPTTTVTSAPTGTVTFNINGNIITINANDQATIDAINAVLAASSTTTAMSTTTTAAPSFGSTVSFNIGGSIITVDAADTATIAAIIQSLNSGSTSSTAISTTSVFPVTATTTTTVVTPVTTTSTVSALATTTTLGSSATSTTSSTSAPNSGSNSISIIININGQFILISASNQSLLLLMQQILQQISSITTITSTTQPPAIPIQPTSAPVQTSGGRPGGCMRRHQSSEEVGGNRGLRGGFGFGIGIRGPLGGGFQAGFGGGISNRMPGAGIGAGIRLQGPLGGSIQTGFGAGLGRAHASVSGRTQGAIGSGFRVQGPLGGVVQAGFGAGFRAGGHRHG';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
  }
  
  if ($test_gene->stable_id eq 'test_nontranslating') {
    my $cdna   = 'CTGGCATTGAAATGCTTCATTTCAATATAAAAAGAACATTTGTTTAGTGTATATTATTCAAACTCGTTTCACCAAATAGAGTTTTTGTATATATATACAGCAAATTTCTTTTGAGCGGTACTCAGTACAAAACTAACTTTGTTTCCTTGTACCCGCAGATCGTTACGTCCGACGGTTACCTATTCAACTAAAACAAACGGATCGAACGAACGGACTTAGGAAACATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAATCAAAGTCATTCATAGATAAAAAAAATTCTACTTTTAACTTTTTTCC';
    my $coding = 'CTGGCATTGAAATGCTTCATTTCAATATAAAAAGAACATTTGTTTAGTGTATATTATTCAAACTCGTTTCACCAAATAGAGTTTTTGTATATATATACAGCAAATTTCTTTTGAGCGGTACTCAGTACAAAACTAACTTTGTTTCCTTGTACCCGCAGATCGTTACGTCCGACGGTTACCTATTCAACTAAAACAAACGGATCGAACGAACGGACTTAGGAAACATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAATCAAAGTCATTCATAGATAAAAAAAATTCTACTTTTAACTTTTTTCC';
    my $aa     = 'LALKCFISI*KEHLFSVYYSNSFHQIEFLYIYTANFF*AVLSTKLTLFPCTRRSLRPTVTYSTKTNGSNERT*ETCD*LLH*SCWEFVFHCKPNHNCWEVVEACSELALEVLPVC*ELVYYRVVPQVDY*ELASVLLQEFLVLVY*EREQLVDY*VLV*EQQPGF*ALDYLEQELPVVCLEQVLELLVEFLELDYSEQELPVDCLELEFLERGLPVDC*EQV*ELLVVFLVLVFLVEVLLQLLLPLPQPLLPLPQQLFPPPQRLFPPPQRLFPLPQPLQPRPQ*QQQQRQQRRRLLQPQLLTHQQRQRLLQPQLLTHQQRQSHRPPQEL*HSTLMETL*L*MQMIKQL*MLSTLS*LLLRLLLL*VQQPQQHHPLEAQYHSTLVDLSLLLMLLIQLRLQRSFKVLTVDRLHRRQFQPLPYFL*LQQLRQLLSLQLRLRVLCQLWRQLQLWVPQLLQLHLAHQHPIADRTQFLLLST*TDNLF*YQPPIKAYCF*CNKFCSKFQA*QQLLQQHNLRPFQFNRRQLLYKLVEDVLVDA*DAISLLRKLEAIEDCVVDLDLALVFVDPLVEDFKLVLEVV*ATECLELVSVLAFDCRVL*EEAYRPGLEPV*DVPMPVYLVELKVQLDQGSVCKVL*EAWFKQDLGLDFELVDIDTDNQSHS*IKKILLLTFF';
    
    is($test_gene->biotype, 'nontranslating_CDS', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    is($transcript->biotype, 'nontranslating_CDS', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct (with stop codons) amino acid sequence stored');
  }
  
  if ($test_gene->stable_id eq 'test_no_utr') {
    my $cdna   = 'CTGGCATTGAAATGCTTCATTTCAATATAAAAAGAACATTTGTTTAGTGTATATTATTCAAACTCGTTTCACCAAATAGAGTTTTTGTATATATATACAGCAAATTTCTTTTGAGCGGTACTCAGTACAAAACTAACTTTGTTTCCTTGTACCCGCAGATCGTTACGTCCGACGGTTACCTATTCAACTAAAACAAACGGATCGAACGAACGGACTTAGGAAACATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAATCAAAGTCATTCATAGATAAAAAAAATTCTACTTTTAACTTTTTTCC';
    my $coding = 'ATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAA';
    my $aa     = 'MRLTFTLILLGICIPLQAKPQLLGGGGGLLGTGVGGSTGLLGTGVLSSGTTSGLLGTGLGATTGVLGSGILGTGTTSGLLGTGLGATTGVLGTGLLGTGTTSGLLGTGLGATSGVLGTGLLGTGTTSGLLGTGILGTGTTSGLLGTGIGATSGLLGTGILSGGATTTTAAVTTTTAAATTTAVPTTTTAVPTTTTAVPTTTTAAATTTVTATTTTTTTTIAPATTTNAPTTTTIAPATTTNAPTTTVTSAPTGTVTFNINGNIITINANDQATIDAINAVLAASSTTTAMSTTTTAAPSFGSTVSFNIGGSIITVDAADTATIAAIIQSLNSGSTSSTAISTTSVFPVTATTTTTVVTPVTTTSTVSALATTTTLGSSATSTTSSTSAPNSGSNSISIIININGQFILISASNQSLLLLMQQILQQISSITTITSTTQPPAIPIQPTSAPVQTSGGRPGGCMRRHQSSEEVGGNRGLRGGFGFGIGIRGPLGGGFQAGFGGGISNRMPGAGIGAGIRLQGPLGGSIQTGFGAGLGRAHASVSGRTQGAIGSGFRVQGPLGGVVQAGFGAGFRAGGHRHG';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');    
  }
  
  if ($test_gene->stable_id eq 'test_no_exon') {
    my $cdna   = 'CTGGCATTGAAATGCTTCATTTCAATATAAAAAGAACTTTGTTTCCTTGTACCCGCAGATCGTTACGTCCGACGGTTACCTATTCAACTAAAACAAACGGATCGAACGAACGGACTTAGGAAACATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAATCAAAGTCATTCATAGATAAAAAAAATTCTACTTTTAACTTTTTTCC';
    my $coding = 'ATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAA';
    my $aa     = 'MRLTFTLILLGICIPLQAKPQLLGGGGGLLGTGVGGSTGLLGTGVLSSGTTSGLLGTGLGATTGVLGSGILGTGTTSGLLGTGLGATTGVLGTGLLGTGTTSGLLGTGLGATSGVLGTGLLGTGTTSGLLGTGILGTGTTSGLLGTGIGATSGLLGTGILSGGATTTTAAVTTTTAAATTTAVPTTTTAVPTTTTAVPTTTTAAATTTVTATTTTTTTTIAPATTTNAPTTTTIAPATTTNAPTTTVTSAPTGTVTFNINGNIITINANDQATIDAINAVLAASSTTTAMSTTTTAAPSFGSTVSFNIGGSIITVDAADTATIAAIIQSLNSGSTSSTAISTTSVFPVTATTTTTVVTPVTTTSTVSALATTTTLGSSATSTTSSTSAPNSGSNSISIIININGQFILISASNQSLLLLMQQILQQISSITTITSTTQPPAIPIQPTSAPVQTSGGRPGGCMRRHQSSEEVGGNRGLRGGFGFGIGIRGPLGGGFQAGFGGGISNRMPGAGIGAGIRLQGPLGGSIQTGFGAGLGRAHASVSGRTQGAIGSGFRVQGPLGGVVQAGFGAGFRAGGHRHG';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 4, 'run method: four exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
  }
  
  if ($test_gene->stable_id eq 'test_nonzero_start_cds_phase') {
    my $cdna   = 'GCGTAACTCAGAAAATCAAGTAAACGGGCTTTCTAGCTCAAAGGAAAATACACAAGACAACGTCACGTGTTCGGTTCCGAGCCGCTTGAGCTCGCCTAGAATAGAAACGGTAAAATCGCTTGATCGCAGGTTTTGGAAATTATTGAGCAAACGACGGCGTGGAAGTTTCCAAGGAGAACTTACAACTGGTGCATCGTAGAATACTTTAATAACTACACGATTTCATTTAGAAACATAAATCACAGATGTACTCTGCACTGTTCTAGATTCAAATAGATACTGCAATATGAGGAACGCTCGGAGAGAATTAAGCTCAAACCTTAAAAAGATTAAAAATTTTAATTAACTCAGTCTTCAAACATTATCAATATGATTTAAATAAAGTCCTTTTGTCTAAAATTAGATAAAACATTATTTTAACGTAGGAAGTGCACTGTCTTTTTGCGATATTTTAGTAGAGGTGATGCCCGGGTCCTATTGCTGTATTCTTTAGAAATAATATATTAATGGTGAAAAATATGCTGGCGATCGTATTTCGATTTGATCGATGTATAAAAGA';
    my $coding = 'CGTAACTCAGAAAATCAAGTAAACGGGCTTTCTAGCTCAAAGGAAAATACACAAGACAACGTCACGTGTTCGGTTCCGAGCCGCTTGAGCTCGCCTAGAATAGAAACGGTAAAATCGCTTGATCGCAGGTTTTGGAAATTATTGAGCAAACGACGGCGTGGAAGTTTCCAAGGAGAACTTACAACTGGTGCATCGTAG';
    my $aa     = 'RNSENQVNGLSSSKENTQDNVTCSVPSRLSSPRIETVKSLDRRFWKLLSKRRRGSFQGELTTGAS';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 2, 'run method: two exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
  }
  
  if ($test_gene->stable_id eq 'test_overlapping_exons') {
    my $cdna   = 'GCGTAACTCAGAAAATCAAGTAAACGGGCTTTCTAGCTCAAAGGAAAATACACAAGACAACGTCACGTGTTCGGTTCCGAGCCGCTTGAGCTCGCCTAGAATAGAAACGGTAAAATCGCTTGATCGCAGGTTTTGGAAATTATTGAGCAAACGACGGCGTGGAAGTTTCCAAGGAGAACTTACAACTGGTGCATCGTAGAATACTTTAATAACTACACGATTTCATTTAGAAACATAAATCACAGATGTACTCTGCACTGTTCTAGATTCAAATAGATACTGCAATATGAGGAACGCTCGGAGAGAATTAAGCTCAAACCTTAAAAAGATTAAAAATTTTAATTAACTCAGTCTTCAAACATTATCAATATGATTTAAATAAAGTCCTTTTGTCTAAAATTAGATAAAACATTATTTTAACGTAGGAAGTGCACTGTCTTTTTGCGATATTTTAGTAGAGGTGATGCCCGGGTCCTATTGCTGTATTCTTTAGAAATAATATATTAATGGTGAAAAATATGCTGGCGATCGTATTTCGATTTGATCGATGTATAAAAGA';
    my $coding = 'CGTAACTCAGAAAATCAAGTAAACGGGCTTTCTAGCTCAAAGGAAAATACACAAGACAACGTCACGTGTTCGGTTCCGAGCCGCTTGAGCTCGCCTAGAATAGAAACGGTAAAATCGCTTGATCGCAGGTTTTGGAAATTATTGAGCAAACGACGGCGTGGAAGTTTCCAAGGAGAACTTACAACTGGTGCATCGTAG';
    my $aa     = 'RNSENQVNGLSSSKENTQDNVTCSVPSRLSSPRIETVKSLDRRFWKLLSKRRRGSFQGELTTGAS';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
  }
  
  if ($test_gene->stable_id eq 'test_negative_strand_overlapping_exons') {
    my $cdna   = 'CTGGCATTGAAATGCTTCATTTCAATATAAAAAGAACATTTGTTTAGTGTATATTATTCAAACTCGTTTCACCAAATAGAGTTTTTGTATATATATACAGCAAATTTCTTTTGAGCGGTACTCAGTACAAAACTAACTTTGTTTCCTTGTACCCGCAGATCGTTACGTCCGACGGTTACCTATTCAACTAAAACAAACGGATCGAACGAACGGACTTAGGAAACATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAATCAAAGTCATTCATAGATAAAAAAAATTCTACTTTTAACTTTTTTCC';
    my $coding = 'ATGCGACTAACTTTTACATTGATCTTGCTGGGAATTTGTATTCCATTGCAAGCCAAACCACAACTGTTGGGAGGTGGTGGAGGCTTGCTCGGAACTGGCGTTGGAGGTTCTACCGGTTTGCTAGGAACTGGTGTACTATCGAGTGGTACCACAAGTGGATTATTAGGAACTGGCCTCGGTGCTACTACAGGAGTTCTTGGTTCTGGTATACTAGGAACGGGAACAACTAGTGGATTATTAGGTACTGGTTTAGGAGCAACAACCGGGGTTCTAGGCACTGGACTACTTGGAACAGGAACTACCAGTGGTCTGCTTGGAACAGGTCTTGGAGCTACTAGTGGAGTTCTTGGAACTGGACTACTCGGAACAGGAACTACCAGTGGACTGCTTGGAACTGGAATTCTTGGAACGGGGACTACCAGTGGACTGTTAGGAACAGGTATAGGAGCTACTAGTGGTCTTCTTGGTACTGGTATTCTTAGTGGAGGTGCTACTACAACTACTGCTGCCGTTACCACAACCACTGCTGCCGCTACCACAACAGCTGTTCCCACCACCACAACGGCTGTTCCCACCACCACAACGGCTGTTCCCACTACCACAACCGCTGCAGCCACGACCACAGTAACAGCAACAACAACGACAACAACGACGACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAACGATTGCTCCAGCCACAACTACTAACGCACCAACAACGACAGTCACATCGGCCCCCACAGGAACTGTAACATTCAACATTAATGGAAACATTATAACTATAAATGCAAATGATCAAGCAACTATAGATGCTATCAACGCTGTCCTAGCTGCTTCTTCGACTACTACTGCTATGAGTACAACAACCACAGCAGCACCATCCTTTGGAAGCACAGTATCATTCAACATTGGTGGATCTATCATTACTGTTGATGCTGCTGATACAGCTACGATTGCAGCGATCATTCAAAGTCTTAACAGTGGATCGACTTCATCGACGGCAATTTCAACCACTTCCGTATTTCCTGTAACTGCAACAACTACGACAACTGTTGTCACTCCAGTTACGACTACGAGTACTGTGTCAGCTTTGGCGACAACTACAACTTTGGGTTCCTCAGCTACTTCAACTACATCTAGCACATCAGCACCCAATAGCGGATCGAACTCAATTTCTATTATTATCAACATAAACGGACAATTTATTTTAATATCAGCCTCCAATCAAAGCTTACTGCTTTTGATGCAACAAATTTTGCAGCAAATTTCAAGCATAACAACAATTACTTCAACAACACAACCTCCGGCCATTCCAATTCAACCGACGTCAGCTCCTGTACAAACTAGTGGAGGACGTCCTGGTGGATGCATGAGACGCCATCAGTCTTCTGAGGAAGTTGGAGGCAATCGAGGACTGCGTGGTGGATTTGGATTTGGCATTGGTATTCGTGGACCCCTTGGTGGAGGATTTCAAGCTGGTTTTGGAGGTGGTATAAGCAACCGAATGCCTGGAGCTGGTATCGGTGCTGGCATTCGACTGCAGGGTCCTTTAGGAGGAAGCATACAGACCGGGTTTGGAGCCGGTTTAGGACGTGCCCATGCCAGTGTATCTGGTAGAACTCAAGGTGCAATTGGATCAGGGTTCCGTGTGCAAGGTCCTCTAGGAGGCGTGGTTCAAGCAGGATTTGGGGCTGGATTTCGAGCTGGTGGACATAGACACGGATAA';
    my $aa     = 'MRLTFTLILLGICIPLQAKPQLLGGGGGLLGTGVGGSTGLLGTGVLSSGTTSGLLGTGLGATTGVLGSGILGTGTTSGLLGTGLGATTGVLGTGLLGTGTTSGLLGTGLGATSGVLGTGLLGTGTTSGLLGTGILGTGTTSGLLGTGIGATSGLLGTGILSGGATTTTAAVTTTTAAATTTAVPTTTTAVPTTTTAVPTTTTAAATTTVTATTTTTTTTIAPATTTNAPTTTTIAPATTTNAPTTTVTSAPTGTVTFNINGNIITINANDQATIDAINAVLAASSTTTAMSTTTTAAPSFGSTVSFNIGGSIITVDAADTATIAAIIQSLNSGSTSSTAISTTSVFPVTATTTTTVVTPVTTTSTVSALATTTTLGSSATSTTSSTSAPNSGSNSISIIININGQFILISASNQSLLLLMQQILQQISSITTITSTTQPPAIPIQPTSAPVQTSGGRPGGCMRRHQSSEEVGGNRGLRGGFGFGIGIRGPLGGGFQAGFGGGISNRMPGAGIGAGIRLQGPLGGSIQTGFGAGLGRAHASVSGRTQGAIGSGFRVQGPLGGVVQAGFGAGFRAGGHRHG';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 5, 'run method: five exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
  }
  
  if ($test_gene->stable_id eq 'test_multiple_transcripts') {
    my $cdna_1   = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGA';
    my $coding_1 = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGA';
    my $aa_1     = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALKPPSALRSGTASRLLANNGISMISSTSQRIGTALGNAGNRMADRPITQHGISGLSTSYGRLGTAVSSNRQIKDKRYWQALLQSKIQEINQETTKILKEKKFLDRERSARKLYEKRVKEAAKELTNLQSTLTSMNLALDNCTSGMTRQHLLNETVALRERNEHIQEQLEVIFKQRQQKDIENKALERNTEQEKNKVIEMINSLPEEDQHKYREYKALSENLRKQNAIYHSQISEMEKQKDRLNTMIMNSQSRSEAHRLKSKLKELLNKRNALREEENNRLSPAQEREKLINDVRSNNQALASIGKQLKIVEDQLIEKKETLQQIDQDLEEGNSERHVKYKELKKRDDVMSAFMDSFKSNMNQEQQCKCIKLIYSKCLIITQSLLCKHVTGIDTLKNQITYAIEQITMQGINMNGMYDAKLEGNGFTSKNDLNSHSGLMKEYKKLGIQLKQLQILEKRTVQQMNSLRQEETEALQSIQKYANLEIPRSEAIAKMNELTTILQEMEDKKRVTENVVDEARNRNHEIKINLKSNDTYRQISHLEDKLIDLMKDNKVLQETVKNIQQVYNAYQFK';
    my $cdna_2   = 'CCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGA';
    my $coding_2 = 'CCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGA';
    my $aa_2     = 'PPSALRSGTASRLLANNGISMISSTSQRIGTALGNAGNRMADRPITQHGISGLSTSYGRLGTAVSSNRQIKDKRYWQALLQSKIQEINQETTKILKEKKFLDRERSARKLYEKRVKEAAKELTNLQSTLTSMNLALDNCTSGMTRQHLLNETVALRERNEHIQEQLEVIFKQRQQKDIENKALERNTEQEKNKVIEMINSLPEEDQHKYREYKALSENLRKQNAIYHSQISEMEKQKDRLNTMIMNSQSRSEAHRLKSKLKELLNKRNALREEENNRLSPAQEREKLINDVRSNNQALASIGKQLKIVEDQLIEKKETLQQIDQDLEEGNSERHVKYKELKKRDDVMSAFMDSFKSNMNQEQQCKCIKLIYSKCLIITQSLLCKHVTGIDTLKNQITYAIEQITMQGINMNGMYDAKLEGNGFTSKNDLNSHSGLMKEYKKLGIQLKQLQILEKRTVQQMNSLRQEETEALQSIQKYANLEIPRSEAIAKMNELTTILQEMEDKKRVTENVVDEARNRNHEIKINLKSNDTYRQISHLEDKLIDLMKDNKVLQETVKNIQQVYNAYQFK';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 2, 'run method: two transcripts stored');
    
    foreach my $transcript (@$transcripts) {
      if ($transcript->stable_id eq 'test_multiple_transcripts-RA') {
        is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
        is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
        is($transcript->seq->seq(), $cdna_1, 'run method: correct cDNA sequence stored');
        is($transcript->translateable_seq(), $coding_1, 'run method: correct coding sequence stored');
        is($transcript->translate->seq(), $aa_1, 'run method: correct amino acid sequence stored');
      }
      if ($transcript->stable_id eq 'test_multiple_transcripts-RB') {
        is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
        is(scalar(@{$transcript->get_all_Exons()}), 2, 'run method: two exons stored');
        is($transcript->seq->seq(), $cdna_2, 'run method: correct cDNA sequence stored');
        is($transcript->translateable_seq(), $coding_2, 'run method: correct coding sequence stored');
        is($transcript->translate->seq(), $aa_2, 'run method: correct amino acid sequence stored');
      }
    }
  }
  
  if ($test_gene->stable_id eq 'test_exon_phase') {
    my $cdna_1   = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGA';
    my $coding_1 = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGA';
    my $aa_1     = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALKPPSALRSGTASRLLANNGISMISSTSQRIGTALGNAGNRMADRPITQHGISGLSTSYGRLGTAVSSNRQIKDKRYWQALLQSKIQEINQETTKILKEKKFLDRERSARKLYEKRVKEAAKELTNLQSTLTSMNLALDNCTSGMTRQHLLNETVALRERNEHIQEQLEVIFKQRQQKDIENKALERNTEQEKNKVIEMINSLPEEDQHKYREYKALSENLRKQNAIYHSQISEMEKQKDRLNTMIMNSQSRSEAHRLKSKLKELLNKRNALREEENNRLSPAQEREKLINDVRSNNQALASIGKQLKIVEDQLIEKKETLQQIDQDLEEGNSERHVKYKELKKRDDVMSAFMDSFKSNMNQEQQCKCIKLIYSKCLIITQSLLCKHVTGIDTLKNQITYAIEQITMQGINMNGMYDAKLEGNGFTSKNDLNSHSGLMKEYKKLGIQLKQLQILEKRTVQQMNSLRQEETEALQSIQKYANLEIPRSEAIAKMNELTTILQEMEDKKRVTENVVDEARNRNHEIKINLKSNDTYRQISHLEDKLIDLMKDNKVLQETVKNIQQVYNAYQFK';
    my $cdna_2   = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATGCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGC';
    my $coding_2 = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATGCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGC';
    my $aa_2     = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALNHLLHCALELQVGSWQIMEFL*FRRHPNALEQL*EMLEMRMADRPITQHGISGLSTSYGRLGTAVSSNRQIKDKR';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 2, 'run method: two transcripts stored');
    
    my %exon_ids;
    foreach my $transcript (@$transcripts) {
      if ($transcript->stable_id eq 'test_exon_phase-RA') {
        foreach my $exon (@{$transcript->get_all_Exons()}) {
          $exon_ids{$exon->stable_id}++;
        }
        is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
        is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
        is($transcript->seq->seq(), $cdna_1, 'run method: correct cDNA sequence stored');
        is($transcript->translateable_seq(), $coding_1, 'run method: correct coding sequence stored');
        is($transcript->translate->seq(), $aa_1, 'run method: correct amino acid sequence stored');
      }
      if ($transcript->stable_id eq 'test_exon_phase-RB') {
        foreach my $exon (@{$transcript->get_all_Exons()}) {
          $exon_ids{$exon->stable_id}++;
        }
        is($transcript->biotype, 'nontranslating_CDS', 'run method: transcript biotype correct for '.$transcript->stable_id);
        is(scalar(@{$transcript->get_all_Exons()}), 3, 'run method: three exons stored');
        is($transcript->seq->seq(), $cdna_2, 'run method: correct cDNA sequence stored');
        is($transcript->translateable_seq(), $coding_2, 'run method: correct coding sequence stored');
        is($transcript->translate->seq(), $aa_2, 'run method: correct amino acid sequence stored');
      }
    }
    is(scalar(keys(%exon_ids)), 6, 'run method: stored exons with same location and different phase as two distinct exons');
  }
  
  if ($test_gene->stable_id eq 'test_polypeptide_with_CDS') {
    my $cdna   = 'TTCGTCTGAGGTATATTATAGGAGATTCCTAAATCTAAATTCAGTTAGAGTGAATAGTTTAGAGTGAAGCGAATCAGTCGCACAAGTGTTAACATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAGAAAGTAGTGCATACAAACTGGTTAGGTTTGATGTGCAACATTCGTACTAAGTTGATAGTGATTTTGTAAATTTCACTAAGATTTCAACGGTTCCCTGTAATATAAAAGTATGAAATTATAGGCATTACTATATTATCAAATCATTGATTAACAAAATCGATAATTTAATATAAAATATAAACCGTTAAATAAAAGATGGTTGTTT';
    my $coding = 'ATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAG';
    my $aa     = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPL';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 5, 'run method: five exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
    is($transcript->translation->stable_id, 'test_polypeptide_with_CDS-polypeptide', 'run method: translation stable_id correct');
  }
  
  if ($test_gene->stable_id eq 'test_polypeptide_no_CDS') {
    my $cdna   = 'TTCGTCTGAGGTATATTATAGGAGATTCCTAAATCTAAATTCAGTTAGAGTGAATAGTTTAGAGTGAAGCGAATCAGTCGCACAAGTGTTAACATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAGAAAGTAGTGCATACAAACTGGTTAGGTTTGATGTGCAACATTCGTACTAAGTTGATAGTGATTTTGTAAATTTCACTAAGATTTCAACGGTTCCCTGTAATATAAAAGTATGAAATTATAGGCATTACTATATTATCAAATCATTGATTAACAAAATCGATAATTTAATATAAAATATAAACCGTTAAATAAAAGATGGTTGTTT';
    my $coding = 'ATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAG';
    my $aa     = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPL';
    
    is($test_gene->biotype, 'protein_coding', 'run method: gene biotype correct for '.$test_gene->stable_id);
    
    my $transcripts = $test_gene->get_all_Transcripts();
    is(scalar(@{$transcripts}), 1, 'run method: one transcript stored');
    
    my $transcript = $$transcripts[0];
    
    is($transcript->biotype, 'protein_coding', 'run method: transcript biotype correct for '.$transcript->stable_id);
    is(scalar(@{$transcript->get_all_Exons()}), 5, 'run method: five exons stored');
    is($transcript->seq->seq(), $cdna, 'run method: correct cDNA sequence stored');
    is($transcript->translateable_seq(), $coding, 'run method: correct coding sequence stored');
    is($transcript->translate->seq(), $aa, 'run method: correct amino acid sequence stored');
    is($transcript->translation->stable_id, 'test_polypeptide_no_CDS-protein', 'run method: translation stable_id correct');
  }
  
}

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

$lg_obj->param('gff3_file', $gff3_file);
$lg_obj->param('ignore_types', ['region', 'chromosome']);
$lg_obj->param('polypeptides', 0);

done_testing();

sub sort_coding {  
  if ($a->strand == -1) {
    return $b->start <=> $a->start;
  } else {
    return $a->start <=> $b->start;
  }
}

sub sort_genomic {  
  return $a->start <=> $b->start;
}
