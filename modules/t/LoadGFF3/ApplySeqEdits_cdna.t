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
use Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits;
use Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);
my $ta     = $dba->get_adaptor('Transcript');

my $gff3_file    = $FindBin::Bin.'/../test-files/agam_seqedits.gff3';
my $fasta_file   = $FindBin::Bin.'/../test-files/agam.fa';
my $genbank_file = $FindBin::Bin.'/../test-files/agam.gbff';
my $protein_file = $FindBin::Bin.'/../test-files/agam_protein.fa';
my $cdna_file    = $FindBin::Bin.'/../test-files/agam_cdna.fa';

my $gene_source = 'Ensembl';
my $logic_name  = 'gff3_test';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(
  add_transcript_seq_edit add_translation_seq_edit extract_edits extract_seq 
  load_fasta seq_edits_from_genbank seq_edits_from_cdna seq_edits_from_protein 
  set_protein_coding);
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

$lg_obj->run();

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $ase_obj = Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits->new;
my $ase_job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$ase_obj->input_job($ase_job_obj);

$ase_obj->param('db_type', 'core');
$ase_obj->param('species', $species);

# seq_edits_from_genbank method
#$ase_obj->seq_edits_from_genbank($dba, $genbank_file);
#$ase_obj->set_protein_coding($dba);

# seq_edits_from_protein method
$testdb->hide($dbtype, qw(transcript_attrib translation_attrib));

{

my $transcript = $ta->fetch_by_stable_id('test_protein_edits-RA');

my $coding_seq_before = 'ATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAGAAAGTAGTGCATACAAACTGGTTAGGTTTGATGTGCAACATTCGTACTAAGTTGATAGTGATTTTGTAAATTTCACTAAGATTTCAACGGTTCCCTGTAATATAAAAGTATGAAATTATAGGCATTACTATATTATCAAATCATTGATTAACAAAATCGATAATTTAA';
my $coding_seq_after  = $coding_seq_before;
my $aa_seq_before     = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPL*KVVHTNWLGLMCNIRTKLIVIL*ISLRFQRFPVI*KYEIIGITILSNH*LTKSII';
my $aa_seq_after      = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPLXKVVHTNWLGLMCNIRTKLIVILYISLRFQRFPVIXKYEIIGITILSNHULTKSII';

my @seq_edits = (
  'amino_acid_sub 309 309 X',
  'amino_acid_sub 332 332 Y',
  'amino_acid_sub 344 344 X',
  '_selenocysteine 358 358 U',
);

is($transcript->translateable_seq(), $coding_seq_before, 'seq_edits_from_protein method (before): coding sequence as expected');
is($transcript->translate->seq(), $aa_seq_before, 'seq_edits_from_protein method (before): amino acid sequence as expected');

$ase_obj->seq_edits_from_protein($dba, $protein_file);

# Need to re-fetch from the database.
$transcript     = $ta->fetch_by_stable_id('test_protein_edits-RA');
my $translation = $transcript->translation();
my $gene        = $transcript->get_Gene();

is($transcript->translateable_seq(), $coding_seq_after, 'seq_edits_from_protein method (after): coding sequence as expected');
is($transcript->translate->seq(), $aa_seq_after, 'seq_edits_from_protein method (after): amino acid sequence as expected');

my $transcript_seq_edits = $transcript->get_all_SeqEdits();
is(scalar(@$transcript_seq_edits), 0, 'seq_edits_from_protein method: no trancript-level seq edits');

my $translation_seq_edits = $translation->get_all_SeqEdits();
is(scalar(@$translation_seq_edits), 4, 'seq_edits_from_protein method: four translation-level seq edits');

my @translation_seq_edits = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$translation_seq_edits;
for (my $i=0; $i<@translation_seq_edits; $i++) {
  my $seq_edit = join(' ',
    $translation_seq_edits[$i]->code,
    $translation_seq_edits[$i]->start,
    $translation_seq_edits[$i]->end,
    $translation_seq_edits[$i]->alt_seq,
  );
  is($seq_edit, $seq_edits[$i], "seq_edits_from_protein method: $seq_edit stored correctly");
}

is($gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (before): gene biotype correct');
is($transcript->biotype, 'nontranslating_CDS', 'set_protein_coding method (before): transcript biotype correct');

$ase_obj->set_protein_coding($dba);

# Need to re-fetch from the database.
$transcript = $ta->fetch_by_stable_id('test_protein_edits-RA');
$gene       = $transcript->get_Gene();

is($gene->biotype, 'protein_coding', 'set_protein_coding method (after): gene biotype correct');
is($transcript->biotype, 'protein_coding', 'set_protein_coding method (after): transcript biotype correct');

}

$testdb->restore($dbtype, qw(transcript_attrib translation_attrib));

#  # seq_edits_from_cdna method
#  $testdb->hide($dbtype, qw(transcript_attrib translation_attrib));
#  
#  {
#  
#  my $transcript = $ta->fetch_by_stable_id('test_cdna_edits-RA');
#  
#  my $coding_seq_before = 'TGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAG';
#  my $coding_seq_after  = 'ATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTGTGTAG';
#  my $aa_seq_before     = '*KKV*HF*PVSLSASKSLALKQKFENFLTLVMLIHRITCDIMHYICLIMIIQTDYDLH*KNVPCQICSFR*RYFLPKNHTFVPRCS*ALIFC*HQRCV*NLCNRLMTIPRVICLF*SRRNTYSITKVVAGT*TKFFIILS*RKNQCTII*PL*SCAIQFVKR*WQMVKVLLLVYGLK*SCAITKFTWANGLSIILNKILPFVGLMFQS*RGRNVARSFRKIK*SYQNLIAE*LKHNYVLRTRKIVP*LSFVNLVRQDHYL*PSVTQFTW*GCLQFILMIVMCKLRCLIRFHHFWIGLKQSYGPILNLC';
#  my $aa_seq_after      = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPV';
#  
#  my @seq_edits = (
#    "_rna_edit 1 1225 $coding_seq_after",
#  );
#  
#  print "cDNA coding start ".$transcript->cdna_coding_start."\n";
#  print "coding region start ".$transcript->coding_region_start."\n";
#  
#  is($transcript->translateable_seq(), $coding_seq_before, 'seq_edits_from_cdna method (before): coding sequence as expected');
#  is($transcript->translate->seq(), $aa_seq_before, 'seq_edits_from_cdna method (before): amino acid sequence as expected');
#  
#  $ase_obj->seq_edits_from_cdna($dba, $cdna_file);
#  
#  # Need to re-fetch from the database.
#  $transcript     = $ta->fetch_by_stable_id('test_cdna_edits-RA');
#  my $translation = $transcript->translation();
#  my $gene        = $transcript->get_Gene();
#  print "cDNA coding start ".$transcript->cdna_coding_start."\n";
#  print "coding region start ".$transcript->coding_region_start."\n";
#  
#  print $transcript->seq->seq."\n";
#  is($transcript->translateable_seq(), $coding_seq_after, 'seq_edits_from_cdna method (after): coding sequence as expected');
#  is($transcript->translate->seq(), $aa_seq_after, 'seq_edits_from_cdna method (after): amino acid sequence as expected');
#  
#  my $transcript_seq_edits = $transcript->get_all_SeqEdits();
#  is(scalar(@$transcript_seq_edits), 1, 'seq_edits_from_cdna method: one trancript-level seq edit');
#  
#  my $translation_seq_edits = $translation->get_all_SeqEdits();
#  is(scalar(@$translation_seq_edits), 0, 'seq_edits_from_protein method: no translation-level seq edits');
#  
#  my @transcript_seq_edits = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$transcript_seq_edits;
#  for (my $i=0; $i<@transcript_seq_edits; $i++) {
#    my $seq_edit = join(' ',
#      $transcript_seq_edits[$i]->code,
#      $transcript_seq_edits[$i]->start,
#      $transcript_seq_edits[$i]->end,
#      $transcript_seq_edits[$i]->alt_seq,
#    );
#    is($seq_edit, $seq_edits[$i], 'seq_edits_from_protein method: seq_edit stored correctly');
#  }
#  
#  }
#  
#  $testdb->restore($dbtype, qw(transcript_attrib translation_attrib));

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

#run with genbank, cdna # cdna should not fix things if genbank already has
#
#run with genbank, protein # protein should not fix things if genbank already has
#
#run with cdna, protein # cdna should not fix things if protein already has
#
#run with genbank, cdna, protein # protein should not fix things if genbank already has; cdna should not fix things if genbank or protein already has

done_testing();
