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

my $gff3_file    = $FindBin::Bin.'/../test-files/LoadGFF3/agam_seqedits.gff3';
my $fasta_file   = $FindBin::Bin.'/../test-files/LoadGFF3/agam.fa';
my $genbank_file = $FindBin::Bin.'/../test-files/LoadGFF3/agam_seqedits.gbff';
my $protein_file = $FindBin::Bin.'/../test-files/LoadGFF3/agam_protein.fa';

my $gene_source = 'Ensembl';
my $logic_name  = 'gff3_test';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(
  add_transcript_seq_edit add_translation_seq_edit
  extract_edits extract_seq
  seq_edits_from_genbank seq_edits_from_protein
  load_fasta set_protein_coding);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

$testdb->hide($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

# Need to load some genes which need seq edits
my $lg_obj  = Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadGFF3->new;
my $lg_job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$lg_obj->input_job($lg_job_obj);

$lg_obj->param('db_type',         'core');
$lg_obj->param('species',         $species);
$lg_obj->param('gff3_file',       $gff3_file);
$lg_obj->param('fasta_file',      $fasta_file);
$lg_obj->param('logic_name',      $logic_name);
$lg_obj->param('gene_source',     $gene_source);
$lg_obj->param('gene_types',      ['gene']);
$lg_obj->param('mrna_types',      ['mRNA']);
$lg_obj->param('exon_types',      ['exon']);
$lg_obj->param('cds_types',       ['CDS']);
$lg_obj->param('utr_types',       ['five_prime_UTR', 'three_prime_UTR']);
$lg_obj->param('ignore_types',    ['region', 'chromosome']);
$lg_obj->param('types_complete',  1);
$lg_obj->param('nontranslating',  'nontranslating_CDS');
$lg_obj->param('polypeptides',    0);
$lg_obj->param('prediction',      0);
$lg_obj->param('min_intron_size', 3);
$lg_obj->param('use_name_field',  'no');
$lg_obj->param('stable_ids',      {});

$lg_obj->run();

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $ase_obj = Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits->new;
my $ase_job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$ase_obj->input_job($ase_job_obj);

# Set and check default parameters.
my $param_defaults = $ase_obj->param_defaults();
$ase_obj->input_job->param_init($param_defaults);
is($ase_obj->param('db_type'), 'core', 'param_defaults method: db_type');

# Species and logic_name are mandatory parameters.
$ase_obj->param('species',    $species);
$ase_obj->param('logic_name', $logic_name);

# seq_edits_from_genbank method
$testdb->hide($dbtype, qw(transcript_attrib translation_attrib));

{
my $transcript_1 = $ta->fetch_by_stable_id('test_gbff_neg_ins-RA');
my $transcript_2 = $ta->fetch_by_stable_id('test_gbff_pos_ins-RA');
my $transcript_3 = $ta->fetch_by_stable_id('test_gbff_pos_del-RA');

my $cdna_seq_1_before   = 'AAAAGAAAAGAAAGAAAAAATCCTTGTTGCGATGCATGCGAAGAATGGTACCATGGAGATTGTATTAATGTGTCAGAAAAGGAAGCCAAGCATATCAAACATTACTACTGTCAGCGATGCAAAGAAGAAGACCCTTCGTTACAGACCGTTTTCCGCTTGGTCCCAGCGCCTGGTCCAATCCCGCTTCCGGAGGAAAAAAAACCTAAAAAGAAAAAGAAGGAAACTCCCATTGGTGGTTCATCTGAGAAGCGTTGCGGTAGCTGCGACGGTTGCCTTGCCGAAAATTGTGGTAACTGTGAAGGTTGTTTGGGTCTCACTAAGAACGGCCGAAAGCAGCGTTGCGATATGCGTTTTTGCAGCAACTCAAGTACCCGGAAGAAGGACCGAGCAGCTATAAAAGCGGGAGGAAGAAACAAAAAGCGAAAACGGTCCGTTACTCCAGAGGTTTTTTGTAACCCTGCATTGGAAGGCGAACGACAATGCTACGGTCCCGGGTGTATTATGACCGCCCGACCTCAATCAAAATACTGTTCGGATGAATGTGGCATGAAGTTAGCAACGAGCCGAATCTATCAGGTTTTACCGCAACGTATTCAAGAGTGGTCTCTTAGCCCAGCTGTTGGCGAAGAACAAAACAAGAAGGCACTTGAATCAATTCGCATGAAGCAAGCGATAGTTAGAGCTACTCTAGCAGAACTGGATAAACGACATGCTGAGTTGGACCTGCTGGTTGAACGTGCAAAGAGGTGCACGTTGGATCCAAATGCATTGGAAAATGCTGATTTGGAGGATGAAATGTCAATGTACTGCATAACGTGTGGTCATGAGATTCATTCCAAAACTGCAATTCGGCATATGGAAAAATGTTTCAACAAATACGAGAGTCAGGCAAGCTTTGGCAGCATTTTCAAGACACGAATTGACGGCAATTCAATGTTCTGTGATTTCTACAATCCGGCTAGCAAGACATACTGCAAACGTTTACGAGTGTTATGCCCAGAACACTGCAAGGACCCAAAAATTAATGATACGGACGTCTGTGGTTGCCCTTTAGTACGCAATGTATTTGAATTAACAGGAGAATTTTGCCGAGCACCAAAGAAATCTTGCTTCAAACATTATGTATGGGAAAAAATACGACGTGCAGAAATTGATCTTGAACGAGTACGACAGTGGTTAAAAATGGACGAACTCGTCGAACAGGAACGTCTAATTCGTCAGGCAATGGCATCACGAGCCGGAGTGCTGGGACTAATGCTGCACTCAACGTACAATCATGAAATCATGGAAAAGTTATGTGCTGGTAAATTTTAG';
my $cdna_seq_1_after    = 'CTGAACCAAGTGAATCTACATTGATAGTCCATTCGTGTATTGCGCTCAGCGTAACTAAACCACATTAAAAATCTAATAAGCACACTATGAGTGAACAAAAGAAAAGAAAGAAAAAATCCAAAGAGGAAATAGCTAAGGAATTCGACTTGCCTGAAAGGAAAAGCAAGATTGCCACTATTTACAAACAAGATGGACAAGCTTACTGCTTGTGTAGATCTTCCGATTCTTCGCGATTCATGATTTGTTGCGATGCATGCGAAGAATGGTACCATGGAGATTGTATTAATGTGTCAGAAAAGGAAGCCAAGCATATCAAACATTACTACTGTCAGCGATGCAAAGAAGAAGACCCTTCGTTACAGACCGTTTTCCGCTTGGTCCCAGCGCCTGGTCCAATCCCGCTTCCGGAGGAAAAAAAACCTAAAAAGAAAAAGAAGGAAACTCCCATTGGTGGTTCATCTGAGAAGCGTTGCGGTAGCTGCGACGGTTGCCTTGCCGAAAATTGTGGTAACTGTGAAGGTTGTTTGGGTCTCACTAAGAACGGCCGAAAGCAGCGTTGCGATATGCGTTTTTGCAGCAACTCAAGTACCCGGAAGAAGGACCGAGCAGCTATAAAAGCGGGAGGAAGAAACAAAAAGCGAAAACGGTCCGTTACTCCAGAGGTTTTTTGTAACCCTGCATTGGAAGGCGAACGACAATGCTACGGTCCCGGGTGTATTATGACCGCCCGACCTCAATCAAAATACTGTTCGGATGAATGTGGCATGAAGTTAGCAACGAGCCGAATCTATCAGGTTTTACCGCAACGTATTCAAGAGTGGTCTCTTAGCCCAGCTGTTGGCGAAGAACAAAACAAGAAGGCACTTGAATCAATTCGCATGAAGCAAGCGATAGTTAGAGCTACTCTAGCAGAACTGGATAAACGACATGCTGAGTTGGACCTGCTGGTTGAACGTGCAAAGAGGTGCACGTTGGATCCAAATGCATTGGAAAATGCTGATTTGGAGGATGAAATGTCAATGTACTGCATAACGTGTGGTCATGAGATTCATTCCAAAACTGCAATTCGGCATATGGAAAAATGTTTCAACAAATACGAGAGTCAGGCAAGCTTTGGCAGCATTTTCAAGACACGAATTGACGGCAATTCAATGTTCTGTGATTTCTACAATCCGGCTAGCAAGACATACTGCAAACGTTTACGAGTGTTATGCCCAGAACACTGCAAGGACCCAAAAATTAATGATACGGACGTCTGTGGTTGCCCTTTAGTACGCAATGTATTTGAATTAACAGGAGAATTTTGCCGAGCACCAAAGAAATCTTGCTTCAAACATTATGTATGGGAAAAAATACGACGTGCAGAAATTGATCTTGAACGAGTACGACAGTGGTTAAAAATGGACGAACTCGTCGAACAGGAACGTCTAATTCGTCAGGCAATGGCATCACGAGCCGGAGTGCTGGGACTAATGCTGCACTCAACGTACAATCATGAAATCATGGAAAAGTTATGTGCTGGTAAATTTTAG';
my $coding_seq_1_before = 'AAAGAAAAGAAAGAAAAAATCCTTGTTGCGATGCATGCGAAGAATGGTACCATGGAGATTGTATTAATGTGTCAGAAAAGGAAGCCAAGCATATCAAACATTACTACTGTCAGCGATGCAAAGAAGAAGACCCTTCGTTACAGACCGTTTTCCGCTTGGTCCCAGCGCCTGGTCCAATCCCGCTTCCGGAGGAAAAAAAACCTAAAAAGAAAAAGAAGGAAACTCCCATTGGTGGTTCATCTGAGAAGCGTTGCGGTAGCTGCGACGGTTGCCTTGCCGAAAATTGTGGTAACTGTGAAGGTTGTTTGGGTCTCACTAAGAACGGCCGAAAGCAGCGTTGCGATATGCGTTTTTGCAGCAACTCAAGTACCCGGAAGAAGGACCGAGCAGCTATAAAAGCGGGAGGAAGAAACAAAAAGCGAAAACGGTCCGTTACTCCAGAGGTTTTTTGTAACCCTGCATTGGAAGGCGAACGACAATGCTACGGTCCCGGGTGTATTATGACCGCCCGACCTCAATCAAAATACTGTTCGGATGAATGTGGCATGAAGTTAGCAACGAGCCGAATCTATCAGGTTTTACCGCAACGTATTCAAGAGTGGTCTCTTAGCCCAGCTGTTGGCGAAGAACAAAACAAGAAGGCACTTGAATCAATTCGCATGAAGCAAGCGATAGTTAGAGCTACTCTAGCAGAACTGGATAAACGACATGCTGAGTTGGACCTGCTGGTTGAACGTGCAAAGAGGTGCACGTTGGATCCAAATGCATTGGAAAATGCTGATTTGGAGGATGAAATGTCAATGTACTGCATAACGTGTGGTCATGAGATTCATTCCAAAACTGCAATTCGGCATATGGAAAAATGTTTCAACAAATACGAGAGTCAGGCAAGCTTTGGCAGCATTTTCAAGACACGAATTGACGGCAATTCAATGTTCTGTGATTTCTACAATCCGGCTAGCAAGACATACTGCAAACGTTTACGAGTGTTATGCCCAGAACACTGCAAGGACCCAAAAATTAATGATACGGACGTCTGTGGTTGCCCTTTAGTACGCAATGTATTTGAATTAACAGGAGAATTTTGCCGAGCACCAAAGAAATCTTGCTTCAAACATTATGTATGGGAAAAAATACGACGTGCAGAAATTGATCTTGAACGAGTACGACAGTGGTTAAAAATGGACGAACTCGTCGAACAGGAACGTCTAATTCGTCAGGCAATGGCATCACGAGCCGGAGTGCTGGGACTAATGCTGCACTCAACGTACAATCATGAAATCATGGAAAAGTTATGTGCTGGTAAATTTTAG';
my $coding_seq_1_after  = 'ATGAGTGAACAAAAGAAAAGAAAGAAAAAATCCAAAGAGGAAATAGCTAAGGAATTCGACTTGCCTGAAAGGAAAAGCAAGATTGCCACTATTTACAAACAAGATGGACAAGCTTACTGCTTGTGTAGATCTTCCGATTCTTCGCGATTCATGATTTGTTGCGATGCATGCGAAGAATGGTACCATGGAGATTGTATTAATGTGTCAGAAAAGGAAGCCAAGCATATCAAACATTACTACTGTCAGCGATGCAAAGAAGAAGACCCTTCGTTACAGACCGTTTTCCGCTTGGTCCCAGCGCCTGGTCCAATCCCGCTTCCGGAGGAAAAAAAACCTAAAAAGAAAAAGAAGGAAACTCCCATTGGTGGTTCATCTGAGAAGCGTTGCGGTAGCTGCGACGGTTGCCTTGCCGAAAATTGTGGTAACTGTGAAGGTTGTTTGGGTCTCACTAAGAACGGCCGAAAGCAGCGTTGCGATATGCGTTTTTGCAGCAACTCAAGTACCCGGAAGAAGGACCGAGCAGCTATAAAAGCGGGAGGAAGAAACAAAAAGCGAAAACGGTCCGTTACTCCAGAGGTTTTTTGTAACCCTGCATTGGAAGGCGAACGACAATGCTACGGTCCCGGGTGTATTATGACCGCCCGACCTCAATCAAAATACTGTTCGGATGAATGTGGCATGAAGTTAGCAACGAGCCGAATCTATCAGGTTTTACCGCAACGTATTCAAGAGTGGTCTCTTAGCCCAGCTGTTGGCGAAGAACAAAACAAGAAGGCACTTGAATCAATTCGCATGAAGCAAGCGATAGTTAGAGCTACTCTAGCAGAACTGGATAAACGACATGCTGAGTTGGACCTGCTGGTTGAACGTGCAAAGAGGTGCACGTTGGATCCAAATGCATTGGAAAATGCTGATTTGGAGGATGAAATGTCAATGTACTGCATAACGTGTGGTCATGAGATTCATTCCAAAACTGCAATTCGGCATATGGAAAAATGTTTCAACAAATACGAGAGTCAGGCAAGCTTTGGCAGCATTTTCAAGACACGAATTGACGGCAATTCAATGTTCTGTGATTTCTACAATCCGGCTAGCAAGACATACTGCAAACGTTTACGAGTGTTATGCCCAGAACACTGCAAGGACCCAAAAATTAATGATACGGACGTCTGTGGTTGCCCTTTAGTACGCAATGTATTTGAATTAACAGGAGAATTTTGCCGAGCACCAAAGAAATCTTGCTTCAAACATTATGTATGGGAAAAAATACGACGTGCAGAAATTGATCTTGAACGAGTACGACAGTGGTTAAAAATGGACGAACTCGTCGAACAGGAACGTCTAATTCGTCAGGCAATGGCATCACGAGCCGGAGTGCTGGGACTAATGCTGCACTCAACGTACAATCATGAAATCATGGAAAAGTTATGTGCTGGTAAATTTTAG';
my $aa_seq_1_before     = 'KEKKEKILVAMHAKNGTMEIVLMCQKRKPSISNITTVSDAKKKTLRYRPFSAWSQRLVQSRFRRKKNLKRKRRKLPLVVHLRSVAVAATVALPKIVVTVKVVWVSLRTAESSVAICVFAATQVPGRRTEQL*KREEETKSENGPLLQRFFVTLHWKANDNATVPGVL*PPDLNQNTVRMNVA*S*QRAESIRFYRNVFKSGLLAQLLAKNKTRRHLNQFA*SKR*LELL*QNWINDMLSWTCWLNVQRGARWIQMHWKMLIWRMKCQCTA*RVVMRFIPKLQFGIWKNVSTNTRVRQALAAFSRHELTAIQCSVISTIRLARHTANVYECYAQNTARTQKLMIRTSVVAL*YAMYLN*QENFAEHQRNLASNIMYGKKYDVQKLILNEYDSG*KWTNSSNRNV*FVRQWHHEPECWD*CCTQRTIMKSWKSYVLVNF';
my $aa_seq_1_after      = 'MSEQKKRKKKSKEEIAKEFDLPERKSKIATIYKQDGQAYCLCRSSDSSRFMICCDACEEWYHGDCINVSEKEAKHIKHYYCQRCKEEDPSLQTVFRLVPAPGPIPLPEEKKPKKKKKETPIGGSSEKRCGSCDGCLAENCGNCEGCLGLTKNGRKQRCDMRFCSNSSTRKKDRAAIKAGGRNKKRKRSVTPEVFCNPALEGERQCYGPGCIMTARPQSKYCSDECGMKLATSRIYQVLPQRIQEWSLSPAVGEEQNKKALESIRMKQAIVRATLAELDKRHAELDLLVERAKRCTLDPNALENADLEDEMSMYCITCGHEIHSKTAIRHMEKCFNKYESQASFGSIFKTRIDGNSMFCDFYNPASKTYCKRLRVLCPEHCKDPKINDTDVCGCPLVRNVFELTGEFCRAPKKSCFKHYVWEKIRRAEIDLERVRQWLKMDELVEQERLIRQAMASRAGVLGLMLHSTYNHEIMEKLCAGKF';

my $cdna_seq_2_before   = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAACCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGAAATAAT';
my $cdna_seq_2_after    = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAANCCACCTTCTGCATTGCGCTCTGGAACTGCAAGTAGGCTCTTGGCAAATAATGGAATTTCTATGATTTCGTCGACATCCCAACGCATTGGAACAGCTTTAGGAAATGCTGGAAATCGAATGGCAGACCGTCCAATTACACAGCATGGCATCAGTGGCCTTTCGACGTCTTACGGACGTCTTGGAACTGCTGTTAGTAGTAATCGTCAAATTAAAGACAAACGCTACTGGCAAGCATTATTACAAAGCAAAATTCAAGAAATCAATCAAGAAACAACAAAGATTTTGAAGGAAAAAAAGTTTCTCGATCGTGAACGATCAGCTAGAAAGCTGTACGAGAAACGTGTAAAAGAAGCTGCAAAAGAATTGACAAACCTACAATCCACCCTAACGTCGATGAATTTAGCTTTGGACAACTGTACATCGGGTATGACACGGCAACATCTATTAAACGAAACGGTAGCATTACGAGAACGCAACGAGCACATTCAAGAACAGTTGGAGGTCATTTTCAAGCAGCGACAACAAAAAGACATCGAAAACAAGGCACTGGAACGAAATACAGAACAGGAGAAAAATAAGGTGATTGAAATGATCAATTCGCTACCAGAAGAGGATCAGCACAAATATCGGGAATACAAAGCTTTGTCGGAAAATCTTCGCAAACAGAATGCCATTTACCATAGCCAGATAAGCGAGATGGAAAAACAAAAGGACAGACTTAACACAATGATTATGAACTCACAATCTCGTTCTGAAGCACATCGTCTAAAATCAAAGTTGAAAGAATTGCTAAACAAACGAAATGCTTTACGAGAAGAGGAAAACAATCGCCTCTCGCCAGCTCAGGAGCGCGAAAAACTAATCAACGATGTACGATCCAACAACCAGGCGCTAGCTAGTATTGGCAAGCAACTAAAGATCGTTGAAGATCAATTGATAGAAAAAAAAGAAACATTACAGCAGATCGACCAAGATCTTGAAGAAGGAAACTCCGAGCGCCATGTAAAGTACAAGGAACTTAAAAAACGAGATGATGTAATGTCAGCATTTATGGACAGTTTTAAAAGCAACATGAATCAGGAACAGCAATGTAAGTGCATTAAATTAATTTACTCAAAGTGTCTTATAATAACGCAATCTCTTCTTTGCAAACATGTTACAGGTATAGATACGCTGAAAAATCAAATCACGTATGCTATTGAGCAAATTACAATGCAAGGAATCAACATGAACGGGATGTATGATGCGAAACTAGAAGGAAACGGATTTACATCTAAGAATGATCTTAATTCTCACTCGGGACTGATGAAGGAATACAAAAAACTTGGCATACAACTCAAGCAACTCCAAATATTGGAAAAACGTACAGTGCAGCAAATGAATTCTCTACGCCAGGAAGAAACAGAAGCTCTGCAAAGTATTCAAAAGTATGCTAACTTGGAAATACCTCGATCCGAGGCAATTGCAAAAATGAATGAGCTTACCACCATTCTGCAAGAGATGGAAGATAAAAAACGTGTTACAGAAAATGTGGTCGATGAAGCTCGCAATAGAAATCATGAGATAAAGATTAATCTGAAAAGTAACGATACATATCGACAAATTTCTCATCTTGAAGATAAGCTGATAGATTTAATGAAGGACAATAAAGTGCTTCAAGAAACAGTGAAAAACATACAGCAGGTATATAATGCATACCAGTTTAAATGAAATAAT';
my $coding_seq_2_before = $cdna_seq_2_before;
my $coding_seq_2_after  = $cdna_seq_2_after;
my $aa_seq_2_before     = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALNHLLHCALELQVGSWQIMEFL*FRRHPNALEQL*EMLEIEWQTVQLHSMASVAFRRLTDVLELLLVVIVKLKTNATGKHYYKAKFKKSIKKQQRF*RKKSFSIVNDQLESCTRNV*KKLQKN*QTYNPP*RR*I*LWTTVHRV*HGNIY*TKR*HYENATSTFKNSWRSFSSSDNKKTSKTRHWNEIQNRRKIR*LK*SIRYQKRISTNIGNTKLCRKIFANRMPFTIAR*ARWKNKRTDLTQ*L*THNLVLKHIV*NQS*KNC*TNEMLYEKRKTIASRQLRSAKN*STMYDPTTRR*LVLASN*RSLKIN**KKKKHYSRSTKILKKETPSAM*STRNLKNEMM*CQHLWTVLKAT*IRNSNVSALN*FTQSVL**RNLFFANMLQV*IR*KIKSRMLLSKLQCKEST*TGCMMRN*KETDLHLRMILILTRD**RNTKNLAYNSSNSKYWKNVQCSK*ILYARKKQKLCKVFKSMLTWKYLDPRQLQK*MSLPPFCKRWKIKNVLQKMWSMKLAIEIMR*RLI*KVTIHIDKFLILKIS**I**RTIKCFKKQ*KTYSRYIMHTSLNEI';
my $aa_seq_2_after      = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALXPPSALRSGTASRLLANNGISMISSTSQRIGTALGNAGNRMADRPITQHGISGLSTSYGRLGTAVSSNRQIKDKRYWQALLQSKIQEINQETTKILKEKKFLDRERSARKLYEKRVKEAAKELTNLQSTLTSMNLALDNCTSGMTRQHLLNETVALRERNEHIQEQLEVIFKQRQQKDIENKALERNTEQEKNKVIEMINSLPEEDQHKYREYKALSENLRKQNAIYHSQISEMEKQKDRLNTMIMNSQSRSEAHRLKSKLKELLNKRNALREEENNRLSPAQEREKLINDVRSNNQALASIGKQLKIVEDQLIEKKETLQQIDQDLEEGNSERHVKYKELKKRDDVMSAFMDSFKSNMNQEQQCKCIKLIYSKCLIITQSLLCKHVTGIDTLKNQITYAIEQITMQGINMNGMYDAKLEGNGFTSKNDLNSHSGLMKEYKKLGIQLKQLQILEKRTVQQMNSLRQEETEALQSIQKYANLEIPRSEAIAKMNELTTILQEMEDKKRVTENVVDEARNRNHEIKINLKSNDTYRQISHLEDKLIDLMKDNKVLQETVKNIQQVYNAYQFK*NN';

my $cdna_seq_3_before   = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAAGTGAGTGCTAAGAACTT';
my $cdna_seq_3_after    = 'ATGGAAGATGTTAACAATAACTCATCGGCACGATCCTATTCAGACCGTCCTGCATCCCGAAGAGGTTTAGGTCCAACTAACAATCTATTTGCTTCAAATAACACAACCACCTCAGCAACGCCTAGCCCTGTTGCCATCATAAGACCCGGAACTGCTTTGAATGAGTGCTAAGAACTT';
my $coding_seq_3_before = $cdna_seq_3_before;
my $coding_seq_3_after  = $cdna_seq_3_after;
my $aa_seq_3_before     = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALK*VLRT';
my $aa_seq_3_after      = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALNEC*EL';

my @seq_edits_1 = (
  '_rna_edit 1 0 CTGAACCAAGTGAATCTACATTGATAGTCCATTCGTGTATTGCGCTCAGCGTAACTAAACCACATTAAAAATCTAATAAGCACACTATGAGTGAAC',
  '_rna_edit 24 23 AAAGAGGAAATAGCTAAGGAATTCGACTTGCCTGAAAGGAAAAGCAAGATTGCCACTATTTACAAACAAGATGGACAAGCTTACTGCTTGTGTAGATCTTCCGATTCTTCGCGATTCATGAT',
);

my @seq_edits_2 = (
  '_rna_edit 162 161 N',
);

my @seq_edits_3 = (
  '_rna_edit 162 162 ',
);

is($transcript_1->seq->seq(), $cdna_seq_1_before, 'seq_edits_from_genbank method (before #1): cdna sequence as expected');
is($transcript_1->translateable_seq(), $coding_seq_1_before, 'seq_edits_from_genbank method (before #1): coding sequence as expected');
is($transcript_1->translate->seq(), $aa_seq_1_before, 'seq_edits_from_genbank method (before #1): amino acid sequence as expected');

is($transcript_2->seq->seq(), $cdna_seq_2_before, 'seq_edits_from_genbank method (before #2): cdna sequence as expected');
is($transcript_2->translateable_seq(), $coding_seq_2_before, 'seq_edits_from_genbank method (before #2): coding sequence as expected');
is($transcript_2->translate->seq(), $aa_seq_2_before, 'seq_edits_from_genbank method (before #2): amino acid sequence as expected');

is($transcript_3->seq->seq(), $cdna_seq_3_before, 'seq_edits_from_genbank method (before #3): cdna sequence as expected');
is($transcript_3->translateable_seq(), $coding_seq_3_before, 'seq_edits_from_genbank method (before #3): coding sequence as expected');
is($transcript_3->translate->seq(), $aa_seq_3_before, 'seq_edits_from_genbank method (before #3): amino acid sequence as expected');

$ase_obj->seq_edits_from_genbank($dba, $genbank_file);

# Need to re-fetch from the database.
$transcript_1 = $ta->fetch_by_stable_id('test_gbff_neg_ins-RA');
$transcript_2 = $ta->fetch_by_stable_id('test_gbff_pos_ins-RA');
$transcript_3 = $ta->fetch_by_stable_id('test_gbff_pos_del-RA');

is($transcript_1->seq->seq(), $cdna_seq_1_after, 'seq_edits_from_genbank method (after #1): cdna sequence as expected');
is($transcript_1->translateable_seq(), $coding_seq_1_after, 'seq_edits_from_genbank method (after #1): coding sequence as expected');
is($transcript_1->translate->seq(), $aa_seq_1_after, 'seq_edits_from_genbank method (after #1): amino acid sequence as expected');

is($transcript_2->seq->seq(), $cdna_seq_2_after, 'seq_edits_from_genbank method (after #2): cdna sequence as expected');
is($transcript_2->translateable_seq(), $coding_seq_2_after, 'seq_edits_from_genbank method (after #2): coding sequence as expected');
is($transcript_2->translate->seq(), $aa_seq_2_after, 'seq_edits_from_genbank method (after #2): amino acid sequence as expected');

is($transcript_3->seq->seq(), $cdna_seq_3_after, 'seq_edits_from_genbank method (after #3): cdna sequence as expected');
is($transcript_3->translateable_seq(), $coding_seq_3_after, 'seq_edits_from_genbank method (after #3): coding sequence as expected');
is($transcript_3->translate->seq(), $aa_seq_3_after, 'seq_edits_from_genbank method (after #3): amino acid sequence as expected');

my $transcript_1_seq_edits = $transcript_1->get_all_SeqEdits();
is(scalar(@$transcript_1_seq_edits), 2, 'seq_edits_from_genbank method: two trancript-level seq edits');

my $translation_1_seq_edits = $transcript_1->translation->get_all_SeqEdits();
is(scalar(@$translation_1_seq_edits), 0, 'seq_edits_from_genbank method: no translation-level seq edits');

my @transcript_1_seq_edits = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$transcript_1_seq_edits;
for (my $i=0; $i<@transcript_1_seq_edits; $i++) {
  my $seq_edit = join(' ',
    $transcript_1_seq_edits[$i]->code,
    $transcript_1_seq_edits[$i]->start,
    $transcript_1_seq_edits[$i]->end,
    $transcript_1_seq_edits[$i]->alt_seq,
  );
  is($seq_edit, $seq_edits_1[$i], "seq_edits_from_genbank method: $seq_edit stored correctly");
}

my $transcript_2_seq_edits = $transcript_2->get_all_SeqEdits();
is(scalar(@$transcript_2_seq_edits), 1, 'seq_edits_from_genbank method: one trancript-level seq edit');

my $translation_2_seq_edits = $transcript_2->translation->get_all_SeqEdits();
is(scalar(@$translation_2_seq_edits), 0, 'seq_edits_from_genbank method: no translation-level seq edits');

my @transcript_2_seq_edits = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$transcript_2_seq_edits;
for (my $i=0; $i<@transcript_2_seq_edits; $i++) {
  my $seq_edit = join(' ',
    $transcript_2_seq_edits[$i]->code,
    $transcript_2_seq_edits[$i]->start,
    $transcript_2_seq_edits[$i]->end,
    $transcript_2_seq_edits[$i]->alt_seq,
  );
  is($seq_edit, $seq_edits_2[$i], "seq_edits_from_genbank method: $seq_edit stored correctly");
}

my $transcript_3_seq_edits = $transcript_3->get_all_SeqEdits();
is(scalar(@$transcript_3_seq_edits), 1, 'seq_edits_from_genbank method: one trancript-level seq edit');

my $translation_3_seq_edits = $transcript_3->translation->get_all_SeqEdits();
is(scalar(@$translation_3_seq_edits), 0, 'seq_edits_from_genbank method: no translation-level seq edits');

my @transcript_3_seq_edits = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$transcript_3_seq_edits;
for (my $i=0; $i<@transcript_3_seq_edits; $i++) {
  my $seq_edit = join(' ',
    $transcript_3_seq_edits[$i]->code,
    $transcript_3_seq_edits[$i]->start,
    $transcript_3_seq_edits[$i]->end,
    $transcript_3_seq_edits[$i]->alt_seq,
  );
  is($seq_edit, $seq_edits_3[$i], "seq_edits_from_genbank method: $seq_edit stored correctly");
}

is($transcript_1->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (before #1): gene biotype correct');
is($transcript_1->biotype, 'nontranslating_CDS', 'set_protein_coding method (before #1): transcript biotype correct');

is($transcript_2->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (before #2): gene biotype correct');
is($transcript_2->biotype, 'nontranslating_CDS', 'set_protein_coding method (before #2): transcript biotype correct');

is($transcript_3->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (before #3): gene biotype correct');
is($transcript_3->biotype, 'nontranslating_CDS', 'set_protein_coding method (before #3): transcript biotype correct');

$ase_obj->set_protein_coding($dba, $logic_name);

# Need to re-fetch from the database.
$transcript_1 = $ta->fetch_by_stable_id('test_gbff_neg_ins-RA');
$transcript_2 = $ta->fetch_by_stable_id('test_gbff_pos_ins-RA');
$transcript_3 = $ta->fetch_by_stable_id('test_gbff_pos_del-RA');

is($transcript_1->get_Gene->biotype, 'protein_coding', 'set_protein_coding method (after #1): gene biotype correct');
is($transcript_1->biotype, 'protein_coding', 'set_protein_coding method (after #1): transcript biotype correct');

is($transcript_2->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (after #2): gene biotype correct');
is($transcript_2->biotype, 'nontranslating_CDS', 'set_protein_coding method (after #2): transcript biotype correct');

is($transcript_3->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (after #3): gene biotype correct');
is($transcript_3->biotype, 'nontranslating_CDS', 'set_protein_coding method (after #3): transcript biotype correct');
}

$testdb->restore($dbtype, qw(transcript_attrib translation_attrib));

# seq_edits_from_protein method
$testdb->hide($dbtype, qw(transcript_attrib translation_attrib));

{
my $transcript_1 = $ta->fetch_by_stable_id('test_protein_edits-RA');
my $transcript_2 = $ta->fetch_by_stable_id('test_protein_edits_fail-RA');

my $coding_seq_1_before = 'ATGAAAAAAGGTGTAACATTTCTAGCCTGTTTCATTGTCTGCCTCAAAGTCGCTTGCTCTGAAGCAGAAGTTCGAAAACTTTTTAACATTAGTCATGTTAATTCATCGGATTACATGCGATATCATGCATTACATTTGTTTAATAATGATCATCCAAACCGACTACGACCTGCATTAAAAAAATGTCCCATGTCAAATATGCTCTTTCCGGTGAAGATATTTTCTACCGAAGAACCATACTTTTGTTCCGCGGTGTTCATAAGCGCTGATTTTTTGCTAGCACCAGCGATGTGTCTGAAACTTATGCAACCGGTTGATGACCATCCCTCGAGTCATATGTTTGTTTTAATCGAGGCGGAACACGTATTCTATTACGAAGGTGGTCGCCGGTACATAAACAAAATTTTTTATCATCCTAAGCTAGAGGAAGAACCAGTGTACCATAATCTAGCCGTTGTGAAGTTGCGCAATCCAATTCGTGAAACGGTGATGGCAAATGGTCAAAGTATTGTTGCTTGTCTATGGTCTGAAATAAAGCTGCGCAATAACAAAGTTTACTTGGGCGAATGGTTTAAGTATCATCCTGAACAAAATCCTGCCTTTCGTTGGCTTGATGTTCCAGTCATAACGAGGAAGGAATGTCGCGAGGAGCTTTCGAAAAATAAAGTAATCATACCAGAATTTGATCGCGGAGTAGCTGAAACACAACTATGTGTTAAGGACAAGAAAAATAGTACCATGATTGAGTTTTGTGAACCTCGTTCGTCAGGACCATTATTTATGACCCTCGGTAACACAGTTTACGTGGTAGGGATGCCTACAGTTCATATTGATGATTGTAATGTGCAAATTGAGGTGTTTAATCAGGTTTCATCATTTTTGGATTGGATTGAAGCAATCGTATGGCCCCATTTTGAACCTCTGTAGAAAGTAGTGCATACAAACTGGTTAGGTTTGATGTGCAACATTCGTACTAAGTTGATAGTGATTTTGTAAATTTCACTAAGATTTCAACGGTTCCCTGTAATATAAAAGTATGAAATTATAGGCATTACTATATTATCAAATCATTGATTAACAAAATCGATAATTTAA';
my $coding_seq_1_after  = $coding_seq_1_before;
my $aa_seq_1_before     = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPL*KVVHTNWLGLMCNIRTKLIVIL*ISLRFQRFPVI*KYEIIGITILSNH*LTKSII';
my $aa_seq_1_after      = 'MKKGVTFLACFIVCLKVACSEAEVRKLFNISHVNSSDYMRYHALHLFNNDHPNRLRPALKKCPMSNMLFPVKIFSTEEPYFCSAVFISADFLLAPAMCLKLMQPVDDHPSSHMFVLIEAEHVFYYEGGRRYINKIFYHPKLEEEPVYHNLAVVKLRNPIRETVMANGQSIVACLWSEIKLRNNKVYLGEWFKYHPEQNPAFRWLDVPVITRKECREELSKNKVIIPEFDRGVAETQLCVKDKKNSTMIEFCEPRSSGPLFMTLGNTVYVVGMPTVHIDDCNVQIEVFNQVSSFLDWIEAIVWPHFEPLXKVVHTNWLGLMCNIRTKLIVILYISLRFQRFPVIXKYEIIGITILSNHULTKSII';

my $aa_seq_2_before = 'HEKRCNISSLFHCLPQSRLL*SRSSKTF*H*SC*FIGLHAISCITFV***SSKPTTTCIKKMSHVKYALSGEDIFYRRTILLFRGVHKR*FFASTSDVSETYATG**PSLESYVCFNRGGTRILLRRWSPVHKQNFLSS*ARGRTSVP*SSRCEVAQSNS*NGDGKWSKYCCLSMV*NKAAQ*QSLLGRMV*VSS*TKSCLSLA*CSSHNEEGMSRGAFEK*SNHTRI*SRSS*NTTMC*GQEK*YHD*VL*TSFVRTIIYDPR*HSLRGRDAYSSY**L*CAN*GV*SGFIIFGLD*SNRMAPF*TSVESSAYKLVRFDVQHSY*VDSDFVNFTKISTVPCNIKV*NYRHYYIIKSLINKIDNL';
my $aa_seq_2_after  = $aa_seq_2_before;

my @seq_edits_1 = (
  'amino_acid_sub 309 309 X',
  'amino_acid_sub 332 332 Y',
  'amino_acid_sub 344 344 X',
  '_selenocysteine 358 358 U',
);

is($transcript_1->translateable_seq(), $coding_seq_1_before, 'seq_edits_from_protein method (before): coding sequence as expected');
is($transcript_1->translate->seq(), $aa_seq_1_before, 'seq_edits_from_protein method (before): amino acid sequence as expected');

is($transcript_2->translate->seq(), $aa_seq_2_before, 'seq_edits_from_protein method (before): amino acid sequence as expected');

$ase_obj->seq_edits_from_protein($dba, $logic_name, $protein_file);

# Need to re-fetch from the database.
$transcript_1 = $ta->fetch_by_stable_id('test_protein_edits-RA');
$transcript_2 = $ta->fetch_by_stable_id('test_protein_edits_fail-RA');

is($transcript_1->translateable_seq(), $coding_seq_1_after, 'seq_edits_from_protein method (after): coding sequence as expected');
is($transcript_1->translate->seq(), $aa_seq_1_after, 'seq_edits_from_protein method (after): amino acid sequence as expected');

my $transcript_1_seq_edits = $transcript_1->get_all_SeqEdits();
is(scalar(@$transcript_1_seq_edits), 0, 'seq_edits_from_protein method: no trancript-level seq edits');

my $translation_1_seq_edits = $transcript_1->translation->get_all_SeqEdits();
is(scalar(@$translation_1_seq_edits), 4, 'seq_edits_from_protein method: four translation-level seq edits');

my @translation_1_seq_edits = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$translation_1_seq_edits;
for (my $i=0; $i<@translation_1_seq_edits; $i++) {
  my $seq_edit = join(' ',
    $translation_1_seq_edits[$i]->code,
    $translation_1_seq_edits[$i]->start,
    $translation_1_seq_edits[$i]->end,
    $translation_1_seq_edits[$i]->alt_seq,
  );
  is($seq_edit, $seq_edits_1[$i], "seq_edits_from_protein method: $seq_edit stored correctly");
}

is($transcript_2->translate->seq(), $aa_seq_2_after, 'seq_edits_from_protein method (after): amino acid sequence as expected');

my $transcript_2_seq_edits = $transcript_2->get_all_SeqEdits();
is(scalar(@$transcript_2_seq_edits), 0, 'seq_edits_from_protein method: no trancript-level seq edits');

my $translation_2_seq_edits = $transcript_2->translation->get_all_SeqEdits();
is(scalar(@$translation_2_seq_edits), 0, 'seq_edits_from_protein method: no translation-level seq edits');

# set_protein_coding method
is($transcript_1->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (before): gene biotype correct');
is($transcript_1->biotype, 'nontranslating_CDS', 'set_protein_coding method (before): transcript biotype correct');

$ase_obj->set_protein_coding($dba, $logic_name);

# Need to re-fetch from the database.
$transcript_1 = $ta->fetch_by_stable_id('test_protein_edits-RA');
$transcript_2 = $ta->fetch_by_stable_id('test_protein_edits_fail-RA');

is($transcript_1->get_Gene->biotype, 'protein_coding', 'set_protein_coding method (after): gene biotype correct');
is($transcript_1->biotype, 'protein_coding', 'set_protein_coding method (after): transcript biotype correct');

is($transcript_2->get_Gene->biotype, 'nontranslating_CDS', 'set_protein_coding method (after): gene biotype correct');
is($transcript_2->biotype, 'nontranslating_CDS', 'set_protein_coding method (after): transcript biotype correct');
}

$testdb->restore($dbtype, qw(transcript_attrib translation_attrib));

# run method
$testdb->hide($dbtype, qw(transcript_attrib translation_attrib));

$ase_obj->param('genbank_file', $genbank_file);
$ase_obj->param('protein_fasta_file', $protein_file);

{
my $aa_seq_1 = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALXPPSALRSGTASRLLANNGISMISSTSQRIGTALGNAGNRMADRPITQHGISGLSTSYGRLGTAVSSNRQIKDKRYWQALLQSKIQEINQETTKILKEKKFLDRERSARKLYEKRVKEAAKELTNLQSTLTSMNLALDNCTSGMTRQHLLNETVALRERNEHIQEQLEVIFKQRQQKDIENKALERNTEQEKNKVIEMINSLPEEDQHKYREYKALSENLRKQNAIYHSQISEMEKQKDRLNTMIMNSQSRSEAHRLKSKLKELLNKRNALREEENNRLSPAQEREKLINDVRSNNQALASIGKQLKIVEDQLIEKKETLQQIDQDLEEGNSERHVKYKELKKRDDVMSAFMDSFKSNMNQEQQCKCIKLIYSKCLIITQSLLCKHVTGIDTLKNQITYAIEQITMQGINMNGMYDAKLEGNGFTSKNDLNSHSGLMKEYKKLGIQLKQLQILEKRTVQQMNSLRQEETEALQSIQKYANLEIPRSEAIAKMNELTTILQEMEDKKRVTENVVDEARNRNHEIKINLKSNDTYRQISHLEDKLIDLMKDNKVLQETVKNIQQVYNAYQFKUNN';
my $aa_seq_2 = 'MEDVNNNSSARSYSDRPASRRGLGPTNNLFASNNTTTSATPSPVAIIRPGTALNECXEL';
my $aa_seq_3 = 'MSEQKKRKKKSKEEIAKEFDLPERKSKIATIYKQDGQAYCLCRSSDSSRFMICCDACEEWYHGDCINVSEKEAKHIKHYYCQRCKEEDPSLQTVFRLVPAPGPIPLPEEKKPKKKKKETPIGGSSEKRCGSCDGCLAENCGNCEGCLGLTKNGRKQRCDMRFCSNSSTRKKDRAAIKAGGRNKKRKRSVTPEVFCNPALEGERQCYGPGCIMTARPQSKYCSDECGMKLATSRIYQVLPQRIQEWSLSPAVGEEQNKKALESIRMKQAIVRATLAELDKRHAELDLLVERAKRCTLDPNALENADLEDEMSMYCITCGHEIHSKTAIRHMEKCFNKYESQASFGSIFKTRIDGNSMFCDFYNPASKTYCKRLRVLCPEHCKDPKINDTDVCGCPLVRNVFELTGEFCRAPKKSCFKHYVWEKIRRAEIDLERVRQWLKMDELVEQERLIRQAMASRAGVLGLMLHSTYNHEIMEKLCAGKF';
my $aa_seq_4 = 'MXXD';
my $aa_seq_5 = 'MKC';
my $aa_seq_6 = 'MXN';
my $aa_seq_7 = 'MXM';
my $aa_seq_8 = 'KXXXE';
my $aa_seq_9 = 'KKE';
my $aa_seq_10 = 'KXKE';
my $aa_seq_11 = 'EXE';

$ase_obj->run();

my $transcript_1  = $ta->fetch_by_stable_id('test_gbff_pos_ins-RA');
my $transcript_2  = $ta->fetch_by_stable_id('test_gbff_pos_del-RA');
my $transcript_3  = $ta->fetch_by_stable_id('test_gbff_neg_ins-RA');
my $transcript_4  = $ta->fetch_by_stable_id('test_gbff_pos_ins_ins-RA');
my $transcript_5  = $ta->fetch_by_stable_id('test_gbff_pos_del_del-RA');
my $transcript_6  = $ta->fetch_by_stable_id('test_gbff_pos_ins_del-RA');
my $transcript_7  = $ta->fetch_by_stable_id('test_gbff_pos_del_ins-RA');
my $transcript_8  = $ta->fetch_by_stable_id('test_gbff_neg_ins_ins-RA');
my $transcript_9  = $ta->fetch_by_stable_id('test_gbff_neg_del_del-RA');
my $transcript_10  = $ta->fetch_by_stable_id('test_gbff_neg_ins_del-RA');
my $transcript_11  = $ta->fetch_by_stable_id('test_gbff_neg_del_ins-RA');

is($transcript_1->translate->seq(), $aa_seq_1, 'run method: amino acid sequence as expected');
is($transcript_2->translate->seq(), $aa_seq_2, 'run method: amino acid sequence as expected');
is($transcript_3->translate->seq(), $aa_seq_3, 'run method: amino acid sequence as expected');
is($transcript_4->translate->seq(), $aa_seq_4, 'run method: amino acid sequence as expected');
is($transcript_5->translate->seq(), $aa_seq_5, 'run method: amino acid sequence as expected');
is($transcript_6->translate->seq(), $aa_seq_6, 'run method: amino acid sequence as expected');
is($transcript_7->translate->seq(), $aa_seq_7, 'run method: amino acid sequence as expected');
is($transcript_8->translate->seq(), $aa_seq_8, 'run method: amino acid sequence as expected');
is($transcript_9->translate->seq(), $aa_seq_9, 'run method: amino acid sequence as expected');
is($transcript_10->translate->seq(), $aa_seq_10, 'run method: amino acid sequence as expected');
is($transcript_11->translate->seq(), $aa_seq_11, 'run method: amino acid sequence as expected');

my $transcript_1_seq_edits = $transcript_1->get_all_SeqEdits();
is(scalar(@$transcript_1_seq_edits), 1, 'run method: one trancript-level seq edit');

my $translation_1_seq_edits = $transcript_1->translation->get_all_SeqEdits();
is(scalar(@$translation_1_seq_edits), 1, 'run method: one translation-level seq edit');

my $transcript_2_seq_edits = $transcript_2->get_all_SeqEdits();
is(scalar(@$transcript_2_seq_edits), 1, 'run method: one trancript-level seq edit');

my $translation_2_seq_edits = $transcript_2->translation->get_all_SeqEdits();
is(scalar(@$translation_2_seq_edits), 1, 'run method: one translation-level seq edit');

my $transcript_3_seq_edits = $transcript_3->get_all_SeqEdits();
is(scalar(@$transcript_3_seq_edits), 2, 'run method: two trancript-level seq edits');

my $translation_3_seq_edits = $transcript_3->translation->get_all_SeqEdits();
is(scalar(@$translation_3_seq_edits), 0, 'run method: no translation-level seq edits');

my $transcript_4_seq_edits = $transcript_4->get_all_SeqEdits();
is(scalar(@$transcript_4_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_5_seq_edits = $transcript_5->get_all_SeqEdits();
is(scalar(@$transcript_5_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_6_seq_edits = $transcript_6->get_all_SeqEdits();
is(scalar(@$transcript_6_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_7_seq_edits = $transcript_7->get_all_SeqEdits();
is(scalar(@$transcript_7_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_8_seq_edits = $transcript_8->get_all_SeqEdits();
is(scalar(@$transcript_8_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_9_seq_edits = $transcript_9->get_all_SeqEdits();
is(scalar(@$transcript_9_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_10_seq_edits = $transcript_10->get_all_SeqEdits();
is(scalar(@$transcript_10_seq_edits), 2, 'run method: two trancript-level seq edits');

my $transcript_11_seq_edits = $transcript_11->get_all_SeqEdits();
is(scalar(@$transcript_11_seq_edits), 2, 'run method: two trancript-level seq edits');

$ase_obj->set_protein_coding($dba, $logic_name);

# Need to re-fetch from the database.
$transcript_1 = $ta->fetch_by_stable_id('test_gbff_pos_ins-RA');
$transcript_2 = $ta->fetch_by_stable_id('test_gbff_pos_del-RA');
$transcript_3 = $ta->fetch_by_stable_id('test_gbff_neg_ins-RA');
$transcript_4 = $ta->fetch_by_stable_id('test_gbff_pos_ins_ins-RA');
$transcript_5 = $ta->fetch_by_stable_id('test_gbff_pos_del_del-RA');
$transcript_6 = $ta->fetch_by_stable_id('test_gbff_pos_ins_del-RA');
$transcript_7 = $ta->fetch_by_stable_id('test_gbff_pos_del_ins-RA');
$transcript_8 = $ta->fetch_by_stable_id('test_gbff_neg_ins_ins-RA');
$transcript_9 = $ta->fetch_by_stable_id('test_gbff_neg_del_del-RA');
$transcript_10 = $ta->fetch_by_stable_id('test_gbff_neg_ins_del-RA');
$transcript_11 = $ta->fetch_by_stable_id('test_gbff_neg_del_ins-RA');

is($transcript_1->get_Gene->biotype, 'protein_coding', 'run method: gene biotype correct');
is($transcript_1->biotype, 'protein_coding', 'run method: transcript biotype correct');

is($transcript_2->get_Gene->biotype, 'protein_coding', 'run method: gene biotype correct');
is($transcript_2->biotype, 'protein_coding', 'run method: transcript biotype correct');

is($transcript_3->get_Gene->biotype, 'protein_coding', 'run method: gene biotype correct');
is($transcript_3->biotype, 'protein_coding', 'run method: transcript biotype correct');

}

$testdb->restore($dbtype, qw(transcript_attrib translation_attrib));

$testdb->restore($dbtype, 
  qw(exon exon_transcript gene meta_coord transcript translation));

done_testing();
