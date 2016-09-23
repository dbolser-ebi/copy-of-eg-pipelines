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

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Hive::AnalysisJob;
use Bio::EnsEMBL::EGPipeline::Xref::LoadUniParc;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my %uniparc_db =
  (
    -host   => 'mysql-eg-pan-prod.ebi.ac.uk',
    -port   => 4276,
    -user   => 'ensro',
    -dbname => 'uniparc',
  );
my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%uniparc_db);

my $aa   = $dba->get_adaptor('Analysis');
my $dbea = $dba->get_adaptor('DBEntry');
my $ta   = $dba->get_adaptor('Translation');

my $module_name    = 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniParc';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @super_methods  = qw(external_db_reset external_db_update remove_xrefs);
my @module_methods = qw(add_xrefs add_xref search_for_upi md5_checksum);
can_ok($module_name, @hive_methods);
can_ok($module_name, @super_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj     = Bio::EnsEMBL::EGPipeline::Xref::LoadUniParc->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# These are the modules param_defaults.
$obj->param('species', $species);
$obj->param('db_type', 'core');
$obj->param('logic_name', 'xrefchecksum');
$obj->param('external_db', 'UniParc');
$obj->param('uniparc_db', \%uniparc_db);

my $analysis    = $aa->fetch_by_logic_name($obj->param('logic_name'));
my $external_db = $obj->param('external_db');

my $translation = $ta->fetch_by_stable_id('AGAP004700-PA');

# search_for_upi method
my $upi = $obj->search_for_upi($uniparc_dba, $translation);
is($upi, 'UPI000069B0ED', 'search_for_upi method: correct UPI value');

$translation->seq($translation->seq.'XXXXXXXXXX');
my $upi_invalid_seq = $obj->search_for_upi($uniparc_dba, $translation);
ok(!defined($upi_invalid_seq), 'search_for_upi method: undefined UPI value');

# add_xref method
$testdb->hide($dbtype, 'object_xref', 'xref');

my $xref = $obj->add_xref($upi, $analysis, $external_db);
isa_ok($xref, 'Bio::EnsEMBL::DBEntry');
is($xref->primary_id, $upi,         'add_xref method: primary_id set');
is($xref->display_id, $upi,         'add_xref method: display_id set');
is($xref->dbname,     $external_db, 'add_xref method: external_db set');
is($xref->info_type, 'CHECKSUM',    'add_xref method: info_type set');

# storing the xref
my $ignore_release = 1;
$dbea->store($xref, $translation->dbID(), 'Translation', $ignore_release);

is_rows(1, $dba, 'xref', 'WHERE dbprimary_acc = ?', [$upi]);
is_rows(1, $dba,
  'xref INNER JOIN object_xref USING (xref_id)',
  'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND ensembl_object_type = "Translation"',
  [$upi, $translation->dbID(), $analysis->dbID()]
);

$testdb->restore($dbtype, 'object_xref', 'xref');

# run method
my $uniparcs = object_xrefs();
is(scalar(@$uniparcs), 158, 'run method: correct number of UniParc xrefs');

$obj->run();

$uniparcs = object_xrefs();
is(scalar(@$uniparcs), 160, 'run method: correct number of UniParc xrefs');

done_testing();

sub object_xrefs {
  my @uniparcs;
  my $translations = $ta->fetch_all();

  foreach my $translation (@$translations) {
    push @uniparcs, @{$translation->get_all_DBEntries('UniParc')};
  }

  return \@uniparcs;
}
