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
use Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my %uniprot_db =
  (
    -driver => 'Oracle',
    -host   => 'ora-dlvm5-026.ebi.ac.uk',
    -port   => 1521,
    -user   => 'spselect',
    -pass   => 'spselect',
    -dbname => 'SWPREAD',
  );

my $uniprot_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%uniprot_db);

my $aa   = $dba->get_adaptor('Analysis');
my $dbea = $dba->get_adaptor('DBEntry');
my $ta   = $dba->get_adaptor('Translation');

my $module_name    = 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @super_methods  = qw(external_db_reset external_db_update remove_xrefs);
my @module_methods = qw(add_xrefs add_xref get_go_for_uniprot);
can_ok($module_name, @hive_methods);
can_ok($module_name, @super_methods);
can_ok($module_name, @module_methods);

my $uniprot_edbs = {reviewed   => 'Uniprot/SWISSPROT', unreviewed => 'Uniprot/SPTREMBL'};

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj     = Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtGO->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('db_type'), 'core',                            'param_defaults method: db_type');
is($obj->param('logic_name'), 'xrefuniprot',                  'param_defaults method: logic_name');
is($obj->param('external_db'), 'GO',                          'param_defaults method: external_db');
is_deeply($obj->param('uniprot_external_dbs'), $uniprot_edbs, 'param_defaults method: uniprot_external_dbs');
is($obj->param('replace_all'), 0,                             'param_defaults method: replace_all');

# Mandatory parameters.
$obj->param('species', $species);
$obj->param('uniprot_db', \%uniprot_db);

my $analysis    = $aa->fetch_by_logic_name($obj->param('logic_name'));
my $external_db = $obj->param('external_db');

my $translation = $ta->fetch_by_stable_id('AGAP004700-PA');

my $ac = 'A0NFG9';
my $uniprot_xrefs = $dbea->fetch_all_by_Translation($translation, 'Uniprot/SPTREMBL');
my $uniprot_xref;
foreach my $xref (@$uniprot_xrefs) {
  if ($xref->primary_id eq $ac) {
    $uniprot_xref = $xref;
    last;
  }
}

# get_go_for_uniprot method
my $gos = $obj->get_go_for_uniprot($uniprot_dba, $ac);
ok(scalar(@$gos) >= 1, 'get_go_for_uniprot method: xrefs returned');

# add_xref method
$testdb->hide($dbtype, 'object_xref', 'xref');

foreach my $go (@$gos) {
  my $xref        = $obj->add_xref($go, $uniprot_xref, $analysis, $external_db);
  my $term        = $$go[0];
  my $description = $$go[1];
  my $evidence    = $$go[2];

  isa_ok($xref, 'Bio::EnsEMBL::OntologyXref');
  is($xref->primary_id, $term,         'add_xref method: primary_id set');
  is($xref->display_id, $term,         'add_xref method: display_id set');
  is($xref->description, $description, 'add_xref method: description set');
  is($xref->dbname, 'GO',              'add_xref method: external_db set');
  is($xref->info_type, 'DEPENDENT',    'add_xref method: info_type set');

  my $linkage_info = $xref->get_all_linkage_info();
  is(scalar(@$linkage_info), 1,                        'add_xref method: one linkage');
  is($$linkage_info[0][0], 'IEA',                      'add_xref method: linkage_type set');
  is($$linkage_info[0][1]->primary_id, $ac,            'add_xref method: master primary_id correct');
  is($$linkage_info[0][1]->dbname, 'Uniprot/SPTREMBL', 'add_xref method: master external_db correct');

  my $ignore_release = 1;
  $dbea->store($xref, $translation->dbID(), 'Translation', $ignore_release);

  is_rows(1, $dba, 'xref', 'WHERE dbprimary_acc = ? ', [$term]);
  is_rows(1, $dba,
    'xref INNER JOIN object_xref USING (xref_id)',
    'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND ensembl_object_type = "Translation"',
    [$term, $translation->dbID(), $analysis->dbID()]
  );
  is_rows(1, $dba,
    'xref INNER JOIN object_xref USING (xref_id) INNER JOIN ontology_xref USING (object_xref_id)',
    'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND source_xref_id =? AND linkage_type = ? AND ensembl_object_type = "Translation"',
    [$term, $translation->dbID(), $analysis->dbID(), $uniprot_xref->dbID(), $evidence]
  );
}

$testdb->restore($dbtype, 'object_xref', 'xref');

# run method
my $go_xrefs = object_xrefs();
is($$go_xrefs{IEA}{'Interpro'}, 545,         'run method: correct number of GO xrefs');
is($$go_xrefs{IEA}{'Uniprot/SWISSPROT'}, 63, 'run method: correct number of GO xrefs');
is($$go_xrefs{IEA}{'Uniprot/SPTREMBL'}, 202, 'run method: correct number of GO xrefs');
is($$go_xrefs{IEA}{'null'}, 6,               'run method: correct number of GO xrefs');
is($$go_xrefs{ISS}{'Uniprot/SWISSPROT'}, 2,  'run method: correct number of GO xrefs');
is($$go_xrefs{IBA}{'Uniprot/SWISSPROT'}, 9,  'run method: correct number of GO xrefs');
is($$go_xrefs{IBA}{'Uniprot/SPTREMBL'}, 137, 'run method: correct number of GO xrefs');

$obj->run();

$go_xrefs = object_xrefs();
is($$go_xrefs{IEA}{'Interpro'}, 545,                    'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IEA}{'Uniprot/SWISSPROT'}, '>=', 114, 'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IEA}{'Uniprot/SPTREMBL'}, '>=', 401,  'run method: correct number of GO xrefs');
is($$go_xrefs{IEA}{'null'}, 6,                          'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{ISS}{'Uniprot/SWISSPROT'}, '>=', 4,   'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IBA}{'Uniprot/SWISSPROT'}, '>=', 20,  'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IBA}{'Uniprot/SPTREMBL'}, '>=', 286,  'run method: correct number of GO xrefs');

# Delete all object_xrefs and ontology_xrefs with matching external_db
$obj->param('replace_all', 1);
$obj->run();

$go_xrefs = object_xrefs();
ok(!defined($$go_xrefs{IEA}{'Interpro'}),              'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IEA}{'Uniprot/SWISSPROT'}, '>=', 51, 'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IEA}{'Uniprot/SPTREMBL'}, '>=', 204, 'run method: correct number of GO xrefs');
ok(!defined($$go_xrefs{IEA}{'null'}),                  'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{ISS}{'Uniprot/SWISSPROT'}, '>=', 2,  'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IBA}{'Uniprot/SWISSPROT'}, '>=', 11, 'run method: correct number of GO xrefs');
cmp_ok($$go_xrefs{IBA}{'Uniprot/SPTREMBL'}, '>=', 152, 'run method: correct number of GO xrefs');

done_testing();

sub object_xrefs {
  my %go_xrefs;
  my $translations = $ta->fetch_all();

  foreach my $translation (@$translations) {
    my $go_xrefs = $translation->get_all_DBEntries('GO');
    foreach my $go_xref (@$go_xrefs) {
      my $linkages = $go_xref->get_all_linkage_info();
      foreach my $linkage (@$linkages) {
        my $dbname = defined($$linkage[1]) ? $$linkage[1]->dbname : 'null';
        $go_xrefs{$$linkage[0]}{$dbname}++;
      }
    }
  }

  return \%go_xrefs;
}
