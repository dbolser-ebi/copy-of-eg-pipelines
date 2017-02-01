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
use Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtXrefs;

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
my $ga   = $dba->get_adaptor('Gene');
my $ta   = $dba->get_adaptor('Translation');

my $module_name    = 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtXrefs';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @super_methods  = qw(external_db_reset external_db_update remove_xrefs);
my @module_methods = qw(add_xrefs add_xref get_xrefs_for_uniprot);
can_ok($module_name, @hive_methods);
can_ok($module_name, @super_methods);
can_ok($module_name, @module_methods);

my $uniprot_edbs = {reviewed => 'Uniprot/SWISSPROT', unreviewed => 'Uniprot/SPTREMBL'};

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj     = Bio::EnsEMBL::EGPipeline::Xref::LoadUniProtXrefs->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('db_type'), 'core',                            'param_defaults method: db_type');
is($obj->param('logic_name'), 'xrefuniprot',                  'param_defaults method: logic_name');
is_deeply($obj->param('external_dbs'), {},                    'param_defaults method: external_dbs');
is_deeply($obj->param('uniprot_external_dbs'), $uniprot_edbs, 'param_defaults method: uniprot_external_dbs');
is($obj->param('replace_all'), 0,                             'param_defaults method: replace_all');

# Mandatory parameters.
$obj->param('species', $species);
$obj->param('uniprot_db', \%uniprot_db);
$obj->param('external_dbs', {});
$obj->param('uniprot_external_dbs', {reviewed => 'Uniprot/SWISSPROT', unreviewed => 'Uniprot/SPTREMBL'});

my $analysis     = $aa->fetch_by_logic_name($obj->param('logic_name'));
my $external_dbs = $obj->param('external_dbs');

my $translation = $ta->fetch_by_stable_id('AGAP004700-PA');

my $ac = 'A0NFG9';

# get_xrefs_for_uniprot method
my $transitive_xrefs = $obj->get_xrefs_for_uniprot($uniprot_dba, $ac);
ok(scalar(@$transitive_xrefs) >= 1, 'get_xrefs_for_uniprot method: xrefs returned');

# add_xref method
$testdb->hide($dbtype, 'object_xref', 'xref');

$external_dbs = {
  ArrayExpress => 'ArrayExpress',
  EMBL         => 'EMBL',
};

foreach my $transitive_xref (@$transitive_xrefs) {
  my $xref         = $obj->add_xref($transitive_xref, $analysis, $external_dbs);
  my $xref_source  = $$transitive_xref[0];
  my $primary_id   = $$transitive_xref[1];
  my $secondary_id = $$transitive_xref[2];

  if ($xref_source eq 'EMBL') {
    ok(defined($xref), "add_xref method: xref returned for $xref_source");
    isa_ok($xref, 'Bio::EnsEMBL::DBEntry');
    is($xref->primary_id, $secondary_id, 'add_xref method: primary_id set');
    is($xref->display_id, $secondary_id, 'add_xref method: display_id set');
    ok(!defined($xref->description),     'add_xref method: no description');
    is($xref->dbname, 'protein_id',      'add_xref method: external_db set');
    is($xref->info_type, 'DEPENDENT',    'add_xref method: info_type set');

    my $ignore_release = 1;
    $dbea->store($xref, $translation->dbID(), 'Translation', $ignore_release);

    is_rows(1, $dba, 'xref', 'WHERE dbprimary_acc = ? ', [$secondary_id]);
    is_rows(1, $dba,
      'xref INNER JOIN object_xref USING (xref_id)',
      'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND ensembl_object_type = "Translation"',
      [$secondary_id, $translation->dbID(), $analysis->dbID()]
    );
  } else {
    ok(!defined($xref), "add_xref method: xref not returned for $xref_source");
  }
}

$external_dbs = {
  GeneID => 'EntrezGene',
  RefSeq => 'RefSeq_peptide',
};

foreach my $transitive_xref (@$transitive_xrefs) {
  my $xref         = $obj->add_xref($transitive_xref, $analysis, $external_dbs);
  my $xref_source  = $$transitive_xref[0];
  my $primary_id   = $$transitive_xref[1];
  my $secondary_id = $$transitive_xref[2];

  if ($xref_source eq 'GeneID') {
    my $db_name = $$external_dbs{$xref_source};

    ok(defined($xref), "add_xref method: xref returned for $xref_source");
    isa_ok($xref, 'Bio::EnsEMBL::DBEntry');
    is($xref->primary_id, $primary_id, 'add_xref method: primary_id set');
    is($xref->display_id, $primary_id, 'add_xref method: display_id set');
    ok(!defined($xref->description),   'add_xref method: no description');
    is($xref->dbname, $db_name,        'add_xref method: external_db set');
    is($xref->info_type, 'DEPENDENT',  'add_xref method: info_type set');

    my $ignore_release = 1;
    $dbea->store($xref, $translation->transcript->get_Gene->dbID(), 'Gene', $ignore_release);

    is_rows(1, $dba, 'xref', 'WHERE dbprimary_acc = ? ', [$primary_id]);
    is_rows(1, $dba,
      'xref INNER JOIN object_xref USING (xref_id)',
      'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND ensembl_object_type = "Gene"',
      [$primary_id, $translation->transcript->get_Gene->dbID(), $analysis->dbID()]
    );
  }

  if ($xref_source eq 'RefSeq') {
    my $db_name = $$external_dbs{$xref_source};

    ok(defined($xref), "add_xref method: xref returned for $xref_source");
    isa_ok($xref, 'Bio::EnsEMBL::DBEntry');
    is($xref->primary_id, $primary_id, 'add_xref method: primary_id set');
    is($xref->display_id, $primary_id, 'add_xref method: display_id set');
    ok(!defined($xref->description),   'add_xref method: no description');
    is($xref->dbname, $db_name,        'add_xref method: external_db set');
    is($xref->info_type, 'DEPENDENT',  'add_xref method: info_type set');

    my $ignore_release = 1;
    $dbea->store($xref, $translation->dbID(), 'Translation', $ignore_release);

    is_rows(1, $dba, 'xref', 'WHERE dbprimary_acc = ? ', [$primary_id]);
    is_rows(1, $dba,
      'xref INNER JOIN object_xref USING (xref_id)',
      'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND ensembl_object_type = "Translation"',
      [$primary_id, $translation->dbID(), $analysis->dbID()]
    );
  }
}

$testdb->restore($dbtype, 'object_xref', 'xref');

# run method
$obj->param('external_dbs',
  {
    ArrayExpress => 'ArrayExpress',
    EMBL         => 'EMBL',
    GeneID       => 'EntrezGene',
    MEROPS       => 'MEROPS',
    RefSeq       => 'RefSeq_peptide',
  }
);

my $xrefs = object_xrefs();

is($$xrefs{xrefuniprot}{'EMBL'}, 7,                      'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefuniprot}{'EntrezGene'}),         'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'MEROPS'}, 9,                    'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'protein_id'}, 306,              'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefuniprot}{'RefSeq_peptide'}),     'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'STRING'}, 102,                  'run method: correct number of xrefs');

is($$xrefs{xrefexonerateprotein}{'EMBL'}, 50,            'run method: correct number of xrefs');
is($$xrefs{xrefexoneratedna}{'EntrezGene'}, 132,         'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'MEROPS'}, 1,           'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'protein_id'}, 51,      'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'RefSeq_peptide'}, 157, 'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'STRING'}),    'run method: correct number of xrefs');

is_rows(0, $dba,
  'dependent_xref INNER JOIN object_xref USING (object_xref_id)',
  'WHERE analysis_id = ?',
  [$analysis->dbID()]
);

$obj->run();

$xrefs = object_xrefs();

is($$xrefs{xrefuniprot}{'EMBL'}, 10,                     'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'EntrezGene'}, 10,               'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'MEROPS'}, 14,                   'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'protein_id'}, 513,              'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'RefSeq_peptide'}, 174,          'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'STRING'}, 102,                  'run method: correct number of xrefs');

is($$xrefs{xrefexoneratedna}{'EntrezGene'}, 132,         'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'EMBL'}, 50,            'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'MEROPS'}, 1,           'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'protein_id'}, 51,      'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'STRING'}),    'run method: correct number of xrefs');
is($$xrefs{xrefexonerateprotein}{'RefSeq_peptide'}, 157, 'run method: correct number of xrefs');

is_rows(524, $dba,
  'dependent_xref INNER JOIN object_xref USING (object_xref_id)',
  'WHERE analysis_id = ?',
  [$analysis->dbID()]
);

# Delete all object_xrefs and ontology_xrefs with matching external_db
$obj->param('replace_all', 1);
$obj->run();

$xrefs = object_xrefs();

is($$xrefs{xrefuniprot}{'EMBL'}, 3,                           'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'EntrezGene'}, 135,                   'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'MEROPS'}, 5,                         'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'protein_id'}, 207,                   'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'RefSeq_peptide'}, 174,               'run method: correct number of xrefs');
is($$xrefs{xrefuniprot}{'STRING'}, 102,                       'run method: correct number of xrefs');

ok(!defined($$xrefs{xrefexoneratedna}{'EntrezGene'}),         'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'EMBL'}),           'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'MEROPS'}),         'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'protein_id'}),     'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'STRING'}),         'run method: correct number of xrefs');
ok(!defined($$xrefs{xrefexonerateprotein}{'RefSeq_peptide'}), 'run method: correct number of xrefs');

is_rows(524, $dba,
  'dependent_xref INNER JOIN object_xref USING (object_xref_id)',
  'WHERE analysis_id = ?',
  [$analysis->dbID()]
);

done_testing();

sub object_xrefs {
  my %xrefs;
  my $genes = $ga->fetch_all();

  foreach my $gene (@$genes) {
    my $xrefs = $gene->get_all_DBLinks();
    foreach my $xref (@$xrefs) {
      my $logic_name = defined($xref->analysis) ? $xref->analysis->logic_name : 'null';
      $xrefs{$logic_name}{$xref->dbname}++;
    }
  }

  return \%xrefs;
}
