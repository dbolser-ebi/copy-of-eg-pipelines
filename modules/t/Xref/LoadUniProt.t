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
use Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my %uniparc_db =
  (
    -driver => 'Oracle',
    -host   => 'ora-vm-004.ebi.ac.uk',
    -port   => 1551,
    -user   => 'uniparc_read',
    -pass   => 'uniparc',
    -dbname => 'UAPRO',
  );

my %uniprot_db =
  (
    -driver => 'Oracle',
    -host   => 'ora-dlvm5-026.ebi.ac.uk',
    -port   => 1521,
    -user   => 'spselect',
    -pass   => 'spselect',
    -dbname => 'SWPREAD',
  );

my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%uniparc_db);
my $uniprot_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%uniprot_db);

my $aa   = $dba->get_adaptor('Analysis');
my $dbea = $dba->get_adaptor('DBEntry');
my $ga   = $dba->get_adaptor('Gene');
my $ta   = $dba->get_adaptor('Translation');

my $module_name    = 'Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @super_methods  = qw(external_db_reset external_db_update remove_xrefs);
my @module_methods = qw(add_xrefs add_xref get_uniprot_for_upi
                        set_descriptions remove_descriptions
                        set_gene_names remove_gene_names);
can_ok($module_name, @hive_methods);
can_ok($module_name, @super_methods);
can_ok($module_name, @module_methods);

my $external_dbs = {reviewed => 'Uniprot/SWISSPROT', unreviewed => 'Uniprot/SPTREMBL'};

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $obj     = Bio::EnsEMBL::EGPipeline::Xref::LoadUniProt->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $obj->param_defaults();
$obj->input_job->param_init($param_defaults);
is($obj->param('db_type'), 'core',                      'param_defaults method: db_type');
is($obj->param('logic_name'), 'xrefuniparc',            'param_defaults method: logic_name');
is_deeply($obj->param('external_dbs'), $external_dbs,   'param_defaults method: external_dbs');
is($obj->param('replace_all'), 0,                       'param_defaults method: replace_all');
is_deeply($obj->param('description_source'), [],        'param_defaults method: description_source');
is($obj->param('overwrite_description'), 0,             'param_defaults method: overwrite_description');
is_deeply($obj->param('description_blacklist'), [],     'param_defaults method: description_blacklist');
is_deeply($obj->param('gene_name_source'), [],          'param_defaults method: gene_name_source');
is($obj->param('overwrite_gene_name'), 0,               'param_defaults method: overwrite_gene_name');
is($obj->param('uniparc_external_db'), 'UniParc',       'param_defaults method: uniparc_external_db');
is($obj->param('uniprot_gn_external_db'), 'Uniprot_gn', 'param_defaults method: uniprot_gn_external_db');

# Mandatory parameters.
$obj->param('species', $species);
$obj->param('uniparc_db', \%uniparc_db);
$obj->param('uniprot_db', \%uniprot_db);

my $analysis = $aa->fetch_by_logic_name($obj->param('logic_name'));

my $tax_id = $dba->get_MetaContainer->get_taxonomy_id();

# Unreviewed UniProt record with no name or (valid) description
{
my $upi = 'UPI000069B0ED';
my $translation = $ta->fetch_by_stable_id('AGAP004700-PA');

# get_uniprot_for_upi method
my $uniprots = $obj->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi);
is(scalar(@$uniprots), 1, 'get_uniprot_for_upi method: one uniprot record');

my %uniprot = %{$$uniprots[0]};
my @synonyms = @{$uniprot{synonyms}};
is($uniprot{ac}, 'A0NFG9', 'get_uniprot_for_upi method: correct accession');
is($uniprot{name}, 'A0NFG9_ANOGA', 'get_uniprot_for_upi method: correct name');
is($uniprot{type}, 'unreviewed', 'get_uniprot_for_upi method: correct accession');
is($uniprot{description}, 'AGAP004700-PA', 'get_uniprot_for_upi method: correct description');
is($uniprot{version}, 1, 'get_uniprot_for_upi method: correct version');
ok(!defined($uniprot{gene_name}), 'get_uniprot_for_upi method: no gene name');
is(scalar(@synonyms), 1, 'get_uniprot_for_upi method: one synonym');
is($synonyms[0], 'AgaP_AGAP004700', 'get_uniprot_for_upi method: correct synonym');

# Deleting the description is done in the add_xrefs method
if ($uniprot{description} eq $translation->stable_id) {
  delete $uniprot{description};
}

## add_xref method
$testdb->hide($dbtype, 'dependent_xref', 'object_xref', 'xref');

my $xref = $obj->add_xref(\%uniprot, $analysis, $external_dbs);
my $external_db = $$external_dbs{'unreviewed'};

isa_ok($xref, 'Bio::EnsEMBL::DBEntry');
is($xref->primary_id,  $uniprot{ac},      'add_xref method: primary_id set');
is($xref->display_id,  $uniprot{name},    'add_xref method: display_id set');
ok(!defined($xref->description),          'add_xref method: no description');
is($xref->version,     $uniprot{version}, 'add_xref method: version set');
is($xref->dbname,      $external_db,      'add_xref method: external_db set');
is($xref->info_type,  'DEPENDENT',        'add_xref method: info_type set');

## storing the xref
my $ignore_release = 1;
$dbea->store($xref, $translation->dbID(), 'Translation', $ignore_release);

is_rows(1, $dba, 'xref', 'WHERE dbprimary_acc = ? ', [$uniprot{ac}]);
is_rows(1, $dba,
  'xref INNER JOIN object_xref USING (xref_id)',
  'WHERE dbprimary_acc = ? AND ensembl_id = ? AND analysis_id = ? AND ensembl_object_type = "Translation"',
  [$uniprot{ac}, $translation->dbID(), $analysis->dbID()]
);

$testdb->restore($dbtype, 'dependent_xref', 'object_xref', 'xref');
}

# Unreviewed UniProt record with name, no description
{
my $upi = 'UPI0001538F73';
my $translation = $ta->fetch_by_stable_id('AGAP004716-PA');

# get_uniprot_for_upi method
my $uniprots = $obj->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi);
is(scalar(@$uniprots), 1, 'get_uniprot_for_upi method: one uniprot record');

my %uniprot = %{$$uniprots[0]};
my @synonyms = @{$uniprot{synonyms}};
is($uniprot{ac}, 'A7UUW5', 'get_uniprot_for_upi method: correct accession');
is($uniprot{name}, 'A7UUW5_ANOGA', 'get_uniprot_for_upi method: correct name');
is($uniprot{type}, 'unreviewed', 'get_uniprot_for_upi method: correct accession');
is($uniprot{description}, 'Gustatory receptor', 'get_uniprot_for_upi method: correct description');
is($uniprot{version}, 1, 'get_uniprot_for_upi method: correct version');
is($uniprot{gene_name}, 'GPRGR57', 'get_uniprot_for_upi method: correct gene name');
is(scalar(@synonyms), 1, 'get_uniprot_for_upi method: one synonym');
is($synonyms[0], 'AgaP_AGAP004716', 'get_uniprot_for_upi method: correct synonym');
## add_xref method ignores unreviewed result if not in external_dbs
my $filtered_external_dbs = {reviewed => 'Uniprot/SWISSPROT'};

my $xref = $obj->add_xref(\%uniprot, $analysis, $filtered_external_dbs);
ok(!defined($xref), 'add_xref method: no result for unreviewed data');
}

# Description blacklisting
{
my $upi = 'UPI0001538F8D';

$obj->param('description_blacklist', ['Eukaryotic']);
my $uniprots = $obj->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi);
my %uniprot = %{$$uniprots[0]};
is($uniprot{description}, 'Eukaryotic translation initiation factor 3 subunit C', 'description blacklist: no partial match');

$obj->param('description_blacklist', ['Eukaryotic.*']);
$uniprots = $obj->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi);
%uniprot = %{$$uniprots[0]};
ok(!defined($uniprot{description}), 'description blacklist: regex match');

$obj->param('description_blacklist', ['Uncharacterized protein', 'Eukaryotic.*']);
$uniprots = $obj->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi);
%uniprot = %{$$uniprots[0]};
ok(!defined($uniprot{description}), 'description blacklist: regex match');

$obj->param('description_blacklist', []);
}

# Reviewed UniProt record with name and description
{
my $upi = 'UPI0001538F8D';
my $translation = $ta->fetch_by_stable_id('AGAP004725-PA');

# get_uniprot_for_upi method
my $uniprots = $obj->get_uniprot_for_upi($uniparc_dba, $uniprot_dba, $tax_id, $upi);
is(scalar(@$uniprots), 1, 'get_uniprot_for_upi method: one uniprot record');

my %uniprot = %{$$uniprots[0]};
my @synonyms = @{$uniprot{synonyms}};
is($uniprot{ac}, 'Q7PMU8', 'get_uniprot_for_upi method: correct accession');
is($uniprot{name}, 'EIF3C_ANOGA', 'get_uniprot_for_upi method: correct name');
is($uniprot{type}, 'reviewed', 'get_uniprot_for_upi method: correct accession');
is($uniprot{description}, 'Eukaryotic translation initiation factor 3 subunit C', 'get_uniprot_for_upi method: correct description');
is($uniprot{version}, 3, 'get_uniprot_for_upi method: correct version');
is($uniprot{gene_name}, 'eIF3-S8', 'get_uniprot_for_upi method: no gene name');
is(scalar(@synonyms), 1, 'get_uniprot_for_upi method: one synonym');
is($synonyms[0], 'AGAP004725', 'get_uniprot_for_upi method: correct synonym');

$testdb->hide($dbtype, 'dependent_xref', 'object_xref', 'xref');

# add and store xref
my $ignore_release = 1;
my $xref = $obj->add_xref(\%uniprot, $analysis, $external_dbs);
$dbea->store($xref, $translation->dbID(), 'Translation', $ignore_release);

# overwriting descriptions
is_rows(1, $dba,
  'gene',
  'WHERE stable_id = ? AND description = ?',
  ['AGAP004725', 'eukaryotic translation initiation factor 3 subunit C [Source:VB Community Annotation;Acc:AGAP004725]']
);

$obj->param('description_source', ['reviewed']);
$obj->set_descriptions($dba, $analysis, $external_dbs);
is_rows(1, $dba,
  'gene',
  'WHERE stable_id = ? AND description = ?',
  ['AGAP004725', 'eukaryotic translation initiation factor 3 subunit C [Source:VB Community Annotation;Acc:AGAP004725]']
);
is_rows(0, $dba,
  'gene',
  'WHERE stable_id = ? AND description = ?',
  ['AGAP004725', 'Eukaryotic translation initiation factor 3 subunit C [Source:UniProtKB/Swiss-Prot;Acc:Q7PMU8]']
);

$obj->param('overwrite_description', 1);
$obj->set_descriptions($dba, $analysis, $external_dbs);
is_rows(0, $dba,
  'gene',
  'WHERE stable_id = ? AND description = ?',
  ['AGAP004725', 'eukaryotic translation initiation factor 3 subunit C [Source:VB Community Annotation;Acc:AGAP004725]']
);
is_rows(1, $dba,
  'gene',
  'WHERE stable_id = ? AND description = ?',
  ['AGAP004725', 'Eukaryotic translation initiation factor 3 subunit C [Source:UniProtKB/Swiss-Prot;Acc:Q7PMU8]']
);

$obj->param('description_source', []);
$obj->param('overwrite_description', 0);

# remove_descriptions method
is_rows(114, $dba, 'gene', 'WHERE description IS NOT NULL');
is_rows(1, $dba, 'gene', 'WHERE description RLIKE "Source:UniProtKB/Swiss-Prot"');
is_rows(0, $dba, 'gene', 'WHERE description RLIKE "Source:UniProtKB/TrEMBL"');
foreach my $external_db (values %$external_dbs) {
  $obj->remove_descriptions($dba, $external_db);
}
is_rows(113, $dba, 'gene', 'WHERE description IS NOT NULL');
is_rows(0, $dba, 'gene', 'WHERE description RLIKE "Source:UniProtKB/Swiss-Prot"');
is_rows(0, $dba, 'gene', 'WHERE description RLIKE "Source:UniProtKB/TrEMBL"');

# set_descriptions method
$obj->set_descriptions($dba, $analysis, $external_dbs);
is_rows(0, $dba, 'gene', 'WHERE description RLIKE "Source:UniProtKB/Swiss-Prot"');

$obj->param('description_source', ['reviewed']);
$obj->set_descriptions($dba, $analysis, $external_dbs);
is_rows(1, $dba, 'gene', 'WHERE description RLIKE "Source:UniProtKB/Swiss-Prot"');

# Restore original description, in case the change affects later tests
my $sql = 'UPDATE gene SET description = "eukaryotic translation initiation factor 3 subunit C [Source:VB Community Annotation;Acc:AGAP004725]" WHERE stable_id="AGAP004725";';
my $sth = $dba->dbc->db_handle->prepare($sql);
$sth->execute();

$obj->param('description_source', []);

# overwriting gene names
# Add a display_xref, in order to test the overwriting functionality
$sql = 'UPDATE gene, xref SET display_xref_id=xref_id WHERE stable_id="AGAP004725" AND dbprimary_acc="Q7PMU8";';
$sth = $dba->dbc->db_handle->prepare($sql);
$sth->execute();

is_rows(1, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id',
  'WHERE stable_id = ? AND display_label = ?',
  ['AGAP004725', 'EIF3C_ANOGA']
);

$obj->param('gene_name_source', ['reviewed']);
$obj->set_gene_names($dba, $analysis, $uniprots);
is_rows(1, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id',
  'WHERE stable_id = ? AND display_label = ?',
  ['AGAP004725', 'EIF3C_ANOGA']
);
is_rows(0, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id',
  'WHERE stable_id = ? AND display_label = ?',
  ['AGAP004725', 'eIF3-S8']
);

$obj->param('overwrite_gene_name', 1);
$obj->set_gene_names($dba, $analysis, $uniprots);
is_rows(0, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id',
  'WHERE stable_id = ? AND display_label = ?',
  ['AGAP004725', 'EIF3C_ANOGA']
);
is_rows(1, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id',
  'WHERE stable_id = ? AND display_label = ?',
  ['AGAP004725', 'eIF3-S8']
);

$obj->param('gene_name_source', []);
$obj->param('overwrite_gene_name', 0);

# remove_gene_names method
is_rows(53, $dba, 'gene', 'WHERE display_xref_id IS NOT NULL');
is_rows(1, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id INNER JOIN external_db USING (external_db_id)',
  'WHERE db_name = ?',
  [$obj->param('uniprot_gn_external_db')]
);
$obj->remove_gene_names($dba, $obj->param('uniprot_gn_external_db'));

is_rows(52, $dba, 'gene', 'WHERE display_xref_id IS NOT NULL');
is_rows(0, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id INNER JOIN external_db USING (external_db_id)',
  'WHERE db_name = ?',
  [$obj->param('uniprot_gn_external_db')]
);

# set_gene_names method
$obj->set_gene_names($dba, $analysis, $uniprots);
is_rows(0, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id INNER JOIN external_db USING (external_db_id)',
  'WHERE db_name = ?',
  [$obj->param('uniprot_gn_external_db')]
);

$obj->param('gene_name_source', ['reviewed']);
$obj->set_gene_names($dba, $analysis, $uniprots);
is_rows(1, $dba,
  'gene INNER JOIN xref ON display_xref_id = xref_id INNER JOIN external_db USING (external_db_id)',
  'WHERE db_name = ?',
  [$obj->param('uniprot_gn_external_db')]
);

$obj->param('gene_name_source', []);

# Restore original display_xref, in case the change affects later tests
$sql = 'UPDATE gene SET display_xref_id = 5052195 WHERE stable_id="AGAP004725";';
$sth = $dba->dbc->db_handle->prepare($sql);
$sth->execute();

$testdb->restore($dbtype, 'dependent_xref', 'object_xref', 'xref');
}

# run method
{
my ($reviewed_xrefs, $unreviewed_xrefs) = object_xrefs();
is(scalar(@$reviewed_xrefs), 16,    'run method: correct number of reviewed UniProt xrefs');
is(scalar(@$unreviewed_xrefs), 140, 'run method: correct number of unreviewed UniProt xrefs');
is_rows(0, $dba,
  'dependent_xref INNER JOIN object_xref USING (object_xref_id)',
  'WHERE analysis_id = ?',
  [$analysis->dbID()]
);

# All default parameters
$obj->run();

($reviewed_xrefs, $unreviewed_xrefs) = object_xrefs();
cmp_ok(scalar(@$reviewed_xrefs), '>=', 31,    'run method: correct number of reviewed UniProt xrefs');
cmp_ok(scalar(@$unreviewed_xrefs), '>=', 269, 'run method: correct number of unreviewed UniProt xrefs');
is_rows(145, $dba,
  'dependent_xref INNER JOIN object_xref USING (object_xref_id)',
  'WHERE analysis_id = ?',
  [$analysis->dbID()]
);

# Delete all object_xrefs with matching external_db
$obj->param('replace_all', 1);
$obj->run();

($reviewed_xrefs, $unreviewed_xrefs) = object_xrefs();
cmp_ok(scalar(@$reviewed_xrefs), '>=', 15,    'run method: correct number of reviewed UniProt xrefs');
cmp_ok(scalar(@$unreviewed_xrefs), '>=', 129, 'run method: correct number of unreviewed UniProt xrefs');
is_rows(145, $dba,
  'dependent_xref INNER JOIN object_xref USING (object_xref_id)',
  'WHERE analysis_id = ?',
  [$analysis->dbID()]
);

# Add descriptions
my $descriptions = descriptions();
is($$descriptions{'null'}, 52,                      'run method: correct number of description sources');
ok(!exists($$descriptions{'no source'}),            'run method: correct number of description sources');
ok(!exists($$descriptions{'UniProtKB/Swiss-Prot'}), 'run method: correct number of description sources');
is($$descriptions{'VB Community Annotation'}, 26,   'run method: correct number of description sources');
is($$descriptions{'VB External Description'}, 63,   'run method: correct number of description sources');

$obj->param('description_source', ['reviewed']);
$obj->run();

$descriptions = descriptions();
cmp_ok($$descriptions{'null'}, '<=', 52,                      'run method: correct number of description sources');
ok(!exists($$descriptions{'no source'}),            'run method: correct number of description sources');
ok(!exists($$descriptions{'UniProtKB/Swiss-Prot'}), 'run method: correct number of description sources');
is($$descriptions{'VB Community Annotation'}, 26,   'run method: correct number of description sources');
is($$descriptions{'VB External Description'}, 63,   'run method: correct number of description sources');

# Add descriptions with overwriting
$obj->param('overwrite_description', 1);
$obj->run();

$descriptions = descriptions();
cmp_ok($$descriptions{'null'}, '<=', 52,                    'run method: correct number of description sources');
ok(!exists($$descriptions{'no source'}),                    'run method: correct number of description sources');
cmp_ok($$descriptions{'UniProtKB/Swiss-Prot'}, '>=', 15,    'run method: correct number of description sources');
cmp_ok($$descriptions{'VB Community Annotation'}, '<=', 13, 'run method: correct number of description sources');
cmp_ok($$descriptions{'VB External Description'}, '<=', 61, 'run method: correct number of description sources');

$obj->param('description_source', ['reviewed', 'unreviewed']);
$obj->run();

$descriptions = descriptions();
cmp_ok($$descriptions{'null'}, '<=', 51,                    'run method: correct number of description sources');
ok(!exists($$descriptions{'no source'}),                    'run method: correct number of description sources');
cmp_ok($$descriptions{'UniProtKB/Swiss-Prot'}, '>=', 15,    'run method: correct number of description sources');
cmp_ok($$descriptions{'UniProtKB/TrEMBL'}, '>=', 19,        'run method: correct number of description sources');
cmp_ok($$descriptions{'VB Community Annotation'}, '<=', 7,  'run method: correct number of description sources');
cmp_ok($$descriptions{'VB External Description'}, '<=', 49, 'run method: correct number of description sources');

# Add gene names
my $gene_names = gene_names();
is($$gene_names{'null'}, 113,                   'run method: correct number of gene names');
is($$gene_names{'RFAM'}, 2,                     'run method: correct number of gene names');
is($$gene_names{'TRNASCAN_SE'}, 23,             'run method: correct number of gene names');
is($$gene_names{'VB_Community_Annotation'}, 28, 'run method: correct number of gene names');

$obj->param('gene_name_source', ['reviewed']);
$obj->run();

$gene_names = gene_names();
is($$gene_names{'null'}, 112,                   'run method: correct number of gene names');
is($$gene_names{'RFAM'}, 2,                     'run method: correct number of gene names');
is($$gene_names{'TRNASCAN_SE'}, 23,             'run method: correct number of gene names');
cmp_ok($$gene_names{'Uniprot_gn'}, '>=', 1,     'run method: correct number of gene names');
is($$gene_names{'VB_Community_Annotation'}, 28, 'run method: correct number of gene names');

# Add genes names with overwriting
$obj->param('overwrite_gene_name', 1);
$obj->run();

$gene_names = gene_names();
is($$gene_names{'null'}, 112,                   'run method: correct number of gene names');
is($$gene_names{'RFAM'}, 2,                     'run method: correct number of gene names');
is($$gene_names{'TRNASCAN_SE'}, 23,             'run method: correct number of gene names');
cmp_ok($$gene_names{'Uniprot_gn'}, '>=', 14,    'run method: correct number of gene names');
is($$gene_names{'VB_Community_Annotation'}, 15, 'run method: correct number of gene names');
}

done_testing();

sub object_xrefs {
  my @reviewed;
  my @unreviewed;
  my $translations = $ta->fetch_all();

  foreach my $translation (@$translations) {
    push @reviewed,  @{$translation->get_all_DBEntries('Uniprot/SWISSPROT')};
    push @unreviewed,  @{$translation->get_all_DBEntries('Uniprot/SPTREMBL')};
  }

  return (\@reviewed, \@unreviewed);
}

sub descriptions {
  my %descriptions;
  my $genes = $ga->fetch_all();

  foreach my $gene (@$genes) {
    if (defined $gene->description) {
      if ($gene->description =~ /Source:([^;]+)/) {
        $descriptions{$1}++;
        say $gene->stable_id.': '.$gene->description;
      } else {
        $descriptions{'no source'}++;
      }
    } else {
      $descriptions{'null'}++;
    }
  }

  return (\%descriptions);
}

sub gene_names {
  my %gene_names;
  my $genes = $ga->fetch_all();

  foreach my $gene (@$genes) {
    if (defined $gene->external_db) {
      $gene_names{$gene->external_db}++;
    } else {
      $gene_names{'null'}++;
    }
  }

  return (\%gene_names);
}
