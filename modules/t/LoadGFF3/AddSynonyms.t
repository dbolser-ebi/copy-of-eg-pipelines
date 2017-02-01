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
use Bio::EnsEMBL::EGPipeline::LoadGFF3::AddSynonyms;

my $species  = 'anopheles_gambiae';
my $curr_dir = $FindBin::Bin.'/../';

my $testdb = Bio::EnsEMBL::Test::MultiTestDB->new($species, $curr_dir);
my $dbtype = 'core';
my $dba    = $testdb->get_DBAdaptor($dbtype);

my $fasta_file = $FindBin::Bin.'/../test-files/LoadGFF3/agam.fa';

my $module_name    = 'Bio::EnsEMBL::EGPipeline::LoadGFF3::AddSynonyms';
my @hive_methods   = qw(param_defaults fetch_input run write_output);
my @module_methods = qw(add_synonym_to_db fetch_slices load_fasta slice_synonyms);
can_ok($module_name, @hive_methods);
can_ok($module_name, @module_methods);

# Create an instance of the module; a dummy job object is required to
# prevent errors when the module generates log messages.
my $as_obj  = Bio::EnsEMBL::EGPipeline::LoadGFF3::AddSynonyms->new;
my $job_obj = Bio::EnsEMBL::Hive::AnalysisJob->new;
$as_obj->input_job($job_obj);

# Set and check default parameters.
my $param_defaults = $as_obj->param_defaults();
$as_obj->input_job->param_init($param_defaults);
is($as_obj->param('db_type'),             'core',                     'param_defaults method: db_type');
is($as_obj->param('synonym_external_db'), 'ensembl_internal_synonym', 'param_defaults method: synonym_external_db');

# fetch_slices method
my %slices = $as_obj->fetch_slices($dba);
is(scalar(keys %slices), 3, 'fetch_slice method: hash with three keys returned');
ok(exists($slices{'2L'}), 'fetch_slice method: slice for chr 2L');
ok(exists($slices{'CM000356.1'}), 'fetch_slice method: slice for chr 2L synonym (CM000356.1)');
ok(exists($slices{'Mt'}), 'fetch_slice method: slice for chr Mt');
isa_ok($slices{'2L'}, 'Bio::EnsEMBL::Slice', 'fetch_slice method: hash value');

# load_fasta method
my %fasta = $as_obj->load_fasta($fasta_file);
is(scalar(keys %fasta), 2, 'load_fasta method: hash with two keys returned');
ok(exists($fasta{'2L'}), 'load_fasta method: sequence for chr 2L');
ok(exists($fasta{'L20934'}), 'load_fasta method: sequence for chr Mt synonym (L20934)');
is(length($fasta{'2L'}), 4000000, 'load_fasta method: 4Mb of 2L sequence');
is(length($fasta{'L20934'}), 15363, 'load_fasta method: 15363b of L20934 sequence');
ok($fasta{'2L'} !~ /\s/m, 'load_fasta method: no whitespace in sequence');

# Synonym found and added successfully.
$testdb->hide($dbtype, 'seq_region_synonym');
$as_obj->slice_synonyms($dba, \%slices, \%fasta);
is_rows(1, $dba, 'seq_region_synonym', 'where synonym = ? ', ['L20934']);
$testdb->restore($dbtype, 'seq_region_synonym');

# Refetch slices, because the Mt one will have been updated with a synonym.
%slices = $as_obj->fetch_slices($dba);

# Synonym cannot be added if external_db is not in the database.
$as_obj->param('synonym_external_db', 'mystery');
dies_ok { $as_obj->slice_synonyms($dba, \%slices, \%fasta) } 'Synonym requires extant external_db';
$as_obj->param('synonym_external_db', 'ensembl_internal_synonym');

# Synonym not added unless sequence matches exactly.
$fasta{'L20934'} .= 'A';
$as_obj->slice_synonyms($dba, \%slices, \%fasta);
is_rows(0, $dba, 'seq_region_synonym', 'where synonym = ? ', ['L20934']);
$fasta{'L20934'} =~ s/A$//;

done_testing();
