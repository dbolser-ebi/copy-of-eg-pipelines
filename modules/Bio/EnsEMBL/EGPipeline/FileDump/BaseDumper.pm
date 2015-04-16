=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);

sub param_defaults {
  return {
    'db_type'            => 'core',
    'eg_dir_structure'   => 0,
    'eg_filename_format' => 0,
  };
}

sub fetch_input {
  my ($self) = @_;
  my $out_file = $self->param('out_file');
  
  if (!defined $out_file) {
    $out_file = $self->generate_filename();
    $self->param('out_file', $out_file);
    $self->param('out_files', [$out_file]);
  }
  open(my $fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
  $self->param('out_fh', $fh);
}

sub generate_filename {
  my ($self) = @_;
  
  my $species            = $self->param('species');
  my $file_type          = $self->param('file_type');
  my $pipeline_dir       = $self->param('pipeline_dir');
  my $eg_dir_structure   = $self->param('eg_dir_structure');
  my $eg_filename_format = $self->param('eg_filename_format');
  
  if ($eg_dir_structure) {
    my $division = $self->get_division();
    $pipeline_dir = catdir($pipeline_dir, $division, $file_type, $species);
    $self->param('pipeline_dir', $pipeline_dir);
  }
  make_path($pipeline_dir);
  
  my $filename;
  if ($eg_filename_format) {
    $filename = $self->generate_eg_filename();
  } else {
    $filename = $self->generate_vb_filename();
  }
  
  return catdir($pipeline_dir, $filename);
}

sub generate_eg_filename {
  my ($self) = @_;
  
  my $species  = $self->param('species');
  my $file_type = $self->param('file_type');
  
  my $dba = $self->core_dba;
  my $dbname = $dba->dbc->dbname();
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  my ($eg_version) = $dbname =~ /([^_]+)_[^_]+_[^_]+$/;
  my $filename = ucfirst($species).".$assembly.$eg_version.$file_type";
  
  return $filename;
}

sub generate_vb_filename {
  my ($self) = @_;
  
  my $species  = $self->param('species');
  my $data_type = $self->param('data_type');
  my $file_type = $self->param('file_type');
  
  $species =~ s/_/-/;
  my $dba = $self->core_dba;
  my $strain = $dba->get_MetaContainer()->single_value_by_key('species.strain');
  my $gene_set = $dba->get_MetaContainer()->single_value_by_key('genebuild.version');
  my $filename = ucfirst($species).'-'.$strain.'_'.uc($data_type).'_'."$gene_set.$file_type";
  
  return $filename;
}

sub get_division {
  my ($self) = @_;
  
  my $dba = $self->core_dba;  
  my $division;
  if ($dba->dbc->dbname() =~ /(\w+)\_\d+_collection_/) {
    $division = $1;
  } else {
    $division = $dba->get_MetaContainer->get_division();
    $division = lc($division);
    $division =~ s/ensembl//;
  }
  return $division;
}

sub has_chromosome {
  my ($self, $dba) = @_;
  my $helper = $dba->dbc->sql_helper();
  my $sql = q{
    SELECT COUNT(*) FROM
    coord_system cs INNER JOIN
    seq_region sr USING (coord_system_id) INNER JOIN
    seq_region_attrib sa USING (seq_region_id) INNER JOIN
    attrib_type at USING (attrib_type_id)
    WHERE cs.species_id = ?
    AND at.code = 'karyotype_rank'
  };
  my $count = $helper->execute_single_result(-SQL => $sql, -PARAMS => [$dba->species_id()]);
  
  $dba->dbc->disconnect_if_idle();
  
  return $count;
}

sub write_output {
  my ($self) = @_;
  
  foreach my $out_file (@{$self->param('out_files')}) {
    $self->dataflow_output_id({out_file => $out_file}, 1);
  }
}

1;
