=head1 LICENSE

Copyright [1999-2014] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::ProteinSimilarity::AnalysisSetup

=head1 DESCRIPTION

Backup tables; if analysis exists, delete any rows connected to
it, and the analysis itself; create the analysis anew. By default
will operate on both core and otherfeature databases.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::AnalysisSetup;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Hive::Process);
use Bio::EnsEMBL::Analysis;

sub param_defaults {
  return {
    'db_type'            => 'core',
    'linked_tables'      => [],
    'delete_existing'    => 0,
    'db_backup_required' => 1,
    'production_lookup'  => 1,
    'production_db_name' => 'ensembl_production',
    'db'                 => undef,
    'db_version'         => undef,
    'db_file'            => undef,
    'program'            => undef,
    'program_version'    => undef,
    'program_file'       => undef,
    'parameters'         => undef,
    'module'             => undef,
    'module_version'     => undef,
    'gff_source'         => undef,
    'gff_feature'        => undef,
    'description'        => undef,
    'display_label'      => undef,
    'displayable'        => 1,
  };
}

sub fetch_input {
  my $self = shift @_;
  
  my $logic_name = $self->param_required('logic_name');
  $self->param('program', $logic_name) unless $self->param_is_defined('program');
  
  # Expect to be passed a db dump
  if ($self->param('db_backup_required')) {
    my $db_backup_file = $self->param_required('db_backup_file');
    
    if (!-e $db_backup_file) {
      $self->throw("Database backup file '$db_backup_file' does not exist");
    }
  }
  
}

sub write_output {
  my $self = shift @_;
  
  my $output_id = {
    'db_type'      => $self->param('db_type'),
    'logic_name'   => $self->param('logic_name'),
    'program_file' => $self->param('program_file'),
    'parameters'   => $self->param('parameters'),
  };
  $self->dataflow_output_id($output_id, 1);
  
}

sub run {
  my $self = shift @_;
  my $species = $self->param('species');
  
  my $aa = Bio::EnsEMBL::Registry->get_adaptor($species, $self->param('db_type'), 'Analysis');
  my $dbh = $aa->dbc->db_handle;
  my $analysis = $aa->fetch_by_logic_name($self->param('logic_name'));
  
  if (defined $analysis) {
    # Dump database here...
    
    
    my $analysis_id = $analysis->dbID;
    foreach my $table (@{$self->param('feature_tables')}) {
      my $sql = "DELETE FROM $table WHERE analysis_id = $analysis_id";
      my $sth = $dbh->prepare($sql) or throw("Failed to delete rows using '$sql': ".$dbh->errstr);
      $sth->execute or throw("Failed to delete rows using '$sql': ".$sth->errstr);
    }
    $aa->remove($analysis);
  }
  
  $aa->store($self->create_analysis);
  
}

sub create_analysis {
  my ($self) = @_;
  
  my $analysis = Bio::EnsEMBL::Analysis->new(
    -logic_name      => $self->param('logic_name'),
    -db              => $self->param('db'),
    -db_version      => $self->param('db_version'),
    -db_file         => $self->param('db_file'),
    -program         => $self->param('program'),
    -program_version => $self->param('program_version'),
    -program_file    => $self->param('program_file'),
    -parameters      => $self->param('parameters'),
    -module          => $self->param('module'),
    -module_version  => $self->param('module_version'),
    -gff_source      => $self->param('gff_source'),
    -gff_feature     => $self->param('gff_feature'),
    -description     => $self->param('description'),
    -display_label   => $self->param('display_label'),
    -displayable     => $self->param('displayable'),
  );
  
  return $analysis;
}

1;
