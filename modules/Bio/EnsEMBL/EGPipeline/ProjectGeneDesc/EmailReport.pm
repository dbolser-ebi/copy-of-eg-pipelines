=head1 LICENSE

Copyright [1999-2015] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::EmailReport

=head1 DESCRIPTION

Email a summary of gene descriptions that were projected.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProjectGeneDesc::EmailReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

sub fetch_input {
  my ($self) = @_;
  
  my $log_file     = $self->param_required('log_file');
  my $store_data   = $self->param_required('store_data');
  my $from_species = $self->param_required('source');
  my $to_species   = $self->param_required('species');
  
  my $reports;
  if (!$store_data) {
    $reports .= "Pipeline was run with a flag to prevent storing anything in the database.\n";
  } else {
    $reports .= $self->description_summary($from_species, $to_species);
  }
  $reports .= "Detailed log file: $log_file.\n\n";
  
  $self->param('text', $reports);
}

sub description_summary {
  my ($self, $from_species, $to_species) = @_;
  
  my $from_species_text = ucfirst($from_species);
  $from_species_text =~ s/_/ /g;
  my $desc_match = "Projected from $from_species_text";
  
  my $dbh = $self->core_dbh();
  my $sql = 'SELECT COUNT(*) AS Total FROM gene WHERE description LIKE "%?%"';
  my $sth = $dbh->prepare($sql);
  $sth->execute($desc_match);
  
  my $title = "Number of projected descriptions from $from_species to $to_species:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);    
}

1;
