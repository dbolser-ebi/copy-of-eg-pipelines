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

Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAGenesReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAGenesReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

sub fetch_input {
  my ($self) = @_;
  my $species    = $self->param_required('species');
  my $logic_name = $self->param_required('logic_name');
  
  my $dbh = $self->core_dbh();
  
  my $reports = $self->text_summary($dbh, $logic_name);
  $reports   .= $self->report_one($dbh, $logic_name);
  
  $self->param('text', $reports);
}

sub text_summary {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      COUNT(*) AS count_of_genes
    FROM  
      gene INNER JOIN 
      analysis USING (analysis_id)
    WHERE logic_name = ? 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $results = $sth->fetchall_arrayref();
  my $count   = $$results[0][0];
  
  my $summary =
    "The RNA genes pipeline has completed, having generated $count genes.\n".
    "The genes are linked to the '$logic_name' analysis.\n";
  
  return $summary;
}

sub report_one {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      biotype, 
      COUNT(*) AS count_of_genes
    FROM  
      gene INNER JOIN 
      analysis USING (analysis_id)
    WHERE logic_name = ? 
    GROUP BY 
      biotype 
    ORDER BY 
      biotype 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = "Biotype summary:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

1;
