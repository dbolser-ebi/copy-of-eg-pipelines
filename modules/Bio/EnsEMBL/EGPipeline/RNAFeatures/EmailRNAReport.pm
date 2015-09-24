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

Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $species      = $self->param_required('species');
  my $run_cmscan   = $self->param_required('run_cmscan');
  my $run_trnascan = $self->param_required('run_trnascan');
  
  my $dba = $self->get_DBAdaptor($self->param_required('db_type'));
  my $dbh = $dba->dbc->db_handle();
  
  my $text;
  
  if ($run_cmscan) {
    my $rfam_logic_name   = $self->param_required('rfam_logic_name');
    my $cmscan_cm_file    = $self->param_required('cmscan_cm_file');
    my $cmscan_logic_name = $self->param_required('cmscan_logic_name');
    
    my $logic_name;
    if (! exists $$cmscan_cm_file{$species} && ! exists $$cmscan_cm_file{'all'}) {
      $logic_name = $rfam_logic_name;
    } elsif (exists $$cmscan_logic_name{$species}) {
      $logic_name = $$cmscan_logic_name{$species};
    } elsif (exists $$cmscan_logic_name{'all'}) {
      $logic_name = $$cmscan_logic_name{'all'};
    } else {
      $logic_name = 'cmscan_custom';
    }
    
    my ($unique, $low_conf, $high_conf) = $self->report_counts($dbh, $logic_name);
    my $summary = $self->report_summary($dbh, $logic_name);
    
    my $text = 
      "The RNA Features pipeline has completed for $species, having aligned ".
      "$unique distinct covariance models with cmscan.\n\n".
      "There are $low_conf low confidence (1e-3 >= e-value > 1e-6) alignments, ".
      "and $high_conf high confidence (e-value <= 1e-6) alignments.\n\n$summary";
  }
  
  $self->param('text', $text);
}

sub report_counts {
  my ($self, $dbh, $logic_name) = @_;
  
  my $unique_sql =
    'SELECT COUNT(distinct hit_name) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'";';
  my ($unique) = $dbh->selectrow_array($unique_sql);
  
  my $low_conf_sql =
    'SELECT COUNT(*) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'" AND evalue <= 1e-3 AND evalue > 1e-6;';
  my ($low_conf) = $dbh->selectrow_array($low_conf_sql);
  
  my $high_conf_sql =
    'SELECT COUNT(*) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'" AND evalue <= 1e-6;';
  my ($high_conf) = $dbh->selectrow_array($high_conf_sql);
  
  return ($unique, $low_conf, $high_conf);
}

sub report_summary {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT
      LEFT(
        SUBSTRING(
          external_data,
          INSTR(
            external_data,
            'Biotype='
          )+8
        ),
        INSTR(
          SUBSTRING(
            external_data,
            INSTR(
              external_data,
              'Biotype='
            )+8
          ),
          ';'
        )-1
      ) AS biotype,
      COUNT(*) AS count_of_alignments
    FROM dna_align_feature
    GROUP BY biotype
    ORDER BY biotype
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  
  my $title = "Summary of RNA gene biotypes:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

1;
