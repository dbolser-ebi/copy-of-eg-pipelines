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
  
  my $counts;
  my @logic_names;
  
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
    push @logic_names, $logic_name;
    
    my ($unique, $low_conf, $high_conf) = $self->cmscan_counts($dbh, $logic_name);
    
    $counts .= 
      "The pipeline aligned $unique distinct covariance models with cmscan ".
      "(logic_name: $logic_name).\n\n".
      "There are $low_conf low confidence (1e-3 >= e-value > 1e-6) alignments, ".
      "and $high_conf high confidence (e-value <= 1e-6) alignments.\n\n";
  }
  
  if ($run_trnascan) {
    my $logic_name = 'trnascan';
    push @logic_names, $logic_name;
    
    my ($unique, $low_conf, $high_conf) = $self->trnascan_counts($dbh, $logic_name);
    
    $counts .= 
      "The pipeline aligned $unique distinct tRNA models with tRNAscan-SE ".
      "(logic_name: $logic_name).\n\n".
      "There are $low_conf low confidence (COVE score < 40) alignments, ".
      "and $high_conf high confidence (COVE score >= 40) alignments.\n\n";
  }
  
  my $summary = $self->report_summary($dbh, \@logic_names);
  my $text =
    "The RNA Features pipeline has completed for $species.\n\n".
    "$counts\n\n".
    "$summary\n\n";
  
  $self->param('text', $text);
}

sub cmscan_counts {
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

sub trnascan_counts {
  my ($self, $dbh, $logic_name) = @_;
  
  my $unique_sql =
    'SELECT COUNT(distinct hit_name) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'";';
  my ($unique) = $dbh->selectrow_array($unique_sql);
  
  my $low_conf_sql =
    'SELECT COUNT(*) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'" AND score < 40;';
  my ($low_conf) = $dbh->selectrow_array($low_conf_sql);
  
  my $high_conf_sql =
    'SELECT COUNT(*) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'" AND score >= 40;';
  my ($high_conf) = $dbh->selectrow_array($high_conf_sql);
  
  return ($unique, $low_conf, $high_conf);
}

sub report_summary {
  my ($self, $dbh, $logic_names) = @_;
  
  my $logic_name_list = "'" . join("','", @$logic_names) . "'";
  
  my $sql = "
    SELECT
      LEFT(
        SUBSTRING(
          external_data,
          INSTR(
            external_data,
            \"'Biotype' => \"
          )+14
        ),
        INSTR(
          SUBSTRING(
            external_data,
            INSTR(
              external_data,
            \"'Biotype' => \"
            )+14
          ),
          \"'\"
        )-1
      ) AS biotype,
      COUNT(*) AS count_of_alignments
    FROM dna_align_feature
    INNER JOIN analysis USING (analysis_id)
    WHERE logic_name in ($logic_name_list)
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
