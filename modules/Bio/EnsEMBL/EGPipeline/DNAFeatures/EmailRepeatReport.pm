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

Bio::EnsEMBL::EGPipeline::DNAFeatures::EmailRepeatReport

=head1 DESCRIPTION

Run a few useful queries on the repeat features, for a given species.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::DNAFeatures::EmailRepeatReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

sub fetch_input {
  my ($self) = @_;
  my $species = $self->param_required('species');
  
  my $dbh = $self->core_dbh();
  
  my $seq_region_length_sql = '
    SELECT sum(length) FROM
      seq_region INNER JOIN
      seq_region_attrib USING (seq_region_id) INNER JOIN
      attrib_type USING (attrib_type_id)
    WHERE code = "toplevel"
  ;';
  my ($seq_region_length) = $dbh->selectrow_array($seq_region_length_sql);
  
  my $reports = $self->text_summary($dbh, $seq_region_length);
  $reports .= $self->report_one($dbh, $seq_region_length);
  $reports .= $self->report_two($dbh, $seq_region_length);
  
  $self->param('text', $reports);
}

sub text_summary {
  my ($self, $dbh, $seq_region_length) = @_;
  
  my $sql = "
    SELECT 
      display_label, 
      logic_name, 
      COUNT(*) AS count_of_features, 
      SUM(seq_region_end - seq_region_start+1) AS total_repeat_length, 
      ROUND(SUM(seq_region_end - seq_region_start+1) / $seq_region_length, 4) AS repeat_coverage 
    FROM  
      analysis LEFT OUTER JOIN 
      analysis_description USING (analysis_id) INNER JOIN 
      repeat_feature USING (analysis_id) 
    GROUP BY 
      display_label, 
      logic_name 
    ORDER BY 
      display_label 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  
  my @results_text;
  my $results = $sth->fetchall_arrayref();
  foreach my $result (@$results) {
    my ($display_label, $logic_name, $count, $length, $coverage) = @$result;
    $length = sprintf('%.0f', $length/1000000);
    $coverage = sprintf('%.1f', $coverage*100);
    if ($logic_name =~ /^repeatmask/) {
      my $library_text;
      if ($logic_name =~ /customlib/) {
        $library_text = "a custom";
      } else {
        $library_text = "the RepBase";
      }
      push @results_text,
        "$count RepeatMasker features (with $library_text library), ".
        "covering $length Mb ($coverage% of the genome)";
    } else {
      push @results_text,
        "$count $display_label features, ".
        "covering $length Mb ($coverage% of the genome)";
    }
  }
  
  my $text = "Recommended summary text for static content:\n".
    "Repeats were annotated with the <a href='http://ensemblgenomes.org/info/data/repeat_features'>".
    "Ensembl Genomes repeat feature pipeline</a>. There are: ";
  $text .= join("; ", @results_text);
  $text .= ".\n\n";
  
  return $text;
}

sub report_one {
  my ($self, $dbh, $seq_region_length) = @_;
  
  my $sql = "
    SELECT 
      display_label, 
      program, 
      created, 
      logic_name, 
      program_file, 
      parameters, 
      COUNT(*) AS count_of_features, 
      SUM(seq_region_end - seq_region_start+1) AS total_repeat_length, 
      ROUND(SUM(seq_region_end - seq_region_start+1) / $seq_region_length, 4) AS repeat_coverage 
    FROM  
      analysis LEFT OUTER JOIN 
      analysis_description USING (analysis_id) INNER JOIN 
      repeat_feature USING (analysis_id) 
    GROUP BY 
      display_label, 
      program, 
      created, 
      logic_name, 
      program_file, 
      parameters 
    ORDER BY 
      display_label 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  
  my $title = "Simple summary of repeat features:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub report_two {
  my ($self, $dbh, $seq_region_length) = @_;
  
  my $sql = "
    SELECT 
      display_label, 
      program, 
      created, 
      repeat_type,
      repeat_class,
      COUNT(*) AS count_of_features, 
      SUM(seq_region_end - seq_region_start+1) AS total_repeat_length, 
      ROUND(SUM(seq_region_end - seq_region_start+1) / analysis_length, 4) AS repeat_proportion
    FROM  
      analysis LEFT OUTER JOIN 
      analysis_description USING (analysis_id) INNER JOIN 
      repeat_feature USING (analysis_id) INNER JOIN
      repeat_consensus USING (repeat_consensus_id) INNER JOIN
      ( SELECT 
          analysis_id, 
          SUM(seq_region_end - seq_region_start+1) as analysis_length 
        FROM repeat_feature
        GROUP BY analysis_id
      ) as analysis_length_tbl USING (analysis_id)
    GROUP BY 
      display_label, 
      program, 
      created, 
      repeat_type,
      repeat_class
    ORDER BY 
      display_label 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  
  my $title = "Summary of repeat features, by type and class:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

1;
