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

Bio::EnsEMBL::EGPipeline::LoadGFF3::CheckProteinSeqs

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::LoadGFF3::CheckProteinSeqs;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::LoadGFF3::Base');

sub param_defaults {
  my ($self) = @_;
  return {
    db_type => 'core',
  };
}

sub run {
  my ($self) = @_;
  my $db_type            = $self->param_required('db_type');
  my $logic_name         = $self->param('logic_name');
  my $protein_fasta_file = $self->param('protein_fasta_file');
  
  my %protein = $self->load_fasta($protein_fasta_file);
  
  my $dba = $self->get_DBAdaptor($db_type);
  my $ta  = $dba->get_adaptor('Transcript');
  
  my $transcripts = $ta->fetch_all_by_logic_name($logic_name);
  
  foreach my $transcript (@$transcripts) {
    my $translation = $transcript->translation;
    
    if ($translation && exists $protein{$translation->stable_id}) {
      my $file_seq = $protein{$translation->stable_id};
      say $transcript->stable_id;
      my $shift_success = $self->shift_translation_start($dba, $transcript, $file_seq);
      
      if (! $shift_success) {
        my @file_seq = split(//, $file_seq);
        my @atts;
        
        my $db_seq = $translation->seq;
        
        if (length($file_seq) == length($db_seq)) {
          while ((my $pos = index($db_seq, '*')) >= 0) {
            my $amino_acid = $file_seq[$pos];
            $db_seq =~ s/\*/$amino_acid/;
            
            $pos += 1;      
            my $att = $self->add_translation_seq_edit($translation, $pos, $pos, $amino_acid);
            push @atts, $att;
          }
          
          if ($db_seq eq $file_seq) {
            $aa->store_on_Translation($translation, \@atts);
          }
        } else {
          $self->warning('Protein sequence length mismatch for '.$translation->stable_id);
        }
      }
    }
  }
}

sub fetch_input {
  my ($self) = @_;
  my $species    = $self->param_required('species');
  my $db_type    = $self->param_required('db_type');
  my $logic_name = $self->param_required('logic_name');
  
  my $dbh = $self->get_DBAdaptor($db_type)->dbc->db_handle;
  
  my $report = $self->biotype_report($dbh, $logic_name);
  $report   .= $self->seq_edit_tt_report($dbh, $logic_name);
  $report   .= $self->seq_edit_tn_report($dbh, $logic_name);
  
  $self->param('text', $report);
}

sub biotype_report {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      g.biotype AS gene_biotype, 
      t.biotype AS transcript_biotype, 
      COUNT(DISTINCT g.stable_id) AS count_of_genes, 
      COUNT(DISTINCT t.stable_id) AS count_of_transcripts 
    FROM 
      gene g INNER JOIN 
      transcript t USING (gene_id) INNER JOIN 
      analysis a ON g.analysis_id = a.analysis_id 
    WHERE 
      a.logic_name = ? 
    GROUP BY 
      g.biotype, 
      t.biotype 
    ORDER BY 
      g.biotype, 
      t.biotype 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = "Imported genes summarised by biotype:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub seq_edit_tt_report {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      t.biotype, 
      COUNT(DISTINCT t.stable_id) AS count_of_transcripts, 
      COUNT(*) AS count_of_seq_edits 
    FROM 
      analysis a INNER JOIN 
      transcript t USING (analysis_id) INNER JOIN
      transcript_attrib ta USING (transcript_id) INNER JOIN
      attrib_type at USING (attrib_type_id)
    WHERE 
      a.logic_name = ? AND
      at.code = '_rna_edit'
    GROUP BY 
      t.biotype 
    ORDER BY 
      t.biotype 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = "Imported genes requiring transcript-level sequence edits:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub seq_edit_tn_report {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      tt.biotype,
      at.name,
      COUNT(DISTINCT tn.stable_id) AS count_of_translations, 
      COUNT(*) AS count_of_seq_edits 
    FROM 
      analysis a INNER JOIN 
      transcript tt USING (analysis_id) INNER JOIN
      translation tn USING (transcript_id) INNER JOIN
      translation_attrib ta USING (translation_id) INNER JOIN
      attrib_type at USING (attrib_type_id)
    WHERE 
      a.logic_name = ? AND
      at.code IN ('_selenocysteine', 'amino_acid_sub')
    GROUP BY 
      tt.biotype,
      at.name 
    ORDER BY 
      tt.biotype,
      at.name 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = "Imported genes requiring translation-level sequence edits:";
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

1;
