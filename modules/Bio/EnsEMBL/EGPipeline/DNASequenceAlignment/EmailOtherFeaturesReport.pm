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

Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::EmailOtherFeaturesReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::EmailOtherFeaturesReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use Bio::SeqIO;

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type'   => 'otherfeatures',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $species      = $self->param_required('species');
  my $file         = $self->param('file');
  my $file_species = $self->param('file_species');
  my $aligner      = $self->param_required('aligner');
  my $data_type    = $self->param_required('data_type');
  my $logic_name   = $self->param_required('logic_name');
  
  my $dba = $self->get_DBAdaptor($self->param_required('db_type'));
  my $dbh = $dba->dbc->db_handle();
  
  my $sql =
    'SELECT COUNT(distinct hit_name) FROM '.
    'dna_align_feature INNER JOIN analysis USING (analysis_id) '.
    'WHERE logic_name = "'.$logic_name.'";';
  my ($unique_hits) = $dbh->selectrow_array($sql);
  
  my $text = 
    "The DNA Sequence Alignment pipeline has completed for $species, ".
    "using $data_type sequences, and the $aligner aligner. ";
  
  if (exists $$file_species{$species}) {
    push @$file, $$file_species{$species};
  }
  
  if (scalar(@$file)) {
    my $files = join(', ', @$file);
  
    my ($seq_count, $seq_length) = (0, 0);
    foreach my $fasta_file (@$file) {
      my $seqs = Bio::SeqIO->new(-format => 'Fasta', -file => $fasta_file);
      while (my $seq = $seqs->next_seq) {
        $seq_count++;
        $seq_length += $seq->length;
      }
    }
    my $mb_length = sprintf("%.0f", $seq_length/1000000);    
    my $seq_hit_pcage = sprintf("%.0f", ($unique_hits/$seq_count)*100);
    
    $text .=
      "Across all DNA files ($files) there are $seq_count sequences ".
      "with a total length of $mb_length Mb. ".
      "Of that total, $unique_hits unique hits ($seq_hit_pcage%) ".
      "were mapped to the genome.\n\n";
  } else {
    $text .= "$unique_hits unique hits were mapped to the genome.\n\n";
  }
  
  $self->param('text', $text);
}

1;
