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

Bio::EnsEMBL::EGPipeline::Exonerate::EmailExonerateReport

=head1 DESCRIPTION

Run a few useful queries on the results of exonerate, for a given species.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Exonerate::EmailExonerateReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use Bio::SeqIO;

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'otherfeatures',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $species = $self->param_required('species');
  my $seq_file = $self->param('seq_file');
  my $seq_files = $self->param('seq_files');
  
  my $fasta_file;
  if (exists $$seq_files{$species}) {
    $fasta_file = $$seq_files{$species};
  } else {
    if (defined $seq_file) {
      $fasta_file = $seq_file;
    } else {
      $self->throw("No seq_file for $species.");
    }
  }
  
  my $seqs = Bio::SeqIO->new(-format => 'Fasta', -file => $fasta_file);
  my ($seq_count, $seq_length) = (0, 0);
  while (my $seq = $seqs->next_seq) {
    $seq_count++;
    $seq_length += $seq->length;
  }
  my $mb_length = sprintf("%.0f", $seq_length/1000000);
  my $mean_length = sprintf("%.1f", $seq_length/($seq_count*1000));
    
  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  my $dbh = $dba->dbc->db_handle();
  
  my $sql = 'select count(distinct hit_name) from dna_align_feature;';
  my ($unique_hits) = $dbh->selectrow_array($sql);
  
  my $seq_hit_pcage = sprintf("%.0f", ($unique_hits/$seq_count)*100);
  
  my $text = 
    "The exonerate pipeline has completed for $species, ".
    "using the sequence file $fasta_file.\n".
    "That file has $seq_count sequences and a total length of $mb_length Mb, ".
    "giving an average sequence length of $mean_length Kb. ".
    "Of that total, $unique_hits ($seq_hit_pcage%) were mapped to the genome.\n\n".
    "Fond regards,\nThe Exonerate Pipeline\n";
  
  $self->param('text', $text);
}

1;
