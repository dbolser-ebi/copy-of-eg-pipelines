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

Bio::EnsEMBL::EGPipeline::LoadGFF3::Base

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::LoadGFF3::Base;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub fetch_slices {
  my ($self, $dba) = @_;
  
  my %slices;
  
  my $slice_adaptor = $dba->get_adaptor("Slice");
  my $synonym_adaptor = $dba->get_adaptor("SeqRegionSynonym");
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    my $seq_region_name = $slice->seq_region_name;
    $slices{$seq_region_name} = $slice;
    my $synonyms = $synonym_adaptor->get_synonyms($slice->get_seq_region_id);
    foreach my $synonym (@$synonyms) {
      $slices{$synonym->name} = $slices{$seq_region_name};
    }
  }
  
  return %slices;
}

sub load_fasta {
  my ($self, $fasta_file) = @_;
  my %fasta;
  
  local $/;
  open(FASTA, $fasta_file);
  %fasta = <FASTA> =~ />(\S+).*\n([^\>]+)/gm;
  close(FASTA);
  
  foreach my $id (keys %fasta) {
    $fasta{$id} =~ s/\s//gm;
    $fasta{$id} = uc($fasta{$id});
  }
  
  return %fasta;
}

sub set_protein_coding {
  my ($self, $dba, $logic_name) = @_;
  
  my $ga = $dba->get_adaptor('Gene');
  my $ta = $dba->get_adaptor('Transcript');
  
  my $genes = $ga->fetch_all_by_logic_name($logic_name);
  
  foreach my $gene (@$genes) {
    my $nontranslating_transcript = 0;
    my $protein_coding_transcript = 0;
    
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      if ($transcript->translation) {
        if ($transcript->translation->seq =~ /\*/) {
          $transcript->biotype('nontranslating_CDS');
          $nontranslating_transcript++;
        } else {
          $transcript->biotype('protein_coding');
          $protein_coding_transcript++;
        }
        $ta->update($transcript);
      }
    }
    
    if ($protein_coding_transcript) {
      $gene->biotype('protein_coding');
      $ga->update($gene);
    } elsif ($nontranslating_transcript) {
      $gene->biotype('nontranslating_CDS');
      $ga->update($gene);
    }
  }
}

sub update_translation_start {
  my ($self, $dba, $dbid, $start) = @_;
  
  my $sql = 'UPDATE translation SET seq_start = ? WHERE translation_id = ?;';
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($start, $dbid) or $self->throw("Failed to execute: $sql");
}

1;
