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

Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::LoadGFF3::ApplySeqEdits;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::LoadGFF3::Base');

use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  return {
    db_type => 'core',
  };
}

sub run {
  my ($self) = @_;
  my $db_type            = $self->param_required('db_type');
  my $logic_name         = $self->param_required('logic_name');
  my $genbank_file       = $self->param('genbank_file');
  my $protein_fasta_file = $self->param('protein_fasta_file');
  
  my $dba = $self->get_DBAdaptor($db_type);
  
  if ($genbank_file) {
    $self->seq_edits_from_genbank($dba, $genbank_file);
  }
  
  if ($protein_fasta_file) {
    $self->seq_edits_from_protein($dba, $logic_name, $protein_fasta_file);
  }
  
  $self->set_protein_coding($dba, $logic_name);
}

sub seq_edits_from_genbank {
  my ($self, $dba, $genbank_file) = @_;
  
  my $ta = $dba->get_adaptor('Transcript');
  my $aa = $dba->get_adaptor("Attribute");
  
  my $genbank_path = path($genbank_file);
  my $genbank = $genbank_path->slurp;
  $genbank =~ s!\A//\s*!!m;
  $genbank =~ s!//\s*\Z!!m;
  my @genbank = split(m!//\s+!, $genbank);
  
  foreach my $record (@genbank) {
    if ($record =~ /##RefSeq-Attributes-START##.*(assembly gap|frameshifts)/ms) {
      
      # Probably will only have mRNA, but check just in case.
      my ($mol_type) = $record =~ /\s+\/mol_type="([^"]+)"/m;
      next if $mol_type ne 'mRNA';
      
      my ($accession) = $record =~ /^VERSION\s+(\S+)/m;
      my $transcript  = $ta->fetch_by_stable_id($accession);
      
      my $cdna_seq = $self->extract_seq($record);
      
      # Check that the sequences don't already match before trying anything...
      if ($transcript->seq->seq ne $cdna_seq) {
        my $edits = $self->extract_edits($record);
        
        my ($five_prime_adjust, $seq_start) = $self->adjust_five_prime_edit($record, $edits);
        if ($five_prime_adjust) {
          $cdna_seq =~ s/^.{$five_prime_adjust}//;
          
          if (defined $seq_start) {
            $self->update_translation_start($dba, $transcript->translation->dbID, $seq_start);
            $transcript = $ta->fetch_by_stable_id($accession);
          }
        }
        
        if (scalar(@$edits)) {
          my $three_prime_adjust = $self->adjust_three_prime_edit($record, $edits);
          if ($three_prime_adjust) {
            $cdna_seq =~ s/.{$three_prime_adjust}$//;
          }
        }
        
        # Need to track the length of the inserts we're making and subtract
        # that from subsequent edit starts, because the seq_edits all need
        # to be relative to the original transcript...
        my $edit_length = 0;
        my @atts;
        
        foreach my $edit (@$edits) {
          my ($seq_name, $t_start, $t_end) = @$edit;
          
          # This looks weird, but having $to one less than $from
          # is how you indicate a sequence insertion...
          my $from = $t_start - $edit_length;
          my $to   = $from - 1;
          
          my $subseq;
          if ($seq_name =~ /"(N+)"/) {
            $subseq = "$1";
          } else {
            $subseq = substr($cdna_seq, $t_start-1, $t_end-$t_start+1);
          }
          
          my $att = $self->add_transcript_seq_edit($transcript, $from, $to, $subseq);
          push @atts, $att;
          
          $edit_length += length($subseq);
        }
        
        if ($transcript->seq->seq eq $cdna_seq) {
          $aa->store_on_Transcript($transcript, \@atts);
        } else {          
          my $att = $self->add_transcript_seq_edit($transcript, 1, $transcript->length, $cdna_seq);
          $aa->store_on_Transcript($transcript, [$att]);
        }
      }
    }
  }
}

sub extract_edits {
  my ($self, $record) = @_;
  
  my %transcriptomes = ();
  my ($transcriptomes) = $record =~ /\s+transcript\s+sequences*\s+\((.*?)\)/ms;
  if ($transcriptomes) {
    my @transcriptomes = split(/,\s*|\s*and\s*/, $transcriptomes);
    %transcriptomes    = map { $_ => 1 } @transcriptomes;
  }
  
  my ($coords) = $record =~ /^PRIMARY[^\n]+\n(^\s+.+)/ms;
  $coords =~ s/\n^\S.+//ms;
  my @coords   = split(/\n/, $coords);
  my @edits    = ();
  
  # Note that it's important to get these in sequential order,
  # in order for the subsequent calculations to work.
  foreach my $coord (@coords) {
    my ($t_start, $t_end, $seq_name) =
      $coord =~ /^\s+(\d+)\-(\d+)\s+(\S+)\s+\d+\-\d+/;
    
    if (exists $transcriptomes{$seq_name} || $seq_name =~ /"N+"/) {
      push @edits, [$seq_name, $t_start, $t_end];
    }
  }
  
  return \@edits;
}

sub extract_seq {
  my ($self, $record) = @_;
  
  my ($seq) = $record =~ /^ORIGIN[^\n]+\n(.+)/ms;
  $seq =~ s/\n^\S.+//ms;
  $seq =~ s/\d+//gm;
  $seq =~ s/\s+//gm;
  
  return uc($seq);
}

sub adjust_five_prime_edit {
  my ($self, $record, $edits) = @_;
  my $adjustment = 0;
  my $seq_start;
  
  # The e! core code cannot deal with sequence added to the 5' end
  # that includes UTR; so remove the UTR, so that we at least get
  # the protein-coding bit right.
  # Also, if the inserted sequence starts mid-codon, lop off the partial codon.
  my ($cds_start, $cds_end) = $record =~ /^\s+CDS\D+(\d+)\D+(\d+)\s*$/m;
  my ($codon_start) = $record =~ /^\s+\/codon_start=(\d)\s*$/m;
  my $codon_adjustment = $codon_start - 1;
  
  # We assume that we will never have a two chunks of sequence
  # added via transcriptomic data added at the 5' end where the first
  # is UTR only.
  if ($cds_start > 1 && $$edits[0][1] == 1) {
    if ($cds_start <= $$edits[0][2]) {
      $adjustment = $codon_adjustment + $cds_start - 1;
      $seq_start = 1;
      
      foreach my $edit (@$edits) {
        $$edit[1] -= $adjustment;
        $$edit[2] -= $adjustment;
      }
      $$edits[0][1] = 1;
    } else {
      $adjustment = $$edits[0][2];
      
      foreach my $edit (@$edits) {
        $$edit[1] -= $adjustment;
        $$edit[2] -= $adjustment;
      }
      shift @$edits;
      
      my ($accession) = $record =~ /^VERSION\s+(\S+)/m;
      $self->warning("Transcriptomic sequence at 5' end of $accession is UTR only, edit will be ignored.");
    }
  } elsif ($codon_adjustment) {
    $adjustment = $codon_adjustment;
    $seq_start = 1;
    
    foreach my $edit (@$edits) {
      $$edit[1] -= $adjustment;
      $$edit[2] -= $adjustment;
    }
    $$edits[0][1] += $adjustment;
  }
  
  return ($adjustment, $seq_start);
}

sub adjust_three_prime_edit {
  my ($self, $record, $edits) = @_;
  my $adjustment = 0;
  
  # The e! core code cannot deal with sequence added to the 3' end
  # that includes UTR; so remove the UTR, so that we at least get
  # the protein-coding bit right.
  my ($cds_start, $cds_end) = $record =~ /^\s+CDS\D+(\d+)\D+(\d+)\s*$/m;
  
  # We assume that we will never have a two chunks of sequence
  # added via transcriptomic data added at the 3' end where the first
  # is UTR only.
  if ($cds_end < $$edits[-1][2]) {
    if ($cds_end > $$edits[-1][1]) {
      $adjustment = $$edits[-1][2] - $cds_end;
      $$edits[-1][2] = $cds_end;
    } else {
      my ($accession) = $record =~ /^VERSION\s+(\S+)/m;
      pop @$edits;
      $self->warning("Transcriptomic sequence at 3' end of $accession is UTR only, edit will be ignored.");
    }
  }
  
  return $adjustment;
}

sub seq_edits_from_protein {
  my ($self, $dba, $logic_name, $protein_fasta_file) = @_;
  
  my $ta = $dba->get_adaptor('Transcript');
  my $aa = $dba->get_adaptor("Attribute");
  
  my %protein = $self->load_fasta($protein_fasta_file);
  
  my $transcripts = $ta->fetch_all_by_logic_name($logic_name);
  
  my @transcript_ids = map {$_->stable_id} @$transcripts;
  
  foreach my $transcript_id (sort @transcript_ids) {
    my $transcript  = $ta->fetch_by_stable_id($transcript_id);
    my $translation = $transcript->translation;
    
    if ($translation && exists $protein{$translation->stable_id}) {
      my $db_seq   = $translation->seq;
      my $file_seq = $protein{$translation->stable_id};
      $file_seq =~ s/\*$//;
      
      if ($db_seq ne $file_seq) {
        my @file_seq = split(//, $file_seq);
        my @atts;
        
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

sub add_transcript_seq_edit {
  my ($self, $transcript, $start, $end, $seq) = @_;
  
  my $attribute = Bio::EnsEMBL::Attribute->new
  (
    -CODE  => '_rna_edit',
    -VALUE => "$start $end $seq",
  );
  
  # Need to add the attribute in order to force recalculation
  # of the sequence for the benefit of any subsequent edits.
  $transcript->add_Attributes($attribute);
  
  return $attribute;
}

sub add_translation_seq_edit {
  my ($self, $translation, $start, $end, $seq) = @_;
  
  my $code = 'amino_acid_sub';
  if ($seq eq 'U') {
    $code = '_selenocysteine';
  }
  
  my $attribute = Bio::EnsEMBL::Attribute->new
  (
    -CODE  => $code,
    -VALUE => "$start $end $seq",
  );
  
  # Need to add the attribute in order to force recalculation
  # of the sequence for the benefit of any subsequent edits.
  $translation->add_Attributes($attribute);
  
  return $attribute;
}

1;
