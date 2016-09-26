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

Bio::EnsEMBL::EGPipeline::BlastAlignment::ExtractSpeciesUniprot

=head1 DESCRIPTION

Extract UniProt entries for a particular species, from Fasta format file.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::ExtractSpeciesUniprot;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;

sub extract_species {
  my ($self, $file, $output_file, $species) = @_;

  my $seq_out = Bio::SeqIO->new(
    -file   => '>'.$output_file,
    -format => 'fasta',
  );

  $species =~ s/_/ /g;
  $species =~ s/[A-Z]//g;
  $species = ucfirst($species);

  open(F, $file);
  my $seq_in = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );

  while (my $inseq = $seq_in->next_seq) {
    if ($inseq->desc =~ /OS=$species/) {
      my (undef, $primary_id, $secondary_id) = split(/\|/, $inseq->display_id);
      $inseq->display_id($primary_id);

      my ($desc, $version) = $inseq->desc =~ /^(.*)\s+OS=.*SV=(\d+)/;
      $inseq->desc(join('|', map { $_ || '' } ($secondary_id, $desc, $version)));

      $seq_out->write_seq($inseq);
    }
  }
}

1;
