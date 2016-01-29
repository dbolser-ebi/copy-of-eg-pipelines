=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::FileDump::AlignmentDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Registry;
use Bio::PrimarySeq;
use Bio::SeqIO;

use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    seq_type     => 'aa',
  };
}

sub run {
  my ($self) = @_;
  my $seq_type = $self->param_required('seq_type');
  my $start_id = $self->param_required('start_id');
  my $end_id   = $self->param_required('end_id');
  my $tree_dir = $self->param_required('tree_dir');
  
  path($tree_dir)->mkpath;
  
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Multi', 'compara');
  
  my $gta = $dba->get_adaptor("GeneTree");
  my $trees = $gta->generic_fetch(
    "stable_id BETWEEN '$start_id' AND '$end_id'",
  );
  
  foreach my $tree (@$trees) {
    my $out_file = catdir($tree_dir, $tree->stable_id . '.fa');
    my $fasta = Bio::SeqIO->new(-format => 'Fasta', -file => ">$out_file");
    
    foreach my $leaf (@{$tree->get_all_leaves}) {      
      my ($seq, $alphabet);
      if ($seq_type eq 'cdna') {
        $seq = $leaf->alignment_string('cds');
        $alphabet = 'dna';
      } else {
        $seq = $leaf->alignment_string;
        $alphabet = 'protein';
      }
			$seq =~ s/\s+//g;
      
      my $seq_obj = Bio::PrimarySeq->new (
        -seq        => $seq,
        -primary_id => $leaf->stable_id,
        -display_id => $self->header($leaf),
        -alphabet   => $alphabet,
      );
      
      $fasta->write_seq($seq_obj);
    }
    
    $tree->release_tree;
  }
}

sub header {
  my ($self, $node) = @_;
  
  my $id = $node->stable_id;
  
  my $location = join(':',
    $node->dnafrag->name,
    $node->dnafrag_start,
    $node->dnafrag_end,
    $node->dnafrag_strand,
  );
    
  my $description = join('|',
    $node->taxonomy_level,
    $location,
    'gene:'.$node->gene_member->stable_id,
  );
  
  return "$id $description";
}

1;
