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

package Bio::EnsEMBL::EGPipeline::FileDump::TreeDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::Graph::GeneTreePhyloXMLWriter;

use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'tree_format' => 'newick',
    'tree_source' => 'VectorBase',
    'seq_type'    => 'aa',
  };
}

sub run {
  my ($self) = @_;
  my $tree_format = $self->param_required('tree_format');
  my $tree_source = $self->param_required('tree_source');
  my $seq_type    = $self->param_required('seq_type');
  my $start_id    = $self->param_required('start_id');
  my $end_id      = $self->param_required('end_id');
  my $sub_dir     = $self->param_required('sub_dir');
  
  path($sub_dir)->mkpath;
  
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Multi', 'compara');
  
  my $gta = $dba->get_adaptor("GeneTree");
  my $trees = $gta->generic_fetch(
    "stable_id BETWEEN '$start_id' AND '$end_id'",
  );
  
  foreach my $tree (@$trees) {
    if ($tree_format eq 'newick') {
      my $out_file = catdir($sub_dir, $tree->stable_id . '.nh');
      
      my $file = path($out_file);
      $file->spew($tree->newick_format('simple'));
      
    } elsif ($tree_format eq 'xml') {
      my $out_file = catdir($sub_dir, $tree->stable_id . '.xml');
      
      my $writer = Bio::EnsEMBL::Compara::Graph::GeneTreePhyloXMLWriter->new(
        -SOURCE  => $tree_source,
        -ALIGNED => 1,
        -CDNA    => $seq_type eq 'cdna' ? 1 : 0,
        -FILE    => $out_file,
      );
      $writer->write_trees($tree);
      $writer->finish();
      
    } else {
      $self->throw("Unrecognised tree format '$tree_format'");
    }
    
    $tree->release_tree;
  }
}

1;
