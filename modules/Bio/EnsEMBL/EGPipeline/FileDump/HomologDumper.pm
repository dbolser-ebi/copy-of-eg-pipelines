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

package Bio::EnsEMBL::EGPipeline::FileDump::HomologDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::Graph::OrthoXMLWriter;

use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'compara'          => 'multi',
    'homolog_format'   => 'xml',
    'homolog_source'   => 'VectorBase',
    'use_gene_version' => 1,
  };
}

sub run {
  my ($self) = @_;
  my $compara          = $self->param_required('compara');
  my $homolog_format   = $self->param_required('homolog_format');
  my $homolog_source   = $self->param_required('homolog_source');
  my $use_gene_version = $self->param_required('use_gene_version');
  my $release_date     = $self->param_required('release_date');
  my $start_id         = $self->param_required('start_id');
  my $end_id           = $self->param_required('end_id');
  my $tree_dir         = $self->param_required('tree_dir');
  
  path($tree_dir)->mkpath;
  
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($compara, 'compara');
  
  my $ha = $dba->get_adaptor("Homology");
  my $gta = $dba->get_adaptor("GeneTree");
  my $trees = $gta->generic_fetch(
    "stable_id BETWEEN '$start_id' AND '$end_id'",
  );
  
  # Compara db has the geneset start_date, whereas VB uses the version
  my $gene_sets;
  if ($use_gene_version) {
    $gene_sets = $self->gene_sets($dba);
  }
  
  my @out_files;
  
  foreach my $tree (@$trees) {
    my $homologies = $ha->fetch_all_by_tree_node_id($tree->root_id);
    
    if ($homolog_format eq 'xml') {
      my $out_file = catdir($tree_dir, $tree->stable_id . '.xml');
      
      my $writer = Bio::EnsEMBL::Compara::Graph::OrthoXMLWriter->new(
        -SOURCE         => $homolog_source,
        -SOURCE_VERSION => "$homolog_source $release_date",
        -FILE           => $out_file,
      );
      $writer->write_homologies($homologies);
      $writer->finish();
      $writer->handle()->close();
      
      if ($use_gene_version) {
        $self->replace_gene_sets($out_file, $gene_sets);
      }
      
      push @out_files, $out_file;
      
    } else {
      $self->throw("Unrecognised homology format '$homolog_format'");
    }
    
    $tree->release_tree;
  }
  
  $self->param('out_files', \@out_files);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $out_file (@{$self->param('out_files')}) {
    $self->dataflow_output_id({out_file => $out_file}, 2);
  }
}

sub gene_sets {
  my ($self, $dba) = @_;
  
  my %gene_sets;
  my $gdba = $dba->get_adaptor("GenomeDB");
  
  my $all_genome_dbs = $gdba->fetch_all();
  foreach my $genome_db (@$all_genome_dbs) {
    my $species = $genome_db->name;
    my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
    
    my $mc = $core_dba->get_MetaContainer;
    my $assembly = $mc->single_value_by_key('assembly.default');
    my $start_date = $mc->single_value_by_key('genebuild.start_date');
    my $gene_set = $mc->single_value_by_key('genebuild.version');
    
    $gene_sets{"$assembly/$start_date"} = $gene_set;
    
    $core_dba->dbc->disconnect_if_idle();
  }
  
  return \%gene_sets;
}

sub replace_gene_sets {
  my ($self, $out_file, $gene_sets) = @_;
  
  my $file = path($out_file);
  my $data = $file->slurp;
  
  foreach my $id (keys %$gene_sets) {
    my $new_id = $$gene_sets{$id};
    $data =~ s/$id/$new_id/gm;
  }
  
  $file->spew($data);
}

1;
