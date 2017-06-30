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

package Bio::EnsEMBL::EGPipeline::FileDump::TreeFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Registry;
use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'compara'          => 'multi',
    'skip_dumps'       => [],
    'files_per_subdir' => 500,
  };
}

sub write_output {
  my ($self) = @_;
  my $compara          = $self->param_required('compara');
  my $dump_types       = $self->param_required('dump_types');
  my $skip_dumps       = $self->param_required('skip_dumps');
  my $results_dir      = $self->param_required('results_dir');
  my $files_per_subdir = $self->param_required('files_per_subdir');
  
  my %skip_dumps = map { $_ => 1 } @$skip_dumps;
  
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($compara, 'compara');
  
  my $gta = $dba->get_adaptor("GeneTree");
  my $all_trees = $gta->fetch_all(
    -tree_type     => 'tree',
    -clusterset_id => 'default',
  );
  
  my @all_trees = sort { $a->stable_id cmp $b->stable_id } @$all_trees;
  my $start_id = $all_trees[0]->stable_id;
  my $end_id;
  my $counter = 1;
  my @tree_groups;
  
  foreach my $tree (@all_trees) {
    if ($counter > $files_per_subdir) {
      push @tree_groups,
        {
          'start_id' => $start_id,
          'end_id'   => $end_id,
          'id_range' => $self->id_range($start_id, $end_id),
        };
      
      $start_id = $tree->stable_id;
      $counter = 1;
    } else {
      $end_id = $tree->stable_id;
      $counter++;
    }
    
    $tree->release_tree;
  }
  
  foreach my $flow (keys %$dump_types) {
    foreach my $dump_type (@{$$dump_types{$flow}}) {
      if (!exists $skip_dumps{$dump_type}) {
        my $sub_dir = $self->sub_dir($dump_type);
        
        foreach my $tree_group (@tree_groups) {
          my %output_ids = %$tree_group;
          
          my $tree_dir = catdir($results_dir, $sub_dir, $$tree_group{'id_range'});
          path($tree_dir)->mkpath;
          $output_ids{'tree_dir'} = $tree_dir;
          
          $self->dataflow_output_id(\%output_ids, $flow);
        }
        
        $self->dataflow_output_id(
          { out_dir => $results_dir , sub_dir => $sub_dir }, 1);
      }
    }
  }
}

sub id_range {
  my ($self, $start_id, $end_id) = @_;
  
  my @start_id = split(//, $start_id);
  my @end_id = split(//, $end_id);
  
  my @range;
  
  for (my $i=0; $i<scalar(@start_id); $i++) {
    if ($start_id[$i] eq $end_id[$i]) {
      push @range, $start_id[$i];
    } else {
      push @range, @start_id[$i .. (scalar(@start_id) - 1)];
      push @range, '-';
      push @range, @end_id[$i .. (scalar(@end_id) - 1)];
      last;
    }
  }
  
  return join('', @range);
}

sub sub_dir {
  my ($self, $dump_type) = @_;
  my $compara_dumps = $self->param_required('compara_dumps');
  my $release_date  = $self->param_required('release_date');
  return $$compara_dumps{$dump_type}."_$release_date";
}

1;
