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

Bio::EnsEMBL::EGPipeline::BlastAlignment::FilterHits

=head1 DESCRIPTION

Retain the top X hits across a genome, by selecting those with the
best E-value. Note that a single hit might map to multiple locations,
but each of those alignments corresponds to unique protein sequence.
In other words, the 'h' coordinates will not overlap.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::FilterHits;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'db_type'      => 'core',
    'filter_top_x' => 1,
    'blastx_top_x' => 10,
    'blastp_top_x' => 1,
    'create_gff'   => 0,
  };
}

sub run {
  my ($self) = @_;
  my $db_type      = $self->param_required('db_type');
  my $seq_type     = $self->param_required('seq_type');
  my $filter_top_x = $self->param_required('filter_top_x');
  my $blastp_top_x = $self->param_required('blastp_top_x');
  my $blastx_top_x = $self->param_required('blastx_top_x');
  my $logic_name   = $self->param_required('logic_name');
  
  if ($filter_top_x) {
    my $dba = $self->get_DBAdaptor($db_type);
    if ($seq_type eq 'genome') {
      $self->best_in_genome($dba, $blastx_top_x, $logic_name);
    } elsif ($seq_type eq 'proteome') {
      $self->best_in_proteome($dba, $blastp_top_x, $logic_name);
    }
  }
}

sub write_output {
  my ($self) = @_;
  my $create_gff = $self->param_required('create_gff');
  my $logic_name = $self->param_required('logic_name');
  
  if ($create_gff) {
    $self->dataflow_output_id({'logic_name' => [$logic_name]}, 2);
  }
}

sub best_in_genome {
  my ($self, $dba, $top_x, $logic_name) = @_;
  
  my %grouped_features;
  
  my $sa = $dba->get_adaptor("Slice");
  my $pafa = $dba->get_adaptor("ProteinAlignFeature");
  
  my $slices = $sa->fetch_all('toplevel');
  foreach my $slice (@$slices) {
    my %pafs;
    my $pafs = $pafa->fetch_all_by_Slice($slice, $logic_name);
    foreach my $paf (@$pafs) {
      push @{$pafs{$paf->hseqname}}, $paf;
    }
    
    # For hits from the same source on a single seq_region, we want
    # to group them into sets within which the hits do not overlap.
    # This means that short hits that would otherwise be lost due to
    # E-value filtering are retained, by piggy-backing on the higher
    # scoring hits.
    foreach my $hit_name (keys %pafs) {
      my @pafs = sort {$a->p_value <=> $b->p_value or $b->score <=> $a->score} @{$pafs{$hit_name}};
      my @groups = ();
      $self->group_features(\@pafs, \@groups);
      push @{$grouped_features{$hit_name}}, @groups;
    }
  }
  
  $self->remove_features($top_x, $pafa, \%grouped_features);
}

sub best_in_proteome {
  my ($self, $dba, $top_x, $logic_name) = @_;
  
  my %grouped_features;
  
  # The protein feature adaptor has no 'fetch_all' method,
  # so we have to get them by iterating over translations.
  my $ta  = $dba->get_adaptor("Translation");
  my $pfa = $dba->get_adaptor("ProteinFeature");
  
  my $ts = $ta->fetch_all();
  
  my %pfs;
  foreach my $t (@$ts) {
    my $pfs = $pfa->fetch_all_by_translation_id($t->dbID);
    foreach my $pf (@$pfs) {
      if ($pf->analysis->logic_name eq $logic_name) {
        push @{$pfs{$pf->hseqname}}, $pf;
      }
    }
    
    foreach my $hit_name (keys %pfs) {
      my @pfs = sort {$a->p_value <=> $b->p_value or $b->score <=> $a->score} @{$pfs{$hit_name}};    
      my @groups = ();
      $self->group_features(\@pfs, \@groups);
      push @{$grouped_features{$hit_name}}, @groups;
    }
  }
  
  $self->remove_features($top_x, $pfa, \%grouped_features);
}

sub group_features {
  my ($self, $features, $groups) = @_;
  
  while (scalar(@$features)) {
    my $group = [];
    $self->disjoint($features, $features, $group);
    
    push @$groups, $group;
  }
} 

sub disjoint {
  my ($self, $features, $features_subset, $group) = @_;
  
  my $best_hit = shift @$features_subset;
  
  if (scalar(@$group)) {
    my $offset;
    for (my $i=0; $i < scalar(@$features); $i++) {
      if ($$features[$i]->dbID == $best_hit->dbID) {
        $offset = $i;
        last;
      }
    }
    splice(@$features, $offset, 1);
  }
  
  push @$group, $best_hit;
  
  my @disjoint = ();
  if (scalar(@$features_subset)) {
    foreach my $feature (@$features_subset) {
      if (! $self->hit_overlap($best_hit, $feature)) {
        push @disjoint, $feature;
      }
    }
    
    if (scalar(@disjoint)) {
      $self->disjoint($features, \@disjoint, $group);
    }
  }
}

sub remove_features {
  my ($self, $top_x, $adaptor, $grouped_features) = @_;
  
  # Because the features within each group are sorted such that the best
  # hit is first, we can sort the groups by sorting on the first element
  # of each. The top X are then removed from the front of the list of
  # grouped features, and whatever is left hasn't made the cut and is deleted.
  
  foreach my $hit_name (keys %$grouped_features) {
    my @groups = sort {$a->[0]->p_value <=> $b->[0]->p_value or $b->[0]->score <=> $a->[0]->score} @{$$grouped_features{$hit_name}};
    splice(@groups, 0, $top_x);
    foreach my $group (@groups) {
      foreach my $feature (@$group) {
        if ($adaptor->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor')) {
          my $sth = $adaptor->prepare("DELETE FROM protein_feature WHERE protein_feature_id = ?");
          $sth->execute($feature->dbID);
        } else {
          $adaptor->remove($feature);
        }
      }
    }
  }
}

sub hit_overlap {
  my ($self, $feature1, $feature2) = @_;
  
  # The important thing to note here is that we're checking
  # for overlapping _hit_ coordinates, not seq_region coordinates.
  
  my ($s1, $e1) = ($feature1->hstart, $feature1->hend);
  my ($s2, $e2) = ($feature2->hstart, $feature2->hend);
  
	my $overlap = 0;
  
	if ($s1 >= $s2 && $e1 <= $e2) {
		$overlap = 1;
	} elsif ($s1 <= $s2 && $e1 >= $e2) {
		$overlap = 2;
	} elsif ($e1 >= $s2 && $e1 <= $e2) {
		$overlap = 3;
	} elsif ($s1 >= $s2 && $s1 <= $e2) {
		$overlap = 4;
	}
  
	return $overlap;
}

1;
