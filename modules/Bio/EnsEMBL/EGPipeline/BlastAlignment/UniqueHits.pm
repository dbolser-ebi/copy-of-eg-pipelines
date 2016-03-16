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

Bio::EnsEMBL::EGPipeline::BlastAlignment::UniqueHits

=head1 DESCRIPTION

Make hits unique within a seq_region, by retaining only the one with the
best E-value. Note that the best hit might map to multiple locations,
but each of those alignments corresponds to unique protein sequence. In other
words, the 'h' coordinates will not overlap.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::UniqueHits;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'db_type'     => 'core',
    'unique_hits' => 1,
    'create_gff'  => 0,
  };
}

sub run {
  my ($self) = @_;
  my $db_type     = $self->param_required('db_type');
  my $seq_type    = $self->param_required('seq_type');
  my $unique_hits = $self->param_required('unique_hits');
  my $logic_name  = $self->param_required('logic_name');
  
  if ($unique_hits) {
    my $dba = $self->get_DBAdaptor($db_type);
    if ($seq_type eq 'genome') {
      $self->unique_seq_region($dba, $logic_name);
    } elsif ($seq_type eq 'proteome') {
      $self->unique_protein_feature($dba, $logic_name);
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

sub unique_seq_region {
  my ($self, $dba, $logic_name) = @_;
  
  my $sa = $dba->get_adaptor("Slice");
  my $pafa = $dba->get_adaptor("ProteinAlignFeature");

  my $slices = $sa->fetch_all('toplevel');
  foreach my $slice (@$slices) {
    my $pafs = $pafa->fetch_all_by_Slice($slice, $logic_name);
    
    my %pafs;
    foreach my $paf (@$pafs) {
      push @{$pafs{$paf->hseqname}}, $paf;
    }
    
    foreach my $hit_name (keys %pafs) {
      my @pafs = sort {$a->p_value <=> $b->p_value or $b->score <=> $a->score} @{$pafs{$hit_name}};    
      $self->remove_overlapping($pafa, \@pafs);
    }
  }
}

sub unique_protein_feature {
  my ($self, $dba, $logic_name) = @_;
  
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
  }
  
  foreach my $hit_name (keys %pfs) {
    my @pfs = sort {$a->p_value <=> $b->p_value or $b->score <=> $a->score} @{$pfs{$hit_name}};    
    $self->remove_overlapping($pfa, \@pfs);
  }
}

sub remove_overlapping {
  my ($self, $adaptor, $features) = @_;
  
  my $best_hit = shift @$features;
  
  my @non_overlapping = ();
  foreach my $feature (@$features) {
    if ($self->hit_overlap($best_hit, $feature)) {
      if !($best_hit->p_value == $feature->p_value && $best_hit->score == $feature->score) {      
        if ($adaptor->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor')) {
          my $sth = $adaptor->prepare("DELETE FROM protein_feature WHERE protein_feature_id = ?");
          $sth->execute($feature->dbID);
        } else {
          $adaptor->remove($feature);
        }
      }
    } else {
      push @non_overlapping, $feature;
    }
  }
  
  if (scalar(@non_overlapping)) {
    $self->remove_overlapping(\@non_overlapping);
  }
}

sub hit_overlap {
  my ($self, $best_hit, $feature) = @_;
  my ($s1, $e1) = ($best_hit->hstart, $best_hit->hend);
  my ($s2, $e2) = ($feature->hstart, $feature->hend);
  
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
