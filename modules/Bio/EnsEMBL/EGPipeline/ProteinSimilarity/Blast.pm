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

Bio::EnsEMBL::EGPipeline::ProteinSimilarity::Blast

=head1 DESCRIPTION

Run BLAST, using the Ensembl modules for parsing and filtering.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::Blast;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BlastEG;
use Bio::EnsEMBL::Analysis::Tools::FeatureFilter;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use Bio::EnsEMBL::ProteinFeature;

use File::Path qw(make_path remove_tree);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'blast_type'       => 'ncbi',
    'output_regex'     => '^\s*(\w+)',
    'query_type'       => 'dna',
    'database_type'    => 'pep',
    'pvalue_threshold' => 0.01,
  };
}

sub fetch_runnable {
  my $self = shift @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  if ($self->param_required('blast_type') eq 'wu') {
    $ENV{BLASTMAT} = $self->param_required('blast_matrix');
  }
  
  $parameters{'TYPE'}     = $self->param_required('blast_type');
  $parameters{'DATABASE'} = $self->param_required('blast_db');
  $parameters{'PARSER'}   = $self->make_parser();
  $parameters{'FILTER'}   = undef;
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastEG->new
  (
    -query    => $self->param('query'),
    -program  => $self->param('program'),
    -analysis => $self->param('analysis'),
    -datadir  => $self->param('datadir'),
    -bindir   => $self->param('bindir'),
    -libdir   => $self->param('libdir'),
    -workdir  => $self->param('workdir'),
    %parameters,
  );
  
  if ($self->param_required('database_type') eq 'pep') {
    if ($self->param_required('query_type') eq 'pep') {
      $self->param('results_index', 'translation');
      $self->param('save_object_type', 'ProteinFeature');
    } else {
      $self->param('save_object_type', 'ProteinAlignFeature');
    }
  } else {
    $self->param('save_object_type', 'DnaAlignFeature');
  }
  
  return $runnable;
}

sub results_by_index {
  my ($self, $results) = @_;
  my %seqnames;
  
  my @results = split(/Query=\s+/, $results);
  my $header = shift @results;
  foreach my $result (@results) {
    my ($seqname) = $result =~ /^\s*(\S+)/;
    #$result =~ s/\n+\z//m;
    $seqnames{$seqname}{'result'} = "Query= $result";
    $seqnames{$seqname}{'header'} = $header;
  }
    
  return %seqnames;
}

sub filter_output {
  my ($self, $runnable) = @_;
  
  my $filter = Bio::EnsEMBL::Analysis::Tools::FeatureFilter->new
    (
      -min_score => 200,
      -prune => 1,
    );
  
  my $filtered = $filter->filter_results($runnable->output);
  
  # Output is cumulative, so need to manually erase existing results.
  # (Note that calling the runnable's 'output' method will NOT work.
  $runnable->{'output'} = $filtered;
}

sub post_processing {
  my ($self, $runnable) = @_;
  
  my @pfs;
  if ($self->param('database_type') eq 'pep') {
    if ($self->param('query_type') eq 'pep') {
      foreach my $feature (@{$runnable->output}) {
        my $pf = Bio::EnsEMBL::ProteinFeature->new(
          -translation_id => $runnable->query->dbID,
          -seqname        => $runnable->query->dbID,
          -slice          => $runnable->query->transcript->slice,
          -start          => $feature->start,
					-end            => $feature->end,
					-hseqname       => $feature->hseqname,
					-hstart         => $feature->hstart,
					-hend           => $feature->hend,
          -percent_id     => $feature->percent_id,
          -score          => $feature->score,
          -p_value        => $feature->p_value,
          -analysis       => $feature->analysis
        );
        push @pfs, $pf;
      }
    }
  }
  
  # Output is cumulative, so need to manually erase existing results.
  # (Note that calling the runnable's 'output' method will NOT work.
  $runnable->{'output'} = \@pfs;
}

sub make_parser {
  my ($self) = @_;
  
  my $parser = Bio::EnsEMBL::Analysis::Tools::FilterBPlite->new
    (
      -regex          => $self->param('output_regex'),
      -query_type     => $self->param('query_type'),
      -database_type  => $self->param('database_type'),
      -threshold_type => 'PVALUE',
      -threshold      => $self->param('pvalue_threshold'),
    );
  
  return $parser;
}

1;
