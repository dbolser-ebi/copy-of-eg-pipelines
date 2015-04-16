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

package Bio::EnsEMBL::EGPipeline::FileDump::UTR;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Feature');

sub new {
  my ($caller, %params) = @_;
  my $class = ref($caller) || $caller;
  
  my $slice = $params{'slice'};
  my $self = bless
  (
    {
      'slice'           => $slice,
      'seq_region_name' => $slice->seq_region_name,
      'start'           => $params{'start'},
      'end'             => $params{'end'},
      'strand'          => $params{'strand'},
      'source'          => $params{'source'},
      'parent_id'       => $params{'parent_id'},
    },
    $class
  );
  
  if ($params{'utr_type'} eq 'five_prime_UTR') {
    $self->SO_term('SO:0000204');
  } elsif ($params{'utr_type'} eq 'three_prime_UTR') {
    $self->SO_term('SO:0000205');
  }
  
  return $self;
}

sub SO_term {
  my ($self, $so_term) = @_;
  $self->{'SO_term'} = $so_term if defined $so_term;
  return $self->{'SO_term'};
}

sub source {
  my ($self, $source) = @_;
  $self->{'source'} = $source if defined $source;
  return $self->{'source'};
}

sub parent_id {
  my ($self, $parent_id) = @_;
  $self->{'parent_id'} = $parent_id if defined $parent_id;
  return $self->{'parent_id'};
}

sub summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->start;
  $summary{'end'}             = $self->end;
  $summary{'strand'}          = $self->strand;
  $summary{'source'}          = $self->source;
  $summary{'Parent'}          = $self->parent_id;
  
  return \%summary;
}

1;
