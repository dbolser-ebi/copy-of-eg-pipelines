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

package Bio::EnsEMBL::EGPipeline::FileDump::SeqRegion;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Feature');

sub new {
  my ($caller, $slice, $source) = @_;
  my $class = ref($caller) || $caller;
  
  my $self = bless
  (
    {
      'slice'   => $slice,
      'seqname' => $slice->seq_region_name,
      'start'   => $slice->start,
      'end'     => $slice->end,
      'strand'  => 0,
      'source'  => $source || $slice->coord_system->version,
    },
    $class
  );
  
  my $seq_region_type = $slice->coord_system->name;
  if ($seq_region_type eq 'chromosome') {
    $self->SO_term('SO:0000340');
  } elsif ($seq_region_type =~ /^(supercontig|scaffold)$/) {
    $self->SO_term('SO:0000148');
  } elsif ($seq_region_type eq 'contig') {
    $self->SO_term('SO:0000149');
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

sub summary_as_hash {
  my $self = shift;
	my %summary;
  
  my @aliases = map { $_->name } @{$self->slice->get_all_synonyms()};
  
	$summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $self->source;
	$summary{'start'}           = $self->seq_region_start;
	$summary{'end'}             = $self->seq_region_end;
	$summary{'strand'}          = $self->seq_region_strand;
  $summary{'id'}              = $self->seq_region_name;
  $summary{'Is_circular'}     = $self->slice->is_circular ? 'true' : undef;
  $summary{'Alias'}           = \@aliases if scalar(@aliases);
  
	return \%summary;
}

1;
