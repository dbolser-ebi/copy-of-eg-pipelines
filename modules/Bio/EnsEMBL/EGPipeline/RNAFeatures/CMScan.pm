=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::CMScan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::CMScan;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'rfam_trna'      => 0,
    'rfam_rrna'      => 1,
    'rfam_blacklist' => [],
  };
}

sub fetch_runnable {
  my ($self) = @_;
  
  my %parameters;
  if (%{$self->param('parameters_hash')}) {
    %parameters = %{$self->param('parameters_hash')};
  }
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::CMScan->new
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
  
  $self->param('save_object_type', 'DnaAlignFeature');
  
  return $runnable;
}

sub results_by_index {
  my ($self, $results) = @_;
  
  $results =~ s/(^#.*\n^[^#]+^#.*\n)//m;
  my $header = $1;
  
  my @lines = split(/\n/, $results);
  
  my @results;
  my $current;
  foreach my $line (@lines) {
    if ($line =~ /^\S/) {
      push(@results, $current) if $current;
      if ($line =~ /^>>/) {
        $current = "$line\n";
      } else {
        $current = undef;
      }
    } elsif ($current) {
      $current .= "$line\n";
    }
  }
  push(@results, $current) if $current;
  
  my %seqnames;
  foreach my $result (@results) {
    my ($seqname) = $result =~ /CS\n.*\n.*\n\s*(\S+)/m;
    $seqnames{$seqname}{'result'} .= "$result\n";
  }
  foreach my $seqname (keys %seqnames) {
    $seqnames{$seqname}{'header'} = $header;
  }
  
  return %seqnames;
}

sub filter_output {
  my ($self, $runnable) = @_;
  
  if ($runnable->db_name eq 'RFAM') {
    my $rfam_trna      = $self->param_required('rfam_trna');
    my $rfam_rrna      = $self->param_required('rfam_rrna');
    my $rfam_blacklist = $self->param_required('rfam_blacklist');
    my %rfam_blacklist = map { $_ => 1 } @$rfam_blacklist;
    
    my @filtered;
    
    foreach my $feature (@{$runnable->output}) {
      my ($acc, $biotype) = $feature->extra_data =~ /Accession=([^;]+).+Biotype=([^;]+)/;
      
      if (exists $rfam_blacklist{$acc}) {
        $self->warning("Skipping blacklisted accession $acc on sequence ".$runnable->query->name);
      } else {
        if (! $rfam_trna && $biotype eq 'tRNA') {
          $self->warning("Skipping tRNA $acc on sequence ".$runnable->query->name);
        } elsif (! $rfam_rrna && $biotype eq 'rRNA') {
          $self->warning("Skipping rRNA $acc on sequence ".$runnable->query->name);
        } else {
          push @filtered, $feature;
        }
      }
    }
    
    # Output is cumulative, so need to manually erase existing results.
    # (Note that calling the runnable's 'output' method will NOT work.
    $runnable->{'output'} = \@filtered;
  }
}

1;
