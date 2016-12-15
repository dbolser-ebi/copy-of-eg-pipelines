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

package Bio::EnsEMBL::EGPipeline::FileDump::BEDDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Bio::EnsEMBL::Utils::IO::BEDSerializer;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'feature_type'       => ['Transcript'],
    'data_type'          => 'basefeatures',
    'file_type'          => 'bed',
    'per_chromosome'     => 0,
    'logic_name'         => [],
  };
}

sub run {
  my ($self) = @_;
  my $species        = $self->param_required('species');
  my $db_type        = $self->param_required('db_type');
  my $out_file       = $self->param_required('out_file');
  my $feature_types  = $self->param_required('feature_type');
  my $per_chromosome = $self->param_required('per_chromosome');
  my $logic_names    = $self->param_required('logic_name');
  
  my $reg = 'Bio::EnsEMBL::Registry';
  
  my $sa = $reg->get_adaptor($species, $db_type, 'Slice');
  my $slices = $sa->fetch_all('toplevel');
  
  open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
  my $serializer = Bio::EnsEMBL::Utils::IO::BEDSerializer->new($out_fh);
  
  my %adaptors;
  foreach my $feature_type (@$feature_types) {
    $adaptors{$feature_type} = $reg->get_adaptor($species, $db_type, $feature_type);
  }
  
  my %chr;
  my $has_chromosomes = $self->has_chromosomes($sa);
  
  foreach my $slice (@$slices) {    
    foreach my $feature_type (@$feature_types) {
      my $features = $self->fetch_features($feature_type, $adaptors{$feature_type}, $logic_names, $slice);
      $serializer->print_feature_list($features);
      
      if ($per_chromosome && $has_chromosomes) {
        my $chr_serializer = $self->chr_serializer($slice, \%chr);
        $chr_serializer->print_feature_list($features);
      }
    }
  }
  
  close($out_fh);
  my $out_files = $self->param('out_files');
  
  foreach my $slice_name (keys %chr) {
    close($chr{$slice_name}{'fh'});
    push @$out_files, $chr{$slice_name}{'file'};
  }
    
  $self->param('out_files', $out_files);
}

sub fetch_features {
  my ($self, $feature_type, $adaptor, $logic_names, $slice) = @_;
  
  my @features;
  if (scalar(@$logic_names) == 0) {
    @features = @{$adaptor->fetch_all_by_Slice($slice)};
  } else {
    foreach my $logic_name (@$logic_names) {
      my $features;
      if ($feature_type eq 'Transcript') {
        $features = $adaptor->fetch_all_by_Slice($slice, 0, $logic_name);
      } else {
        $features = $adaptor->fetch_all_by_Slice($slice, $logic_name);
      }
      push @features, @$features;
    }
  }
  
  return \@features;
}

sub chr_serializer {
  my ($self, $slice, $chr) = @_;
  
  my $out_file = $self->param_required('out_file');
  
  my $slice_name;
  if ($slice->karyotype_rank > 0) {
    $slice_name = 'chromosome.'.$slice->seq_region_name;
  } else {
    $slice_name = 'nonchromosomal';
  }
    
  unless (exists $$chr{$slice_name}) {
    (my $chr_file = $out_file) =~ s/([^\.]+)$/$slice_name.$1/;
    open(my $chr_fh, '>', $chr_file) or $self->throw("Cannot open file $chr_file: $!");
    
    my $chr_serializer = Bio::EnsEMBL::Utils::IO::BEDSerializer->new($chr_fh);
    $chr_serializer->print_main_header([$slice]);
    
    $$chr{$slice_name}{'fh'} = $chr_fh;
    $$chr{$slice_name}{'file'} = $chr_file;
    $$chr{$slice_name}{'serializer'} = $chr_serializer;
  }
  
  return $$chr{$slice_name}{'serializer'};
}

1;
