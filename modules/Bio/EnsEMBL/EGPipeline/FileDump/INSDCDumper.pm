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

package Bio::EnsEMBL::EGPipeline::FileDump::INSDCDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::GenomicFeatureDumper');

use Bio::SeqIO;

use Path::Tiny qw(path);

my @COMMENTS = (
  'This sequence displays annotation from <SOURCE> (<SOURCE_URL>), based on '.
  'underlying annotation from <PROVIDER> (<PROVIDER_URL>).',

  'All feature locations are relative to the first (5\') base of the sequence.'.
	'The sequence presented is always the forward strand of the assembly.',

  'The /gene indicates a unique id for a gene, /note="transcript_id=..." '.
	'a unique id for a transcript, /protein_id a unique id for a peptide '.
	'and /note="exon_id=..." a unique id for an exon. These ids are '.
	'maintained wherever possible between versions.',
);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'insdc_format'   => 'embl',
    'feature_type'   => ['Gene'],
    'data_type'      => 'basefeatures',
    'file_type'      => 'embl.dat',
    'gene_centric'   => 1,
    'per_chromosome' => 0,
    'source'         => 'VectorBase',
    'source_url'     => 'https://www.vectorbase.org',
  };
}

sub run {
  my ($self) = @_;
  my $species        = $self->param_required('species');
  my $db_type        = $self->param_required('db_type');
  my $out_file       = $self->param_required('out_file');
  my $feature_types  = $self->param_required('feature_type');
  my $per_chromosome = $self->param_required('per_chromosome');
  my $insdc_format   = $self->param_required('insdc_format');
  
  my $seqio = Bio::SeqIO->new(-file => ">$out_file", -format => $insdc_format);
  
  my $reg = 'Bio::EnsEMBL::Registry';
  my $sa = $reg->get_adaptor($species, $db_type, 'Slice');
  my $slices = $sa->fetch_all('toplevel');
  
  my %adaptors;
  foreach my $feature_type (@$feature_types) {
    $adaptors{$feature_type} = $reg->get_adaptor($species, $db_type, $feature_type);
  }
  my $mca = $sa->db->get_MetaContainer;
  
  my %chr;
  my $has_chromosomes = $self->has_chromosomes($sa);
  my $species_obj = $self->create_species_obj($mca);
  
  foreach my $slice (@$slices) {
    my $seq_obj = $self->create_seq_obj($slice, $species_obj);
    
    $self->add_secondary_acc($slice, $seq_obj);
    $self->add_annotations($mca, $seq_obj);
    $self->add_slice_feature($slice, $species_obj, $seq_obj);
    
    foreach my $feature_type (@$feature_types) {
      $self->add_features($feature_type, $adaptors{$feature_type}, $slice, $seq_obj);
    }
    $self->add_contig_features($slice, $seq_obj);
    
    $seqio->write_seq($seq_obj);
    
    if ($per_chromosome && $has_chromosomes) {
      my $chr_seqio = $self->chr_seqio($out_file, $slice, \%chr);
      $chr_seqio->write_seq($seq_obj);
    }
  }
  
  my $out_files = $self->param('out_files');
  foreach my $slice_name (keys %chr) {
    push @$out_files, $chr{$slice_name}{'file'};
  }
  $self->param('out_files', $out_files);
}

sub chr_seqio {
  my ($self, $out_file, $slice, $chr) = @_;
  
  my $insdc_format = $self->param_required('insdc_format');
  
  my $slice_name;
  if ($slice->karyotype_rank > 0) {
    $slice_name = 'chromosome.'.$slice->seq_region_name;
  } else {
    $slice_name = 'nonchromosomal';
  }
  
  unless (exists $$chr{$slice_name}) {
    (my $chr_file = $out_file) =~ s/([^\.]+)$/$slice_name.$1/;
    my $chr_seqio = Bio::SeqIO->new(-file => $chr_file, -format => $insdc_format);
    $$chr{$slice_name}{'seqio'} = $chr_seqio;
    $$chr{$slice_name}{'file'} = $chr_file;
  }
  
  return $$chr{$slice_name}{'seqio'};
}

1;
