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

package Bio::EnsEMBL::EGPipeline::FileDump::GFF3Dumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Bio::EnsEMBL::Utils::IO::GFFSerializer;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::RepeatFeature;
use EnsEMBL::REST::EnsemblModel::CDS;
use EnsEMBL::REST::EnsemblModel::ExonTranscript;
use Bio::EnsEMBL::EGPipeline::FileDump::SeqRegion;
use Bio::EnsEMBL::EGPipeline::FileDump::UTR;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'feature_types'    => ['Gene', 'Transcript'],
    'include_scaffold' => 1,
    'data_type'        => 'basefeatures',
    'file_type'        => 'gff3',
    'logic_name'       => undef,
    'per_chromosome'   => 0,
  };
}

sub run {
  my ($self) = @_;
  my $species          = $self->param_required('species');
  my $db_type          = $self->param_required('db_type');
  my $out_file         = $self->param_required('out_file');
  my $out_fh           = $self->param_required('out_fh');
  my $feature_types    = $self->param_required('feature_types');
  my $include_scaffold = $self->param_required('include_scaffold');
  my $logic_name       = $self->param('logic_name');
  my $per_chromosome   = $self->param('per_chromosome');
  
  my $reg = 'Bio::EnsEMBL::Registry';

  my $sa = $reg->get_adaptor($species, $db_type, 'Slice');
  my $slices = $sa->fetch_all('toplevel');
  
  my $oa = $reg->get_adaptor('multi', 'ontology', 'OntologyTerm');
  
  my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($oa, $out_fh);
  $serializer->print_main_header($slices);
  
  my $mc = $self->core_dba->get_MetaContainer();
  my $provider = $mc->single_value_by_key('provider.name');
  my $assembly = $mc->single_value_by_key('assembly.default');
  $serializer->print_metadata("#genome-build $provider $assembly");
  
  my %adaptors;
  foreach my $feature_type (@$feature_types) {
    $adaptors{$feature_type} = $reg->get_adaptor($species, $db_type, $feature_type);
    if ($feature_type eq 'Transcript') {
      $adaptors{'Exon'} = $reg->get_adaptor($species, $db_type, 'Exon');
    }
  }
  
  my %chr;
  my $has_chromosome = $self->has_chromosome($sa);
  
  foreach my $slice (@$slices) {
    if ($include_scaffold) {
      my $feature = Bio::EnsEMBL::EGPipeline::FileDump::SeqRegion->new($slice, $provider);
      $serializer->print_feature($feature);
    }
    
    foreach my $feature_type (@$feature_types) {
      my $features = $adaptors{$feature_type}->fetch_all_by_Slice($slice, $logic_name);
      if ($feature_type eq 'Transcript') {
        my $exons = $adaptors{'Exon'}->fetch_all_by_Slice($slice);
        $features = $self->transcript_features($features, $exons);
      }
      $serializer->print_feature_list($features);
      
      if ($per_chromosome && $has_chromosome) {
        my $chr_serializer = $self->chr_serializer($out_file, $oa, $slice, \%chr);
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
  
  $self->param('out_files', $out_files)
}

sub transcript_features {
  my ($self, $features, $exons) = @_;
  
  my $exon_features = EnsEMBL::REST::EnsemblModel::ExonTranscript->build_all_from_Exons($exons);
  my $cds_features = EnsEMBL::REST::EnsemblModel::CDS->new_from_Transcripts($features);
  my @utr_features;
  foreach my $transcript (@$features) {
    my $five_prime = $transcript->five_prime_utr_Feature;
    if ($five_prime) {
      my $utrs = $self->utr_features($five_prime, 'five_prime_UTR', $transcript);
      push @utr_features, @$utrs;
    }
    
    my $three_prime  = $transcript->three_prime_utr_Feature;
    if ($three_prime) {
      my $utrs = $self->utr_features($three_prime, 'three_prime_UTR', $transcript);
      push @utr_features, @$utrs;
    }
  }
  
  return [@$features, @$exon_features, @$cds_features, @utr_features];
}

sub utr_features {
  my ($self, $utr, $utr_type, $transcript) = @_;
  my @utrs;
  
  my $exons = $transcript->get_all_Exons();
  
  foreach my $exon (@$exons) {
    my $strand = $exon->strand;
    
    my %params = (
      slice     => $exon->slice,
      strand    => $strand,
      source    => $transcript->source,
      parent_id => $transcript->stable_id,
      utr_type  => $utr_type,
    );
    
    if ($utr->start <= $exon->start && $utr->end >= $exon->end) {
      # Whole exon is UTR
      $params{start} = $exon->start;
      $params{end} = $exon->end;
      
    } else {
      if ($utr_type eq 'five_prime_UTR') {
        if ($strand == -1) {
          if (is_between($utr->start, $exon->start, $exon->end)) {
            $params{start} = $utr->start;
            $params{end} = $exon->end;
          }
        } else {
          if (is_between($utr->end, $exon->start, $exon->end)) {
            $params{start} = $exon->start;
            $params{end} = $utr->end;
          }
        }
      } elsif ($utr_type eq 'three_prime_UTR') {
        if ($strand == -1) {
          if (is_between($utr->end, $exon->start, $exon->end)) {
            $params{start} = $exon->start;
            $params{end} = $utr->end;
          }
        } else {
          if (is_between($utr->start, $exon->start, $exon->end)) {
            $params{start} = $utr->start;
            $params{end} = $exon->end;
          }
        }
      }
    }
    
    if (defined $params{start} && defined $params{end}) {
      my $utr = Bio::EnsEMBL::EGPipeline::FileDump::UTR->new(%params);
      push @utrs, $utr;
    }
  }
  
  return \@utrs;
}

sub is_between {
  my ($a, $x, $y) = @_;
  
  if ($a >= $x && $a <= $y) {
    return 1;
  } else {
    return 0;
  }
}

sub chr_serializer {
  my ($self, $out_file, $oa, $slice, $chr) = @_;
  
  my $slice_name;
  if ($slice->karyotype_rank > 0) {
    $slice_name = 'chromosome.'.$slice->seq_region_name;
  } else {
    $slice_name = 'nonchromosomal';
  }
    
  unless (exists $$chr{$slice_name}) {
    (my $chr_file = $out_file) =~ s/([^\.]+)$/$slice_name.$1/;
    open(my $chr_fh, '>', $chr_file) or $self->throw("Cannot open file $chr_file: $!");
    
    my $chr_serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($oa, $chr_fh);
    $chr_serializer->print_main_header([$slice]);
    
    $$chr{$slice_name}{'fh'} = $chr_fh;
    $$chr{$slice_name}{'file'} = $chr_file;
    $$chr{$slice_name}{'serializer'} = $chr_serializer;
  }
  
  return $$chr{$slice_name}{'serializer'};
}

sub Bio::EnsEMBL::Gene::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $self->source;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'Name'}            = $self->external_name;
  $summary{'biotype'}         = $self->biotype;
  $summary{'description'}     = $self->description;
  $summary{'version'}         = $self->version;
  
  return \%summary;
}

sub Bio::EnsEMBL::Transcript::summary_as_hash {
  my $self = shift;
  my %summary;
  
  my $parent_gene = $self->get_Gene();
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $parent_gene->source;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'Parent'}          = $parent_gene->stable_id;
  $summary{'biotype'}         = $self->biotype;
  $summary{'version'}         = $self->version;
  
  # Add xrefs
  my $xrefs = $self->get_all_xrefs();
  my (@db_xrefs, @go_xrefs);
  foreach my $xref (sort {$a->dbname cmp $b->dbname} @$xrefs) {
    my $dbname = $xref->dbname;
    if ($dbname eq 'GO') {
      push @go_xrefs, "$dbname:".$xref->display_id;
    } else {
      $dbname =~ s/^RefSeq.*/RefSeq/;
      $dbname =~ s/^Uniprot.*/UniProtKB/;
      $dbname =~ s/^protein_id.*/NCBI_GP/;
      push @db_xrefs,"$dbname:".$xref->display_id;
    }
  }
  $summary{'Dbxref'} = \@db_xrefs if scalar(@db_xrefs);
  $summary{'Ontology_term'} = \@go_xrefs if scalar(@go_xrefs);
  
  return \%summary;
}

sub Bio::EnsEMBL::Exon::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'constitutive'}    = $self->is_constitutive;
  
  return \%summary;
}

sub EnsEMBL::REST::EnsemblModel::CDS::summary_as_hash {
  my ($self) = @_;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'Parent'}          = $self->parent_id;
  $summary{'phase'}           = $self->phase;
  $summary{'source'}          = $self->source;
  
  return \%summary;
}

sub Bio::EnsEMBL::RepeatFeature::summary_as_hash {
  my $self = shift;
  my %summary;
  
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  
  my $rc = $self->repeat_consensus;
  $summary{'Name'}  = $rc->name;
  $summary{'type'}  = $rc->repeat_type;
  $summary{'class'} = $rc->repeat_class;
  if ($rc->repeat_consensus =~ /^[^N]\S*/) {
    $summary{'repeat_consensus'} = $rc->repeat_consensus;
  }
  
  return \%summary;
}

1;
