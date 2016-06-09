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

package Bio::EnsEMBL::EGPipeline::FileDump::GenomicFeatureDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Bio::Seq::RichSeq;
use Bio::SeqFeature::Generic;
use Bio::Location::Simple;
use Bio::Location::Split;
use Bio::Species;

use Time::Piece;

my $MOL_TYPE = 'genomic DNA';
my @COMMENTS = ();

sub create_species_obj {
	my ($self, $mca) = @_;
  
  # Cannot use the meta table classifications, because unique key constraints
  # prevent the storage of the same name at different rank levels. This feels
  # like a more robust way of getting the data anyway.
  
  my $taxon_id       = $mca->get_taxonomy_id();
  my $dba            = $self->taxonomy_dba();
  my $nta            = $dba->get_TaxonomyNodeAdaptor();
  my $taxon_node     = $nta->fetch_by_taxon_id($taxon_id);
  my $ancestors      = $nta->fetch_ancestors($taxon_node);
  my @classification = ($taxon_node->name);
  foreach my $ancestor (@$ancestors) {
    if ($ancestor->name !~ /\s/ && $ancestor->name ne 'root') {
      push @classification, $ancestor->name;
    }
  }
  
  my $species_obj = Bio::Species->new(
    -ncbi_taxid     => $taxon_id,
    -classification => \@classification,
  );
  
  my $common_name = $mca->get_common_name();
  if(defined $common_name && $common_name ne $species_obj->binomial) {
    $species_obj->common_name($common_name);
	}
  
  my $strain = $mca->single_value_by_key('species.strain');
  if(defined $strain) {
    $species_obj->variant($strain);
	}
  
	return $species_obj;
}

sub create_seq_obj {
  my ($self, $slice, $species_obj) = @_;
  
  my $seq_obj = Bio::Seq::RichSeq->new(
    -molecule         => $MOL_TYPE,
    -accession_number => $slice->seq_region_name,
    -description      => $self->slice_description($species_obj, $slice),
    -seq              => $slice->seq,
    -length           => $slice->length,
    -is_circular      => $slice->is_circular,
    -species          => $species_obj,
  );
  
  my $t = localtime;
  $seq_obj->add_date(join('-', $t->mday, uc($t->month), $t->year));
  $seq_obj->add_keyword('.');
  
  return $seq_obj;
}

sub slice_description {
  my ($self, $species_obj, $slice) = @_;
  
  my @description = (
    $species_obj->binomial,
    $slice->coord_system->name,
    $slice->seq_region_name,
    $slice->coord_system->version || '',
    $slice->start . '..' . $slice->end,
    'from',
    $self->param_required('source'),
  );
  
  return join(' ', @description);
}

sub add_secondary_acc {
  my ($self, $slice, $seq_obj) = @_;
  
  my @synonyms = @{$slice->get_all_synonyms()};
  foreach my $synonym (@synonyms) {
    $seq_obj->add_secondary_accession($synonym->name);
  }
}

sub add_annotations {
  my ($self, $mca, $seq_obj) = @_;
  my $source       = $self->param_required('source');
  my $source_url   = $self->param_required('source_url');
  my $provider     = $mca->single_value_by_key('provider.name') || $source;
  my $provider_url = $mca->single_value_by_key('provider.url')  || $source_url;
  
  foreach my $comment (@COMMENTS) {
    $comment =~ s/<SOURCE>/$source/gm;
    $comment =~ s/<SOURCE_URL>/$source_url/gm;
    $comment =~ s/<PROVIDER>/$provider/gm;
    $comment =~ s/<PROVIDER_URL>/$provider_url/gm;
    
    my $comment_obj = Bio::Annotation::Comment->new();
    $comment_obj->text($comment);
    $seq_obj->add_Annotation('comment', $comment_obj);
  }
}

sub add_slice_feature {
  my ($self, $slice, $species_obj, $seq_obj) = @_;
  
  my $seq_feature = Bio::SeqFeature::Generic->new(
    -start        => $slice->start,
    -end          => $slice->end,
    -primary_tag  => 'source',
  );
  $seq_feature->add_tag_value('mol_type', $MOL_TYPE);
  $seq_feature->add_tag_value('organism', $species_obj->binomial);
  $seq_feature->add_tag_value('strain', $species_obj->variant);
  $seq_feature->add_tag_value('db_xref', 'taxon:'.$species_obj->ncbi_taxid);
  
  $seq_obj->add_SeqFeature($seq_feature);
}

sub add_features {
  my ($self, $feature_type, $adaptor, $slice, $seq_obj) = @_;
  
  # There's scope here for customised handling of other feature
  # types, but in lieu of that, everything that's not a gene
  # goes into a big generic bucket marked "misc_feature".
  
  if ($feature_type eq 'Gene') {
    $self->add_gene_features($feature_type, $adaptor, $slice, $seq_obj);
  } else {
    $self->add_generic_features($feature_type, $adaptor, $slice, $seq_obj);
  }
}

sub add_contig_features {
  my ($self, $slice, $seq_obj) = @_;
  
  foreach my $segment (@{$slice->project('seqlevel')}) {
    my ($start, $end, $ctg_slice) = @$segment;
    
    my $seq_feature = Bio::SeqFeature::Generic->new(
      -start        => $start,
      -end          => $end,
      -primary_tag  => 'misc_feature',
    );
    $seq_feature->add_tag_value('note', $ctg_slice->name);
    
    $seq_obj->add_SeqFeature($seq_feature);
  }
}

sub add_generic_features {
  my ($self, $feature_type, $adaptor, $slice, $seq_obj) = @_;
  
  my @features = @{$adaptor->fetch_all_by_Slice($slice)};
  
  foreach my $feature (@features) {
    my $seq_feature = Bio::SeqFeature::Generic->new(
      -start        => $feature->start,
      -end          => $feature->end,
      -strand       => $feature->strand,
      -primary_tag  => 'misc_feature',
    );
    if ($feature->display_id) {
      $seq_feature->add_tag_value('note', lc($feature_type)."_id=" . $feature->display_id);
    }
    
    $seq_obj->add_SeqFeature($seq_feature);
  }
}

sub add_gene_features {
  my ($self, $feature_type, $adaptor, $slice, $seq_obj) = @_;
  
  my @genes = @{$adaptor->fetch_all_by_Slice($slice)};
  
  foreach my $gene (@genes) {
    my $gene_xrefs = $gene->get_all_DBEntries;
    
    my $gene_feature = $self->create_gene_feature($gene);
    $self->add_gene_tags($gene, $gene_feature);
    $self->add_xref_tags($gene_xrefs, $gene_feature);
    
    $seq_obj->add_SeqFeature($gene_feature);
    
    my %exon_features;
    my @transcripts = sort {$a->stable_id cmp $b->stable_id} @{ $gene->get_all_Transcripts };
    
    foreach my $transcript (@transcripts) {
      my $transcript_xrefs = $transcript->get_all_DBEntries;
      
      my $transcript_feature = $self->create_transcript_feature($transcript);
      $self->add_transcript_tags($transcript, $transcript_feature);
      $self->add_xref_tags($transcript_xrefs, $transcript_feature);
      
      my $transcript_location = Bio::Location::Split->new();
      foreach my $exon (@{ $transcript->get_all_ExonTranscripts }) {
        my $exon_feature = $self->create_exon_feature($exon);
        $self->add_exon_tags($exon, $exon_feature);
        
        $transcript_location->add_sub_Location($exon_feature->location);
        
        $exon_features{$exon->stable_id} = $exon_feature;
      }
      $transcript_feature->location($transcript_location);
      
      $seq_obj->add_SeqFeature($transcript_feature);
      
      if (defined $transcript->translation) {
        my $cds_xrefs = $transcript->translation->get_all_DBEntries;
        
        my $cds_feature = $self->create_cds_feature($transcript);
        $self->add_cds_tags($transcript, $cds_feature);
        $self->add_xref_tags($cds_xrefs, $cds_feature);
        
        $seq_obj->add_SeqFeature($cds_feature);
      }
    }
    
    foreach my $exon_feature (sort {$a->start <=> $b->start} values %exon_features) {
      $seq_obj->add_SeqFeature($exon_feature);
    }
  }
}

sub create_gene_feature {
  my ($self, $gene) = @_;
  
  my $gene_feature = Bio::SeqFeature::Generic->new(
    -start       => $gene->start,
    -end         => $gene->end,
    -strand      => $gene->strand,
    -primary_tag => 'gene',
  );
  
  return $gene_feature;
}

sub add_gene_tags {
  my ($self, $gene, $gene_feature) = @_;
  
  $gene_feature->add_tag_value('gene', $gene->stable_id);
  if ($gene->display_xref && $gene->biotype eq 'protein_coding') {
    $gene_feature->add_tag_value('locus_tag', $gene->display_xref->display_id);
  }
  if ($gene->description) {
    $gene_feature->add_tag_value('note', $gene->description);
  }
}

sub create_transcript_feature {
  my ($self, $transcript) = @_;
  
  my $primary_tag;
  if ($transcript->biotype eq 'protein_coding') {
    $primary_tag = 'mRNA';
  } else {
    $primary_tag = 'misc_RNA';
  }
  
  my $transcript_feature = Bio::SeqFeature::Generic->new(
    -strand      => $transcript->strand,
    -primary_tag => $primary_tag,
  );
  
  return $transcript_feature;
}

sub add_transcript_tags {
  my ($self, $transcript, $transcript_feature) = @_;
  
  $transcript_feature->add_tag_value('gene', $transcript->get_Gene->stable_id);
  $transcript_feature->add_tag_value('note', 'transcript_id='.$transcript->stable_id);
  if ($transcript->biotype ne 'protein_coding') {
    $transcript_feature->add_tag_value('note', $transcript->biotype);
  }
}

sub create_exon_feature {
  my ($self, $exon) = @_;
  
  my $exon_feature = Bio::SeqFeature::Generic->new(
    -start       => $exon->start,
    -end         => $exon->end,
    -strand      => $exon->strand,
    -primary_tag => 'exon',
  );
  
  return $exon_feature;
}

sub add_exon_tags {
  my ($self, $exon, $exon_feature) = @_;
  
  $exon_feature->add_tag_value('note', 'exon_id='.$exon->stable_id);
}

sub create_cds_feature {
  my ($self, $transcript) = @_;
  
  my $cds_feature = Bio::SeqFeature::Generic->new(
    -strand      => $transcript->strand,
    -primary_tag => 'CDS',
  );
        
  my $cds_location = Bio::Location::Split->new();
        
  foreach my $cds (@{ $transcript->get_all_CDS }) {
    my $location = Bio::Location::Simple->new(
      -start  => $cds->start,
      -end    => $cds->end,
      -strand => $cds->strand,
    );
    $cds_location->add_sub_Location($location);
  }
  $cds_feature->location($cds_location);
  
  return $cds_feature;
}

sub add_cds_tags {
  my ($self, $transcript, $cds_feature) = @_;
  
  my $protein_name;
  my $desc = $transcript->get_Gene->description;
  if ($desc) {
    if ($desc !~ /(projected|hypothetical|putative|uncharacterized|predicted)/i) {
      ($protein_name = $desc) =~ s/\s+\[Source.*//;
    }
  }
  
  $cds_feature->add_tag_value('gene', $transcript->get_Gene->stable_id);
  $cds_feature->add_tag_value('note', 'transcript_id='.$transcript->stable_id);
  $cds_feature->add_tag_value('protein_id', $transcript->translation->stable_id);
  if ($protein_name) {
    $cds_feature->add_tag_value('product', $protein_name);
  }
  $cds_feature->add_tag_value('translation', $transcript->translate->seq);
}

sub add_xref_tags {
  my ($self, $xrefs, $feature) = @_;
  
  foreach my $xref (@$xrefs) {
    my ($dbname, $xref_id) = ($xref->dbname, $xref->display_id);
    $dbname  =~ s/^RefSeq.*/RefSeq/;
    $dbname  =~ s/^Uniprot.*/UniProtKB/;
    $dbname  =~ s/^protein_id.*/NCBI_GP/;
    $xref_id =~ s/^GO://;
    $feature->add_tag_value('db_xref', "$dbname:$xref_id");
  }
}

1;
