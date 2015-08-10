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

package Bio::EnsEMBL::Hive::RunnableDB::ESTPipeline::LoadAlignments;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');
use Bio::EnsEMBL::EGPipeline::Common::SAMParser;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;

use Bio::EnsEMBL::Utils::Exception qw/throw/;

use Data::Dumper;

sub param_defaults {
  my ($self) = @_;
  return {};
}

sub fetch_input {
  my ($self) = @_;
  return;
}

sub run {
  my ($self) = @_;

  my $sam_file = $self->param('sam_file');
  my $calmd_file = $self->param('calmd_file');
  my $load_core  = $self->param("load_core");

  my $group = 'otherfeatures';
  if ($load_core) {
      $group = 'core';
  }

  my $otherfeatures_dbh = $self->core_dbh();
  if (!$load_core) {
      $otherfeatures_dbh = $self->otherfeatures_dbh();
  }
  
  # Get the analysis_id

  my $logic_name = $self->param('logic_name');

  warn("Fetching analysis object for logic_name, $logic_name\n");
  
  my $analysis_adaptor = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), $group, 'Analysis');
  my $analysis_obj = $analysis_adaptor->fetch_by_logic_name($logic_name);
  my $analysis_id = $analysis_obj->dbID();

  # Get the hash seq_region_id / name
  # Only required when loading the data through sql_inserts

  warn("Loading the seq_regions...\n");

  my $seq_href = {};
  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), $group, 'Slice');
  foreach my $slice_obj (@{$slice_adaptor->fetch_all('toplevel')}) {
      $seq_href->{$slice_obj->seq_region_name()} = $slice_obj->get_seq_region_id();
  }

  warn("Built sequence hash with " . keys (%$seq_href) . " sequences\n");

  # Parse the sam file

  my $sam_parser = Bio::EnsEMBL::EGPipeline::Common::SAMParser->new('sam_file' => $calmd_file);
  my $alignments_per_seqs_href = $sam_parser->parse_sam();
  
  warn("Loading the data\n");

  # $self->load_alignments_through_sql_inserts($alignments_per_seqs_href, $seq_href, $analysis_id, $otherfeatures_dbh);
  $self->load_alignments_through_api ($alignments_per_seqs_href, $analysis_obj, $group, $otherfeatures_dbh);

  # Check meta_coord table is accurate

  $self->load_meta_coord($otherfeatures_dbh);

  # Add dna_align_feature.level meta attribute if missing

  $self->add_meta_level_attribute($otherfeatures_dbh);

  warn("Loading done\n");
  
  return;
} ## end sub run

sub write_output {
  my ($self) = @_;
  return;
}


sub load_alignments_through_api {
    my ($self, $alignments_per_seqs_href, $analysis_obj, $group) = @_;
    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), $group, 'Slice');
    my $alignment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), $group, 'DnaAlignFeature');

    foreach my $seq_name (keys %$alignments_per_seqs_href) {
	
	my $slice = $slice_adaptor->fetch_by_region('toplevel',$seq_name);
	
	if (!defined $slice) {
	    throw("Failed to get a slice object for seq_region, $seq_name!\n");
	}
	
	foreach my $alignment_href (@{$alignments_per_seqs_href->{$seq_name}}) {
	    my $read_name  = $alignment_href->{read_name};
	    my $seq_start  = $alignment_href->{seq_start};
	    my $seq_end    = $alignment_href->{seq_end};
	    my $seq_strand = $alignment_href->{seq_strand};
	    my $cigar_line = $alignment_href->{cigar_line};
	    my $perc_id    = $alignment_href->{perc_id};
	    my $score      = $alignment_href->{score};
	    
	    my $hit_start  = $alignment_href->{hit_start};
	    my $hit_end    = $alignment_href->{hit_end};
	    my $hit_strand = $alignment_href->{hit_strand};
	    my $hcoverage  = $alignment_href->{hcoverage};
	    
	    my $evalue = undef;
	    my $external_data = undef;
	    
	    my $external_db_id = 700;
	    
	    # In Ensembl, in the cigar_line, 'N' is not supported, replace them by 'I' !
	    # Swap also 'D' and 'I'

	    $cigar_line =~ s/I/Z/g;
	    $cigar_line =~ s/[D|N]/I/g;
	    $cigar_line =~ s/Z/D/g;
	    
	    if ($hit_strand == -1) {
		$hit_strand = 1;
		$seq_strand = -1;
	    }
	    
	    my $fps_aref = [];
	    
            # Parse SAM cigar into ungapped blocks here
	    
	    my $cigar_line_parsing = $cigar_line;
	    
	    #print STDERR "alignment ref/query: $seq_start, $seq_end, $hit_start, $hit_end\n";
	    #print STDERR "cigar_line, $cigar_line\n";
	    
	    while ($cigar_line_parsing =~ s/(\d+[A-T])//) {
		my $block = $1;
		
		#print STDERR "block, $block\n";
		
		if ($block =~ /(\d+)M/) {
		    
		    my $block_length = $1;
		    
		    #print STDERR "block length, $block_length\n";

		    $seq_end = $seq_start + $block_length - 1;
		    if ($seq_strand == 1) {
			$hit_end = $hit_start + $block_length - 1;
		    }
		    else {
			$hit_start = $hit_end - $block_length + 1;
		    }

		    #print STDERR "Creating block ref start/end - query start/end, $seq_start, $seq_end, $hit_start, $hit_end\n";

		    my $fp = Bio::EnsEMBL::FeaturePair->new();
		    
		    $fp->start($seq_start);  
		    $fp->end($seq_end);  
		    $fp->strand($seq_strand);  
		    $fp->hseqname($read_name);  
		    $fp->hstart($hit_start);  
		    $fp->hend($hit_end);  
		    $fp->hstrand($hit_strand);
		    $fp->hcoverage($hcoverage);
		    $fp->slice($slice);
		    $fp->score($score);
		    $fp->percent_id($perc_id);
		    $fp->external_db_id($external_db_id);
		    $fp->analysis($analysis_obj);
		    
		    push (@$fps_aref, $fp);

		    $seq_start = $seq_end + 1;
		    if ($seq_strand == 1) {
			$hit_start = $hit_end + 1;
		    }
		    else {
			$hit_end = $hit_start - 1;
		    }
		    
		}
		elsif ($block =~ /(\d+)I/) {
		    my $block_length = $1;
		    
		    #print STDERR "block length, $block_length\n";
		    
		    # Insertion on the reference
		    
		    $seq_start += $block_length;
		    
		}
		elsif ($block =~ /(\d+)D/) {
		    my $block_length = $1;
		    
		    #print STDERR "block length, $block_length\n";
		    
		    # Deletion on the reference
		    
		    if ($seq_strand == 1) {
			$hit_start += $block_length;
		    }
		    else {
			$hit_end = $hit_start - $block_length - 1;
		    }

		}
		else {
		    warn("Can not parse this cigar block, $block (block of unknown type)\n");
		}

	    }
	    
	    #print STDERR "Dumping fps, " . Dumper ($fps_aref) . "\n";
	    #print STDERR "Nb fps: " . @$fps_aref . "\n";

	    my $alignment_obj = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => $fps_aref);  # will create the correct cigar string for you

	    # print STDERR "alignment_obj, " . Dumper ($alignment_obj) . "\n";

	    my @alignments = ($alignment_obj);
	    $alignment_adaptor->store(@alignments);    
	}
	
    }
    
}

sub load_alignments_through_sql_inserts {
    my ($self, $alignments_per_seqs_href, $seq_href, $analysis_id, $otherfeatures_dbh) = @_;
    
    #$otherfeatures_dbh->do("LOCK TABLES `dna_align_feature` WRITE;");  
    #$otherfeatures_dbh->do("/*!40000 ALTER TABLE `dna_align_feature` DISABLE KEYS */");

    foreach my $seq_name (keys %$alignments_per_seqs_href) {
	my $seq_id = $seq_href->{$seq_name};
        
	my $index = 1;
	my $insert_sql = "";
	
	foreach my $alignment_href (@{$alignments_per_seqs_href->{$seq_name}}) {
	    my $read_name  = $alignment_href->{read_name};
	    my $seq_start  = $alignment_href->{seq_start};
	    my $seq_end    = $alignment_href->{seq_end};
	    my $seq_strand = $alignment_href->{seq_strand};
	    my $cigar_line = $alignment_href->{cigar_line};
	    my $perc_id    = $alignment_href->{perc_id};
	    my $score      = $alignment_href->{score};
	    
	    my $hit_start  = $alignment_href->{hit_start};
	    my $hit_end    = $alignment_href->{hit_end};
	    my $hit_strand = $alignment_href->{hit_strand};
	    my $hcoverage  = $alignment_href->{hcoverage};
	    
	    my $evalue = undef;
	    my $external_data = undef;
	    
	    my $external_db_id = 700;
	    
	    # In Ensembl, in the cigar_line, 'N' is not supported, replace them by 'I' !
	    
	    $cigar_line =~ s/N/I/g;
	    
	    if ($index == 1) {
		
		$insert_sql .= "INSERT INTO `dna_align_feature` (seq_region_id,seq_region_start,seq_region_end,seq_region_strand,hit_start,hit_end,hit_strand,hit_name,analysis_id,score,perc_ident,cigar_line,external_db_id,hcoverage) VALUES ($seq_id,$seq_start,$seq_end,$seq_strand,$hit_start,$hit_end,$hit_strand,'$read_name',$analysis_id,$score,$perc_id,'$cigar_line',$external_db_id,$hcoverage)";
		
	    }
	    else {
		
		$insert_sql .= ",($seq_id,$seq_start,$seq_end,$seq_strand,$hit_start,$hit_end,$hit_strand,'$read_name',$analysis_id,$score,$perc_id,'$cigar_line',$external_db_id,$hcoverage)";
	    }
	    
	    $index++;
	}
	
	$otherfeatures_dbh->do($insert_sql);
	
    }
    
    $otherfeatures_dbh->do("/*!40000 ALTER TABLE `dna_align_feature` ENABLE KEYS */");
    $otherfeatures_dbh->do("UNLOCK TABLES");

}

1;

