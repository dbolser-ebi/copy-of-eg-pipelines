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

package Bio::EnsEMBL::EGPipeline::FileDump::MAFDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'compara'           => 'multi',
    'stats_base_url'    => 'https://www.vectorbase.org',
    'maf_source'        => 'VectorBase',
    'file_per_chr'      => 1,
    'file_per_scaffold' => 0,
  };
}

sub fetch_input {
  my ($self) = @_;
  
  if (defined $self->param('escape_branch') and 
      $self->input_job->retry_count >= $self->input_job->analysis->max_retry_count) 
  {
    $self->dataflow_output_id($self->input_id, $self->param('escape_branch'));
    $self->input_job->autoflow(0);
    $self->complete_early("Failure probably due to memory limit, retrying with a higher limit.");
  }
}

sub run {
  my ($self) = @_;
  my $compara     = $self->param_required('compara');
  my $mlss_id     = $self->param_required('mlss_id');
  my $aln_dir     = $self->param_required('aln_dir');
  my $ref_species = $self->param_required('ref_species');
  
  path($aln_dir)->mkpath;
  
  my $dba   = Bio::EnsEMBL::Registry->get_DBAdaptor($compara, 'compara');
  my $mlssa = $dba->get_adaptor("MethodLinkSpeciesSet");
  my $gaba  = $dba->get_adaptor("GenomicAlignBlock");
  
  my $mlss = $mlssa->fetch_by_dbID($mlss_id);
  my $gabs = $gaba->fetch_all_by_MethodLinkSpeciesSet($mlss);
  
  my %maf;
  foreach my $gab (@$gabs) {
    $self->maf_format($ref_species, $gab, \%maf);
  }
  
  my @out_files;
  
  if ($self->multiple_maf_files($ref_species)) {
    foreach my $scaffold (sort keys %maf) {
      my $out_file = catdir($aln_dir, "$scaffold.maf");
      open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
      
      print $out_fh $self->maf_header($mlss_id);
      foreach my $start_pos (sort {$a <=> $b} keys %{$maf{$scaffold}}) {
        print $out_fh $maf{$scaffold}{$start_pos};
      }
      print $out_fh "\n";
      close($out_fh);
      
      push @out_files, $out_file;
    }
  } else {
    my $out_file = catdir($aln_dir, 'all_scaffolds.maf');
    open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
    
    print $out_fh $self->maf_header($mlss_id);
    foreach my $scaffold (sort keys %maf) {
      foreach my $start_pos (sort {$a <=> $b} keys %{$maf{$scaffold}}) {
        print $out_fh $maf{$scaffold}{$start_pos};
      }
    }
    print $out_fh "\n";
    close($out_fh);
      
    push @out_files, $out_file;
  }
  
  $self->param('out_files', \@out_files);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $out_file (@{$self->param('out_files')}) {
    $self->dataflow_output_id({out_file => $out_file}, 2);
  }
}

sub maf_header {
  my ($self, $mlss_id) = @_;
  my $stats_base_url = $self->param_required('stats_base_url');
  my $maf_source     = $self->param_required('maf_source');
  my $release_date   = $self->param_required('release_date');
  
  my $url = $stats_base_url.'/mlss.html?mlss='.$mlss_id;
  
	my $header;
	$header  = "##maf version=1\n";
	$header .= "# Alignment method and statistics: $url \n";
	$header .= "# Exported from $maf_source, $release_date \n";
  
  return $header;
}

sub maf_format {
  my ($self, $ref_species, $gab, $maf) = @_;
  
  my $maf_block = "\na score=" . $gab->score . "\n";
  
  my $gas = $gab->get_all_GenomicAligns();
  
  my @sorted_gas;
  my ($ref_dnafrag_name, $ref_dnafrag_start);
  
  foreach my $ga (@$gas) {
    if ($ga->genome_db->name eq $ref_species) {
      unshift @sorted_gas, $ga;
      $ref_dnafrag_name  = $ga->dnafrag->name;
      $ref_dnafrag_start = $ga->dnafrag_start;
    } else {
      push @sorted_gas, $ga;
    }
  }
  
  foreach my $ga (@sorted_gas) {
    # MAF-format has zero-based starts, and if on the -ve strand
    # need to calculate position on revcomp sequence.
    my $start;
    if ($ga->dnafrag_strand == 1) {
      $start = $ga->dnafrag_start - 1;
    } else {
      $start = $ga->dnafrag->length - $ga->dnafrag_end;
    }
    
    $ga->original_sequence($ga->get_Slice->get_repeatmasked_seq(undef,1)->seq);
    
    my @line = (
      $ga->genome_db->assembly . '.' . $ga->dnafrag->name,
      $start,
      $ga->dnafrag_end - $ga->dnafrag_start + 1,
      $ga->dnafrag_strand eq 1 ? '+' : '-',
      $ga->dnafrag->length,
      $ga->aligned_sequence,
    );
    $maf_block .= 's '.join("\t", @line)."\n",
  }
  
  $$maf{$ref_dnafrag_name}{$ref_dnafrag_start} = $maf_block;
}

sub multiple_maf_files {
  my ($self, $ref_species) = @_;
  my $file_per_chr      = $self->param_required('file_per_chr');
  my $file_per_scaffold = $self->param_required('file_per_scaffold');
  
  my $multiple_maf_files = 0;
  if ($file_per_scaffold) {
    $multiple_maf_files = 1;
  } else {
    if ($file_per_chr) {
      my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($ref_species, 'core');
      if ($self->has_chromosomes($core_dba)) {
        $multiple_maf_files = 1;
      }
      $core_dba->dbc->disconnect_if_idle();
    }
  }
  
  return $multiple_maf_files;
}

1;
