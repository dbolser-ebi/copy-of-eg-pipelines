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

package Bio::EnsEMBL::EGPipeline::FileDump::DumpFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type'    => 'core',
    'gene_dumps' => [],
    'skip_dumps' => [],
  };
}
    
sub write_output {
  my ($self) = @_;
  my $db_type    = $self->param_required('db_type');
  my $dump_types = $self->param_required('dump_types');
  my $gene_dumps = $self->param_required('gene_dumps');
  my $skip_dumps = $self->param_required('skip_dumps');
  
  my %gene_dumps = map { $_ => 1 } @$gene_dumps;
  my %skip_dumps = map { $_ => 1 } @$skip_dumps;
  
  my $dba = $self->get_DBAdaptor($db_type);
  my $has_chromosomes = $self->has_chromosomes($dba);
  my $has_genes      = $self->has_genes($dba);
  
  foreach my $flow (keys %$dump_types) {
    foreach my $dump_type (@{$$dump_types{$flow}}) {
      if (!exists $skip_dumps{$dump_type}) {
        if (!exists $gene_dumps{$dump_type} || $has_genes) {
          my %output_ids;
          $output_ids{'gene_centric'} = exists $gene_dumps{$dump_type} ? 1 : 0;
          
          if ($dump_type eq 'fasta_toplevel') {
            $self->fasta_toplevel($has_chromosomes, \%output_ids);
          }
          
          if ($dump_type eq 'fasta_seqlevel') {
            $self->fasta_seqlevel($has_chromosomes, \%output_ids);
          }
          
          if ($dump_type eq 'agp_assembly') {
            $self->agp_assembly($has_chromosomes, \%output_ids);
          }
          
          $self->dataflow_output_id(\%output_ids, $flow);
        }
      }
    }
  }
}

sub fasta_toplevel {
  my ($self, $has_chromosomes, $output_ids) = @_;
  
  if ($has_chromosomes) {
    $$output_ids{dump_level} = 'toplevel';
    $$output_ids{data_type}  = 'chromosomes';
  } else {
    $$output_ids{dump_level} = 'toplevel';
    $$output_ids{data_type}  = 'scaffolds';
  }
}

sub fasta_seqlevel {
  my ($self, $has_chromosomes, $output_ids) = @_;
  
  if ($has_chromosomes) {
    $$output_ids{dump_level} = 'scaffold';
    $$output_ids{data_type}  = 'scaffolds';
  } else {
    $$output_ids{dump_level} = 'contig';
    $$output_ids{data_type}  = 'contigs';
  }
}

sub agp_assembly {
  my ($self, $has_chromosomes, $output_ids) = @_;
  
  if ($has_chromosomes) {
    $$output_ids{seq_level} = 'scaffold';
    $$output_ids{data_type} = 'chromosome2scaffold';
  } else {
    $$output_ids{seq_level} = 'contig';
    $$output_ids{data_type} = 'scaffold2contig';
  }
}

1;
