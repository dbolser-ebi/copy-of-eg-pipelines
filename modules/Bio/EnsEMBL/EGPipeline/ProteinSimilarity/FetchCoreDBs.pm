=head1 LICENSE

Copyright [1999-2014] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::ProteinSimilarity::FetchCoreDBs

=head1 DESCRIPTION

Dump proteomes from core databases into a single Fasta-format file.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::FetchCoreDBs;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::SeqIO;
use File::Path qw(make_path);
use Net::FTP;

sub param_defaults {
  return {
    
  };
}

sub fetch_input {
  my $self = shift @_;
  my $out_dir = $self->param_required('out_dir');
  
  if (!-e $out_dir) {
    warning "Output directory '$out_dir' does not exist. I shall create it.";
    make_path($out_dir) or throw "Failed to create output directory '$out_dir'";
  }
  
  my $fasta_file = "$out_dir/core_dbs.fa";
  $self->param('fasta_file', $fasta_file);
}

sub write_output {
  my $self = shift @_;
  
  my $output_id = {
    'fasta_file' => $self->param('fasta_file'),
  };
  $self->dataflow_output_id($output_id, 1);
}

sub run {
  my $self = shift @_;
  
  
}

1;
