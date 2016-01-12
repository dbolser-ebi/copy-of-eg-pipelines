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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::CreateBigWig;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    'bedtools_dir'  => '/nfs/panda/ensemblgenomes/external/bedtools/bin',
    'ucscutils_dir' => '/nfs/panda/ensemblgenomes/external/ucscutils',
    'clean_up'      => 1,
  };
}

sub run {
	my ($self) = @_;
  
  my $bedtools_dir  = $self->param_required('bedtools_dir');
  my $ucscutils_dir = $self->param_required('ucscutils_dir');
  my $clean_up      = $self->param_required('clean_up');
  my $length_file   = $self->param_required('length_file');
  my $bam_file      = $self->param_required('merged_bam_file');
  
  my $wig_file = "$bam_file.wig";
  (my $bw_file = $bam_file) =~ s/\.bam$/\.bw/;
  if ($bw_file eq $bam_file) {
    $bw_file = "$bam_file.bw";
  }
  
  my $wig_cmd =
    "$bedtools_dir/bedtools genomecov ".
    " -g $length_file ".
    " -ibam $bam_file ".
    " -bg ".
    " -split ".
    " > $wig_file ";
  system($wig_cmd) == 0 || $self->throw("Cannot execute $wig_cmd");
  
  my $bw_cmd =
    "$ucscutils_dir/wigToBigWig ".
    " $wig_file ".
    " $length_file ".
    " $bw_file ";
  system($bw_cmd) == 0 || $self->throw("Cannot execute $bw_cmd");

  if ($clean_up) {
    unlink $wig_file;
  }
  
  $self->param('bw_file', $bw_file);
}

sub write_output {
  my ($self) = @_;
  
  my $dataflow_output = {
    'bw_file' => $self->param('bw_file'),
  };
  
  $self->dataflow_output_id($dataflow_output, 1);
}

1;
