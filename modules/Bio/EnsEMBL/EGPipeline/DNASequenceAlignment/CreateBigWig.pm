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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::CreateBigWig;

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
  my $ini_file      = $self->param_required('ini_file');
  
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
  system($bw_cmd) == 0 || $self->throw("Cannot execute $wig_cmd");

  if ($clean_up) {
    unlink $wig_file;
  }
  
  $self->update_inifile($ini_file, $bam_file, $bw_file);
  
  $self->param('bw_file', $bw_file);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({'bw_file' => $self->param('bw_file')}, 1);
}

sub update_inifile {
  my ($self, $ini_file, $bam_file, $bw_file) = @_;
  
  open (my $read_fh, '<', $ini_file) or die "Failed to open file '$ini_file'";
  my $ini_text;
  {
    local $/;
    $ini_text = <$read_fh>;
  }
  close($read_fh);
  
  $ini_text =~ s/ENSEMBL_INTERNAL_BAM_SOURCES/ENSEMBL_INTERNAL_BIGWIG_SOURCES/m;
  $ini_text =~ s/$bam_file/$bw_file/m;
  $ini_text =~ s/source_type\s*=\s*bam/source_type = bigWig/m;
  
  open (my $write_fh, '>', $ini_file) or die "Failed to open file '$ini_file'";
  print $write_fh $ini_text;
  close($write_fh);
}

1;
