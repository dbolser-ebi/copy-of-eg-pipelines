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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::MergeBam;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::EGPipeline::Common::Aligner;
use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'samtools_dir' => '/nfs/panda/ensemblgenomes/external/samtools',
    'vcf'          => 0,
    'use_csi'      => 0,
    'clean_up'     => 1,
  };
}

sub fetch_input {
	my ($self) = @_;
  
	my $study_ids = $self->param_required('study_id');
  $self->warning('STUDY IDS '.join(',', @$study_ids));
  
	my $merge = $self->param_required('merge');
  # Invert the hash and gather old keys into a new value list.
  my %merge_bam;
  foreach my $bam_file (keys %$merge) {
    my $merge_id = $$merge{$bam_file};
    push @{$merge_bam{$merge_id}}, $bam_file;
  }
  
  $self->param('merge_bam', \%merge_bam);
}
  
sub run {
	my ($self) = @_;
  
  my $work_dir     = $self->param_required('work_directory');
  my $samtools_dir = $self->param_required('samtools_dir');
	my $vcf          = $self->param_required("vcf");
	my $use_csi      = $self->param_required("use_csi");
  my $clean_up     = $self->param_required('clean_up');
	my $merge_bam    = $self->param_required('merge_bam');
  
  my $aligner = Bio::EnsEMBL::EGPipeline::Common::Aligner->new(
    -samtools_dir => $samtools_dir,
  );
  
  my @merged_bam_files;
  
  foreach my $merge_id (keys %$merge_bam) {
    my $merged_bam_file = catdir($work_dir, "$merge_id.bam");
    my @bam_files = @{$$merge_bam{$merge_id}};
    my $size = scalar @bam_files;
    
    if ($size == 1) {
      my ($bam_file) = $bam_files[0];
      rename $bam_file, $merged_bam_file;
    } else {
      $aligner->merge_bam(\@bam_files, $merged_bam_file);
      if ($clean_up) {
        map { unlink $_ } @bam_files;
			}
		}
    
    my $sorted_bam = $aligner->sort_bam($merged_bam_file);
    rename $sorted_bam, $merged_bam_file;
    $aligner->index_bam($merged_bam_file, $use_csi);
    if ($vcf) {
      $aligner->generate_vcf($merged_bam_file);
    }
    
    push @merged_bam_files, $merged_bam_file;
	}
  
  $self->param('merged_bam_files', \@merged_bam_files);
}

sub write_output {
  my ($self) = @_;
  
	my $bigwig = $self->param_required('bigwig');
  my $branch = $bigwig ? 4 : 3;
  
  foreach my $merged_bam_file (@{$self->param('merged_bam_files')}) {
    $self->dataflow_output_id({'merged_bam_file' => $merged_bam_file}, $branch);
  }
}

# Generate filename with assembly and study if poss
#  my $dba = $self->core_dba;
#  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');

1;
