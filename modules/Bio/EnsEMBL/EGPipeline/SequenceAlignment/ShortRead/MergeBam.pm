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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::MergeBam;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::EGPipeline::Common::Aligner;

use File::Path qw(make_path);
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
  
  my $results_dir = $self->param_required('results_dir');
  make_path($results_dir) unless -e $results_dir;
  
  my $dba = $self->core_dba;
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  $self->param('assembly', $assembly);
}

sub run {
  my ($self) = @_;
  
  my $results_dir = $self->param_required('results_dir');
  my $merges      = $self->param_required('merges');
  my $assembly    = $self->param('assembly');
  
  foreach my $merge (@$merges) {
    my $merge_id = $$merge{'merge_id'};
    my $merged_bam_file = catdir($results_dir, "$merge_id\_$assembly.bam");
    $$merge{'merged_bam_file'} = $merged_bam_file;
    
    $self->merge_bam($merge_id, $merged_bam_file);
  }
}

sub write_output {
  my ($self) = @_;
  
  foreach my $merge (@{$self->param_required('merges')}) {
    $self->dataflow_output_id($merge, 1);
  }
}

sub merge_bam {
  my ($self, $merge_id, $merged_bam_file) = @_;
  
  my $samtools_dir = $self->param_required('samtools_dir');
  my $vcf          = $self->param_required("vcf");
  my $use_csi      = $self->param_required("use_csi");
  my $clean_up     = $self->param_required('clean_up');
  
  my $sql = "SELECT bam_file FROM merge_bam WHERE merge_id = ?;";
  my $sth = $self->hive_dbh->prepare($sql);
  $sth->execute($merge_id);
  my @bam_files;
  while (my $results = $sth->fetch) {
    push @bam_files, $$results[0];
  }
  
  my $aligner = Bio::EnsEMBL::EGPipeline::Common::Aligner->new(
    -samtools_dir => $samtools_dir,
  );
  
  if (scalar(@bam_files) == 1) {
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
}

1;
