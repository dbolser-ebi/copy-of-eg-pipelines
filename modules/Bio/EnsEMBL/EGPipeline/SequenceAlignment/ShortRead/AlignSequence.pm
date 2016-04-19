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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::AlignSequence;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    'threads'   => 4,
    'read_type' => 'default',
    'clean_up'  => 1,
  };
}

sub fetch_input {
  my ($self) = @_;
  
  if (defined $self->param('escape_branch') and $self->input_job->retry_count >= $self->input_job->analysis->max_retry_count) {
    $self->dataflow_output_id($self->input_id, $self->param('escape_branch'));
    $self->input_job->autoflow(0);
    $self->complete_early("Failure probably due to memory limit, retrying with a higher limit.");
  }
  
  my $aligner_class = $self->param_required('aligner_class');
  my $aligner_dir   = $self->param_required('aligner_dir');
  my $samtools_dir  = $self->param_required('samtools_dir');
  my $threads       = $self->param_required('threads');
  my $read_type     = $self->param_required('read_type');
  
  my $max_intron_length;
  if ($aligner_class =~ /StarAligner/) {
    $max_intron_length = $self->max_intron_length();
  }
  
  eval "require $aligner_class";
  
  my $aligner_object = $aligner_class->new(
    -aligner_dir       => $aligner_dir,
    -samtools_dir      => $samtools_dir,
    -threads           => $threads,
    -read_type         => $read_type,
    -max_intron_length => $max_intron_length,
  );
  
  $self->param('aligner_object', $aligner_object);
}

sub run {
  my ($self) = @_;
  
  my $aligner     = $self->param_required('aligner_object');
  my $genome_file = $self->param_required('genome_file');
  my $seq_file_1  = $self->param_required('seq_file_1');
  my $seq_file_2  = $self->param('seq_file_2');
  my $clean_up    = $self->param_required('clean_up');
  
  my $sam_file = "$seq_file_1.sam";
  $aligner->align($genome_file, $sam_file, $seq_file_1, $seq_file_2);
  my $bam_file = $aligner->sam_to_bam($sam_file);
  unlink $sam_file if $clean_up;
  
  $self->param('bam_file', $bam_file);
}

sub write_output {
  my ($self) = @_;
  
  my $dataflow_output = {
    'bam_file' => $self->param('bam_file'),
  };
  
  $self->dataflow_output_id($dataflow_output, 1);
}

sub max_intron_length {
  my ($self) = @_;
  
  my $max_intron_length = 0;
  
  my $dba = $self->get_DBAdaptor('core');
  my $ta = $dba->get_adaptor('Transcript');
  my $transcripts = $ta->fetch_all_by_biotype('protein_coding');
  
  foreach my $transcript (@$transcripts) {
    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      if ($intron->length() > $max_intron_length) {
        $max_intron_length = $intron->length();
	    }
    }
    if ($max_intron_length > 25000) {
      $max_intron_length = 25000;
      last;
    }
  }
  
  return $max_intron_length;
}

1;
