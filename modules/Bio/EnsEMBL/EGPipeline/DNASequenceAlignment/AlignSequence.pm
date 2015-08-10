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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::AlignSequence;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);
use File::Spec::Functions qw(catdir);

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
  
  my $mode        = $self->param_required('mode');
  my $work_dir    = $self->param_required('work_directory');
  my $genome_file = $self->param_required('genome_file');
  my $aligner     = $self->param_required('aligner_object');
  my $clean_up    = $self->param_required('clean_up');
  
  my ($seq_file_1, $seq_file_2, $sam_file);
  
  if ($mode eq 'file') {
    $seq_file_1 = $self->param_required('split_file');
    $sam_file = "$seq_file_1.sam";
    
  } elsif ($mode eq 'study') {
    my $run_acc  = $self->param_required('run_acc');
    ($seq_file_1, $seq_file_2) = $self->retrieve_files($work_dir, $run_acc);
    $sam_file = catdir($work_dir, "$run_acc.sam");
    
  } else {
    $self->throw("Unrecognised mode of operation, '$mode'");
    
  }
  
  $aligner->align($genome_file, $sam_file, $seq_file_1, $seq_file_2);
  my $bam_file = $aligner->sam_to_bam($sam_file);
  unlink $sam_file if $clean_up;
  
  $self->param('bam_file', $bam_file);
}


sub write_output {
  my ($self) = @_;
  
  my $dataflow_output = {
    'bam_file' => $self->param('bam_file'),
    'merge_id' => $self->param('merge_id'),
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

sub retrieve_files {
  my ($self, $work_dir, $run_acc) = @_;
  
  my $run_adaptor = get_adaptor('Run');
  my ($run) = @{$run_adaptor->get_by_accession($run_acc)};
  my @files = @{$run->files()};
  
  my $experiment = $run->experiment();
  my $paired = defined $experiment->design()->{LIBRARY_DESCRIPTOR}{LIBRARY_LAYOUT}{PAIRED};
  
  my ($seq_file_1, $seq_file_2);
  
  if ($paired) {
    if (scalar(@files) != 2) {
	    $self->throw("A paired read experiment can only have 2 files");
    }
    
    $seq_file_1 = catdir($work_dir, "$run_acc\_all_1.fastq");
    $seq_file_2 = catdir($work_dir, "$run_acc\_all_2.fastq");
    
    if (-e $seq_file_2) {
      $self->warning("File '$seq_file_2' exists, and will be deleted");
      #unlink $seq_file_2;
    }
  } else {
    $seq_file_1 = catdir($work_dir, "$run_acc\_all.fastq");
  }
    
  if (-e $seq_file_1) {
    $self->warning("File '$seq_file_1' exists, and will be deleted");
    #unlink $seq_file_1;
  }
  
  for my $file (@files) {
    my $file_name = $file->file_name();
    
    if (index($file_name, '.fastq') != -1) {
      my $fq = $file->retrieve($work_dir);
      my $seq_file = $seq_file_1;
      
      if ($paired) {
        if (index($file_name, '_1.fastq') != -1) {
          $seq_file = $seq_file_1;
        } elsif (index($file_name, '_2.fastq') != -1) {
          $seq_file = $seq_file_2;
        } else {
          $self->throw("Cannot process paired end file '$file_name'");
        }
      }
      
      my $cmd = "zcat $fq >> $seq_file";
      system($cmd) == 0 || $self->throw("Cannot execute $cmd");
      #unlink $fq;
    } else {
      $self->throw("Cannot process file '$file_name'");
    }
  }
  
  return ($seq_file_1, $seq_file_2);
}

1;
