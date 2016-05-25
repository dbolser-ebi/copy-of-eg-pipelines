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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SRASeqFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_dir');
  my $run_id   = $self->param_required('run_id');
  
  make_path($work_dir) unless -e $work_dir;
  
  my $run_adaptor = get_adaptor('Run');
  my @runs = @{$run_adaptor->get_by_accession($run_id)};
  
  my @output = ();
  
  foreach my $run (@runs) {
    my ($seq_file_1, $seq_file_2, $sam_file) = $self->retrieve_files($work_dir, $run);
    
    push @output,
      {
        'seq_file_1' => $seq_file_1,
        'seq_file_2' => $seq_file_2,
        'sam_file'   => $sam_file,
      };
  }
  
  $self->param('output', \@output);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $output (@{$self->param('output')}) {
    $self->dataflow_output_id($output, 2);
  }
}

sub retrieve_files {
  my ($self, $work_dir, $run) = @_;
  
  my $run_acc = $run->accession;
  my @files = @{$run->files()};
  
  my $experiment = $run->experiment();
  my $paired = defined $experiment->design()->{LIBRARY_DESCRIPTOR}{LIBRARY_LAYOUT}{PAIRED};
  
  my ($seq_file_1, $seq_file_2);
  my $sam_file = catdir($work_dir, "$run_acc.sam");
  
  # Prepare files names
  if ($paired) {
    if (scalar(@files) != 2) {
	    $self->throw("A paired read experiment can only have 2 files");
    }
    
    $seq_file_1 = catdir($work_dir, "$run_acc\_all_1.fastq");
    $seq_file_2 = catdir($work_dir, "$run_acc\_all_2.fastq");
  } else {
    $seq_file_1 = catdir($work_dir, "$run_acc\_all.fastq");
  }
  
  # Retrieve each file
  FILE: for my $file (@files) {
    my $file_name = $file->file_name();
    
    # Decide the file name to use
    if (index($file_name, '.fastq') != -1) {
      my $seq_file = $seq_file_1;
      
      # Choose a name for a pair
      if ($paired) {
        if (index($file_name, '_1.fastq') != -1) {
          $seq_file = $seq_file_1;
        } elsif (index($file_name, '_2.fastq') != -1) {
          $seq_file = $seq_file_2;
        } else {
          $self->throw("Cannot process paired end file '$file_name'");
        }
      }
      
      # Reuse files if possible
      my $fastq = $seq_file;
      my $fastq_gz = $file->file_name;
      
      if (-s $fastq_gz) {
        if (-s $fastq) {
          # Interrupted unzip
          unlink $fastq;
        } else {
          # Interrupted download
          unlink $fastq_gz;
          $file->retrieve($work_dir);
        }
        $self->_unzip($fastq_gz, $fastq);
      } else {
        if (-s $fastq) {
          # Finished: can reuse
          next FILE;
        } else {
          # All new
          $file->retrieve($work_dir);
          $self->_unzip($fastq_gz, $fastq);
        }
      }
    } else {
      $self->throw("Cannot process file '$file_name'");
    }
  }
  
  return ($seq_file_1, $seq_file_2, $sam_file);
}

sub _unzip {
  my $self = shift;
  my ($fastq_gz, $fastq) = @_;

  my $cmd = "zcat $fastq_gz >> $fastq";
  system($cmd) == 0 || $self->throw("Cannot execute $cmd");
  unlink $fastq_gz;
  
  return $fastq;
}

1;
