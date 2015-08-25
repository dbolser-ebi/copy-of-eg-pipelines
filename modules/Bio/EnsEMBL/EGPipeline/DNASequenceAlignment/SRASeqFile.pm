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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::SRASeqFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);
use File::Spec::Functions qw(catdir);

sub run {
  my ($self) = @_;
  
  my $species     = $self->param_required('species');
  my $work_dir    = $self->param_required('work_directory');
  my $study_acc   = $self->param_required('study_acc');
  my $merge_level = lc($self->param_required('merge_level'));
  
  my @runs = ();
  my $study_adaptor = get_adaptor('Study');
  foreach my $study (@{$study_adaptor->get_by_accession($study_acc)}) {
    foreach my $run (@{$study->runs()}) {
      my $experiment = $run->experiment;
      my $sample = $run->sample;
      
      my ($merge_id, $merge_label);
      if ($merge_level eq 'study') {
        $merge_id    = $study->accession;
        $merge_label = $study->title;
        
      } elsif ($merge_level eq 'taxon') {
        $merge_id    = $sample->taxon ? $sample->taxon->taxon_id : 'unknown';
        $merge_label = $sample->taxon ? $sample->taxon->name : 'unknown species';
        
      } elsif ($merge_level eq 'sample') {
        $merge_id    = $sample->accession;
        $merge_label = $sample->title || $sample->description;
        
      } elsif ($merge_level eq 'experiment') {
        $merge_id    = $experiment->accession;
        $merge_label = $experiment->title;
        
      } else {
        $self->throw("Unrecognised merge level, '$merge_level'");
        
      }
      
      my ($seq_file_1, $seq_file_2, $sam_file) = $self->retrieve_files($work_dir, $run);
      
      push @runs,
        {
          seq_file_1  => $seq_file_1,
          seq_file_2  => $seq_file_2,
          sam_file    => $sam_file,
          merge_id    => $merge_id,
          merge_label => $merge_label,
          study_id    => $study_acc,
        };
    }
  }
  
  $self->param('runs', \@runs);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $dataflow_output (@{$self->param('runs')}) {
    $self->dataflow_output_id($dataflow_output, 2);
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
  
  if ($paired) {
    if (scalar(@files) != 2) {
	    $self->throw("A paired read experiment can only have 2 files");
    }
    
    $seq_file_1 = catdir($work_dir, "$run_acc\_all_1.fastq");
    $seq_file_2 = catdir($work_dir, "$run_acc\_all_2.fastq");
    
    if (-e $seq_file_2) {
      $self->warning("File '$seq_file_2' exists, and will be deleted");
      unlink $seq_file_2;
    }
  } else {
    $seq_file_1 = catdir($work_dir, "$run_acc\_all.fastq");
  }
    
  if (-e $seq_file_1) {
    $self->warning("File '$seq_file_1' exists, and will be deleted");
    unlink $seq_file_1;
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
      unlink $fq;
    } else {
      $self->throw("Cannot process file '$file_name'");
    }
  }
  
  return ($seq_file_1, $seq_file_2, $sam_file);
}

1;
