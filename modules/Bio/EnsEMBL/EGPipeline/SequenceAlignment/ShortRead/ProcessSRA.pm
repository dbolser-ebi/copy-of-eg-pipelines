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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::ProcessSRA;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor taxonomy_adaptor);
use File::Spec::Functions qw(catdir);

sub merge_id {
  my ($self, $merge_level, $study, $run) = @_;
  
  my $merge_id;
  
  my $experiment = $run->experiment;
  my $sample = $run->sample;
  
  if ($merge_level =~ /study/i) {
    $merge_id = $study->accession;
    
  } elsif ($merge_level =~ /taxon/i) {
    $merge_id = $sample->taxon ? $sample->taxon->taxon_id : 'unknown';
    
  } elsif ($merge_level =~ /sample/i) {
    $merge_id = $sample->accession;
    
  } elsif ($merge_level =~ /experiment/i) {
    $merge_id = $experiment->accession;
    
  } else {
    $self->throw("Unrecognised merge level, '$merge_level'");
    
  }
  
  return $merge_id;
}

sub merge_label {
  my ($self, $merge_level, $merge_id) = @_;
  
  my $merge_label;
  
  if ($merge_level =~ /study/i) {
    my $study_adaptor = get_adaptor('Study');
    my @study = @{$study_adaptor->get_by_accession($merge_id)};
    $merge_label = $study[0]->title;
    
  } elsif ($merge_level =~ /taxon/i) {
    if ($merge_id eq 'unknown') {
      $merge_label = 'unknown species';
    } else {
      my $taxonomy_adaptor = taxonomy_adaptor();
      my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
      my $node = $node_adaptor->fetch_by_taxon_id($merge_id);
      $merge_label = $node->name;
    }
    
  } elsif ($merge_level =~ /sample/i) {
    my $sample_adaptor = get_adaptor('Sample');
    my @sample = @{$sample_adaptor->get_by_accession($merge_id)};
    $merge_label = $sample[0]->title || $sample[0]->description;
    
  } elsif ($merge_level =~ /experiment/i) {
    my $experiment_adaptor = get_adaptor('Experiment');
    my @experiment = @{$experiment_adaptor->get_by_accession($merge_id)};
    $merge_label = $experiment[0]->title;
    
  } else {
    $self->throw("Unrecognised merge level, '$merge_level'");
    
  }
  
  return $merge_label;
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
    return ($seq_file_1, $seq_file_2, $sam_file);
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

sub study_ids {
  my ($self, $studies, $merge_level, $merge_id) = @_;
  
  my %study_ids;
  
  my $study_adaptor = get_adaptor('Study');
  foreach my $study_id (@$studies) {
    foreach my $study (@{$study_adaptor->get_by_accession($study_id)}) {
      foreach my $run (@{$study->runs()}) {
        my $run_merge_id = $self->merge_id($merge_level, $study, $run);
        if ($run_merge_id eq $merge_id) {
          $study_ids{$study_id} = $study->title;
        }
      }
    }
  }
  
  return \%study_ids;
}

1;
