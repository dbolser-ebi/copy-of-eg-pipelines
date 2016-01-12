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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteIniFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::ProcessSRA');

use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'bigwig'   => 0,
    'ini_type' => 'rnaseq_align',
  };
}

sub run {
  my ($self) = @_;
  
  my $work_dir = $self->param_required('work_directory');
  my $bigwig   = $self->param_required('bigwig');
  my $ini_type = $self->param_required('ini_type');
  my $mode     = $self->param_required('mode');
  my $merge_id = $self->param_required('merge_id');
  
  my ($name, $caption, $description);
  
  if ($mode eq 'file') {
    ($name, $caption, $description) = ($merge_id, '', '');
  } elsif ($mode eq 'study') {
    ($name, $caption, $description) = $self->study_ini($merge_id);
  } else {
    $self->throw("Unrecognised mode of operation, '$mode'");
  }
  
  my ($header, $data_file, $source_type);
  if ($bigwig) {
    $header = '[ENSEMBL_INTERNAL_BIGWIG_SOURCES]';
    $data_file = $self->param('bw_file');
    $source_type = 'bigWig';
  } else {
    $header = '[ENSEMBL_INTERNAL_BAM_SOURCES]';
    $data_file = $self->param('merged_bam_file');
    $source_type = 'bam';
  }
  
  $header .= "\n$merge_id = $ini_type\n";
  my $body = 
    "[$merge_id]\n".
    "source_name = $name\n".
    "caption     = $caption\n".
    "description = $description\n".
    "source_url  = $data_file\n".
    "source_type = $source_type\n".
    "display     = off\n";
  
  my $ini_file = catdir($work_dir, "/$merge_id.ini");
  open (my $fh, '>', $ini_file) or die "Failed to open file '$ini_file'";
	print $fh "$header\n$body\n";
  close($fh);
  
  $self->param('ini_file', $ini_file);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({'ini_file' => $self->param('ini_file')}, 1);
}

sub study_ini {
  my ($self, $merge_id) = @_;
  
  my $species       = $self->param_required('species');
  my $merge_level   = $self->param_required('merge_level');
  my $studies       = $self->param_required('study');
  my $study_species = $self->param_required('study_species');
  
  if (exists $$study_species{$species}) {
    push @$studies, $$study_species{$species};
  }
  
  my $merge_label = $self->merge_label($merge_level, $merge_id);
  my $study_ids   = $self->study_ids($studies, $merge_level, $merge_id);
  
  my @studies;
  foreach my $study_id (keys %$study_ids) {
    push @studies, $self->ena_link($study_id) . ': ' . $$study_ids{$study_id};
  }
  
  my $name = $merge_label;
  my $caption = "$merge_level $merge_id";
  my $description = ucfirst(lc($merge_level)) . ' ' . $self->ena_link($merge_id);
  if ($merge_level ne 'Study') {
    if (scalar @studies == 1) {
      $description .= ' (Study: ';
    } else {
      $description .= ' (Studies: ';
    }
    $description .= join('; ', @studies).')';
  }
  
  return ($name, $caption, $description);
}

sub ena_link {
  my ($self, $id) = @_;
  return "<a href='http://www.ebi.ac.uk/ena/data/view/$id' target='_blank'>$id</a>"
}

1;
