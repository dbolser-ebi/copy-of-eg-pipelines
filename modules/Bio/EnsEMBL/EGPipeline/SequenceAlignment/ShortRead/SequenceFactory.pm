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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::SequenceFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'reformat_header' => 1,
    'trim_est'        => 1,
    'trimest_exe'     => '/nfs/panda/ensemblgenomes/external/EMBOSS/bin/trimest',
  };
}

sub write_output {
  my ($self) = @_;
  
  my $species         = $self->param_required('species');
  my $files           = $self->param('seq_file');
  my $file_species    = $self->param('seq_file_species');
  my $studies         = $self->param('study');
  my $study_species   = $self->param('study_species');
  my $taxids          = $self->param('taxids');
  my $merge_level     = $self->param_required('merge_level');
  my $data_type       = $self->param('data_type');
  my $reformat_header = $self->param('reformat_header');
  my $trim_est        = $self->param('trim_est');
  my $trimest_exe     = $self->param('trimest_exe');
  
  if (exists $$file_species{$species}) {
    push @$files, $$file_species{$species};
  }
  
  my @all;
  foreach my $fasta_file (@$files) {
    push @all, path($fasta_file)->basename;
  }
  
  foreach my $fasta_file (@$files) {
    if ($data_type eq 'est') {
      if ($trim_est) {
        $fasta_file = $self->trim_est($fasta_file);
      }
      if ($reformat_header) {
        $self->reformat_header($fasta_file);
      }
    }
    
    my $merge_id;
    if (lc($merge_level) eq 'file') {
      $merge_id = path($fasta_file)->basename;
    } else {
      $merge_id = join('_', @all);
    }
    
    my $dataflow_output = {
      'mode'       => 'file',
      'fasta_file' => $fasta_file,
      'merge_id'   => $merge_id,
    };
    $self->dataflow_output_id($dataflow_output, 3);
  }
  
  if (exists $$study_species{$species}) {
    push @$studies, $$study_species{$species};
  }
  
  my $taxid;
  if (exists $$taxids{$species}) {
    $taxid = $$taxids{$species};
  }
  
  foreach my $study_id (@$studies) {
    my $dataflow_output = {
      'mode'     => 'study',
      'study_id' => $study_id,
    };
    if (defined $taxid) {
      $dataflow_output->{taxid} = $taxid;
    }
    
    $self->dataflow_output_id($dataflow_output, 4);
  }
}

sub reformat_header {
  my ($self, $fasta_file) = @_;
  
  my $file = path($fasta_file);
  my $data = $file->slurp;
  $data =~ s/^>gi\|\d+\|gb\|([^\|]+)\|\S+/>$1/g;
  $file->spew($data);
}

sub trim_est {
  my ($self, $fasta_file) = @_;
  
  my $trimest_exe = $self->param('trimest_exe');
  (my $trimmed_file = $fasta_file) =~ s/(\.\w+)$/-trimmed$1/;
  
  my $cmd = "$trimest_exe -seq $fasta_file -out $trimmed_file";
  system($cmd) == 0 || $self->throw("Cannot execute $cmd");
  
  return $trimmed_file;
}

1;
