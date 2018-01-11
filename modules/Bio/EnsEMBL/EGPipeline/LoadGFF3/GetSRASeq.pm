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

package Bio::EnsEMBL::EGPipeline::LoadGFF3::GetSRASeq;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use LWP::Simple;
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    ena_url => 'http://www.ebi.ac.uk/ena/data/view/%s&display=fasta',
  };
}

sub run {
  my ($self) = @_;
  my $work_dir = $self->param_required('work_dir');
  my $ena_url  = $self->param_required('ena_url');
  my $sra_acc  = $self->param_required('sra_accession');
  
  path($work_dir)->mkpath unless -e $work_dir;
  my $seq_file = path($work_dir)->child("$sra_acc.fa");
  
  my $fasta = get sprintf($ena_url, $sra_acc);
  $fasta =~ s/^>ENA\|\w+\|/>/;
  
  $seq_file->spew($fasta);
  
  $self->param('seq_file', $seq_file->canonpath);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({ seq_file => $self->param('seq_file') }, 1);
}

1;
