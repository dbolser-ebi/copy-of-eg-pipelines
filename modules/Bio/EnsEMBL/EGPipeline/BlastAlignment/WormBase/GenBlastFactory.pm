=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::GenBlastFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run{
  my ($self)=@_;

  my $dba = $self->get_DBAdaptor('core');
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name('genblast');

  my @files;
  $self->param('files', \@files);

  # as it will barf the method call below otherwise
  return undef unless $analysis; 

  my $genblast_files = $analysis->db_file;

  # species who don't need genblasting will return NULL
  return undef unless -e $genblast_files;

  $genblast_files=~s/genome\/genome\.fa/split_wormpep/;
  @files = glob("$genblast_files/*");
}

sub write_output{
  my ($self)=@_;
  foreach my $file(@{$self->param('files')}){
    $self->dataflow_output_id({query => $file}, 2)
  }
}

1;
