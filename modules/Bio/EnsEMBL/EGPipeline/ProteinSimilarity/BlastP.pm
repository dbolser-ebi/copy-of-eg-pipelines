=head1 LICENSE

Copyright [1999-2014] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::ProteinSimilarity::BlastP

=head1 DESCRIPTION

.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::BlastP;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my $self = shift @_;
  
  my $proteome_dir = $self->param_required('proteome_dir');
  my $species      = $self->param_required('species');
  $self->param('out_file', "$proteome_dir/blast/$species.out");
  
}

sub write_output {
  my $self = shift @_;
  
}

sub run {
  my $self = shift @_;
  
  my $blast_db      = $self->param_required('blast_db');
  my $proteome_file = $self->param_required('proteome_file');
  my $logic_name    = $self->param_required('logic_name');
  my $program_file  = $self->param_required('program_file');
  my $parameters    = $self->param_required('parameters');
  my $out_file      = $self->param('out_file');
  
  my $out = `$program_file -db $blast_db -query $proteome_file -out $out_file $parameters 2>&1`;
  if ($out =~ /error/mi) {
    throw "Error when executing $program_file:\n$out";
  #} elsif () {
  #  throw "Error when executing $program_file:\n$out";
  }
  
}

1;
