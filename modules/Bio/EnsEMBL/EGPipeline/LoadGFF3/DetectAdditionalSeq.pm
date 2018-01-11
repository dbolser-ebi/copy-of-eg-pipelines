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

Bio::EnsEMBL::EGPipeline::LoadGFF3::DetectAdditionalSeq

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::LoadGFF3::DetectAdditionalSeq;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Path::Tiny qw(path);

sub run {
  my ($self) = @_;
  my $genbank_file = $self->param_required('genbank_file');
  
  my $accesssions = $self->find_accessions($genbank_file);
  
  $self->param('accessions', $accesssions);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $accession (@{ $self->param('accessions') }) {
    $self->dataflow_output_id({ sra_accession => $accession }, 2);
  }
}

sub find_accessions {
  my ($self, $genbank_file) = @_;
  
  my $genbank_path = path($genbank_file);
  my $genbank = $genbank_path->slurp;
  $genbank =~ s!\A//\s*!!m;
  $genbank =~ s!//\s*\Z!!m;
  my @genbank = split(m!//\s+!, $genbank);
  
  # The accessions used for additional sequence in the GenBank file
  # are stored in a predictable way, so can parse them out of the text.
  # Ideally, there would be a less fragile method, but this is the
  # only way I've found that works.
  my @accessions = ();
  
  foreach my $record (@genbank) {
    if ($record =~ /##RefSeq-Attributes-START##.*(assembly gap|frameshifts)/ms) {
      
      # Probably will only have mRNA, but check just in case.
      my ($mol_type) = $record =~ /\s+\/mol_type="([^"]+)"/m;
      next if $mol_type ne 'mRNA';
      
      my ($accessions) = $record =~ /\s+transcript\s+sequences*\s+\((.*?)\)/ms;
      if ($accessions) {
        push @accessions, split(/,\s*|\s*and\s*/, $accessions);
      }
    }
  }
  
  return \@accessions;
}

1;
