
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

=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::Xref::BlastSearch

=head1 DESCRIPTION

Class for finding uniprot matches for a given sequence

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::Needle; 
use base Bio::EnsEMBL::EGPipeline::Xref::RestService;
use warnings;
use strict;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use LWP::UserAgent;
use XML::Simple;
use Data::Dumper;
use Carp;

=head1 CONSTRUCTOR
=head2 new


  Example    : $info = Bio::EnsEMBL::Utils::MetaData::GenomeInfo->new(...);
  Description: Creates a new info object
  Returntype : Bio::EnsEMBL::Utils::MetaData::GenomeInfo
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ( $proto, @args ) = @_;
  my $self = $proto->SUPER::new(@args);
  $self->{url} = 'http://www.ebi.ac.uk/Tools/services/rest/emboss_needle';
  return $self;
}


sub align {
  my ( $self, $sequences ) = @_;
  return $self->run_jobs($sequences,{format=>'pair'},'aln');
} 

sub parse_result {
    my ($self, $results, $id, $results_str) = @_;
    $results->{$id} = {};   
    for my $line (split '\n',$results_str) {
        $line.="\n";
        if($line =~ m/^#/) {
            if($line =~ m/# Identity:\s+(\d+)\/(\d+)\s+\(([0-9.]+)%\)/) {
                $results->{$id}{pos}.=$1;
                $results->{$id}{len}.=$2;
                $results->{$id}{id}.=$3;
            }
            elsif($line =~ m/1:\s+(\S+)/) {
                $results->{$id}{seq1}.=$1;
            }
            elsif($line =~ m/2:\s+(\S+)/) {
                $results->{$id}{seq2}.=$1;
            }
        } else {
            $results->{$id}{aln}.=$line;
        }
    }
    return;
}


1;
