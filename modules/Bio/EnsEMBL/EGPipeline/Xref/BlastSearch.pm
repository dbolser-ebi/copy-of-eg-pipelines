
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

package Bio::EnsEMBL::EGPipeline::Xref::BlastSearch; 
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

#my $srv = "ncbiblast";
my $srv = "wublast";

sub new {
  my ( $proto, @args ) = @_;
  my $self = $proto->SUPER::new(@args);
  $self->{url} = 'http://www.ebi.ac.uk/Tools/services/rest/'.$srv;
  return $self;
}

=head1 METHODS
=head2 logger
  Arg        : (optional) logger to set
  Description: Get logger
  Returntype : logger object reference
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut

sub logger {
  my ($self) = @_;
  if ( !defined $self->{logger} ) {
	$self->{logger} = get_logger();
  }
  return $self->{logger};
}

sub search {
  my ( $self, $seq, $collection, $stype, $program, $opts ) = @_;
  if(!defined $opts) {
      $opts = { 
          exp      => '1e-4'
      }
  }
  if ( ref $seq eq 'HASH' ) {
	my %params = map {
	  $_ => { 
			  sequence => $seq->{$_},
			  database => [$collection],
			  program  => $program,
			  stype    => $stype }
	} keys %$seq;
	return $self->run_blast( \%params, $opts );

  }
  else {
      return
	  $self->run_blast(
						{ 1 => { email    => 'ensgen@ebi.ac.uk',
								 sequence => $seq,
								 database => [$collection],
								 program  => $program,
								 stype    => $stype } } )->{1};
  }
} ## end sub search

sub run_blast {
  my ( $self, $inputs, $opts ) = @_;
  my $results = $self->run_jobs($inputs,$opts);
  while (my ($id,$res) = each %$results) {
      while( my ($hit,$hres) = each %$res) {
          $hres->{query} = $id;
          $hres->{expectation} = $hres->{alignments}{alignment}->[0]->{expectation}; 
      }
  }

  return $results;

} ## end sub run_blast

sub parse_result {
    my ($self, $results, $id, $results_str) = @_;
    # parse results as XML
    my $pres = XMLin($results_str)->{SequenceSimilaritySearchResult};
# service treats single and multiple hits differently in format, so need to check
    if ( $pres->{hits}{total} == 1 ) {
        $results->{$id} =
        { $pres->{hits}->{hit}->{id} => $pres->{hits}->{hit} };
    }
    else {
        $results->{$id} = $pres->{hits}->{hit};
    }    
    return;
}


1;
