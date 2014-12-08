
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

Bio::EnsEMBL::EGPipeline::Xref::UniProtSequenceSearch

=head1 DESCRIPTION

Class for finding uniprot matches for a given sequence

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::UniProtSequenceSearch;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::IdentityXref;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Data::Dumper;

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
  my $class = ref($proto) || $proto;
  my $self = bless( {}, $class );
  $self->{logger} = get_logger();
  ( $self->{uniprot_dba} ) = rearrange( ['UNIPROT_DBA'], @args );
  $self->{blast} = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
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

sub get_taxid {
  my ( $self, $id ) = @_;
  return
	$self->{uniprot_dba}->dbc()->sql_helper()->execute_single_result(
					  -SQL => 'select tax_id from dbentry where name=?',
					  -PARAMS => [$id] );
}

sub search {
  my ( $self, $seq, $collection, $taxids ) = @_;
  # blast results - hash of results by UniProt ID
  my $results =
	$self->{blast}->search( $seq, $collection, 'protein', 'blastp' );
  if ( ref $seq eq 'HASH' ) {
	my $output = {};
	while ( my ( $id, $res ) = each %$results ) {
	  $output->{$id} = $self->parse( $res, $taxids );
	}
	return $output;
  }
  else {
	return $self->parse( $results, $taxids );
  }
}

sub parse {
  my ( $self, $results, $taxids ) = @_;
  if ( defined $taxids && scalar(@$taxids) > 0 ) {
	my %taxids = map { $_ => 1 } @$taxids;
	# filter using UniProt to find taxid
	for my $id ( keys %$results ) {
	  my $taxid = $self->get_taxid($id);
	  if ( !defined $taxids{$taxid} ) {
		delete $results->{$id};
	  }
	  else {
		$results->{$id}{taxid} = $taxid;
	  }
	}
  }
  # return IdentityXref objects
  my @results = grep {$_->ensembl_identity()>=90} map {
	my $result = $results->{$_};
	my $ali    = $result->{alignments}->{alignment};
	Bio::EnsEMBL::IdentityXref->new(
	   -XREF_IDENTITY    => $ali->{identity},
	   -ENSEMBL_IDENTITY => $ali->{identity},
	   -SCORE            => $ali->{score},
	   -EVALUE           => $ali->{expectation},
	   -CIGAR_LINE       => $ali->{pattern},
	   -XREF_START       => $ali->{matchSeq}->{start},
	   -XREF_END         => $ali->{matchSeq}->{end},
	   -ENSEMBL_START    => $ali->{querySeq}->{start},
	   -ENSEMBL_END      => $ali->{querySeq}->{end},
	   -PRIMARY_ID       => $result->{ac},
	   -DISPLAY_ID       => $result->{ac},
	   -DESCRIPTION      => $result->{description},
	   -DBNAME => ( $result->{database} eq 'TR' ) ? 'Uniprot/SPTREMBL' :
		 'Uniprot/SWISSPROT' );
  } keys %$results;
  return \@results;
} ## end sub parse

1;
