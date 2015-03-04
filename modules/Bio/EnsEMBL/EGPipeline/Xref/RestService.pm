
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

Bio::EnsEMBL::EGPipeline::Xref::RestService

=head1 DESCRIPTION

Class for finding uniprot matches for a given sequence

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::RestService;
use warnings;
use strict;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use List::MoreUtils qw(natatime);
use LWP::UserAgent;
use XML::Simple;
use Data::Dumper;
use Carp;

my $JOBN = 30;

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
  $self->{ua}     = LWP::UserAgent->new;
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

sub run_jobs {
  my ( $self, $inputs, $opts, $renderer ) = @_;
  my $results = {};

  my $submit_url = $self->{url}.'/run/';
  $opts->{email} ||= $ENV{USER}.'@ebi.ac.uk';
  $opts->{format} ||= 'xml';
  $renderer ||= 'xml';
  $self->logger()->debug("Processing total of ".scalar(keys %$inputs));
  my $it = natatime($JOBN,keys %$inputs);
  
  while(my @inputs_keys = $it->()) {

      $self->logger()->debug("Processing chunk of ".scalar(@inputs_keys));

      print Dumper($inputs);
      my %inputs_chunk = map { $_ => $inputs->{$_} } @inputs_keys;
      print Dumper(%inputs_chunk);
      my $job_ids = {};
      while ( my ( $id, $params ) = each %inputs_chunk ) {
          while(my($k,$v) = each %$opts) {
              $params->{$k} = $v;
          }
          # submit job using post
          $self->logger()->debug("Submitting job for $id to $submit_url");
          $job_ids->{$id} = $self->post( $submit_url, $params );
      }
      
      my $status_url = $self->{url}.'/status/';
      my $statuses = {};

      $self->logger()->debug("Awaiting completion");
      while () {
          while ( my ( $id, $job_id ) = each %$job_ids ) {
              if ( !defined $statuses->{$job_id} ) {
                  my $status = $self->get( $status_url . $job_id );
                  if ( $status ne 'RUNNING' ) {
                      $statuses->{$job_id} = $status;
                  }
              }
          }
          last if ( scalar( keys %$statuses ) == scalar( keys %$job_ids ) );
      }
      $self->logger()->debug("All jobs complete - parsing results");

      my $results_url = $self->{url}.'/result/';
      for my $id ( keys %inputs_chunk ) {
          my $job_id = $job_ids->{$id};
          my $status = $statuses->{$job_id};
          if ( $status eq 'FINISHED' ) {
              
              # get results
              my $results_str = $self->get( $results_url . $job_id . '/'. $renderer );
              
              $self->parse_result($results,$id,$results_str);
          }
          else {
              croak "Job $job_id completed with status $status";
          }       
      }
  }
  $self->logger()->debug("Returning ".scalar(keys %$results)." results");
  return $results;

} ## end sub run_blast

sub get {
  my ( $self, $url ) = @_;
  my $response = $self->{ua}->get($url);
  return $self->handle_response($response);
}

sub post {
  my ( $self, $url, $params ) = @_;
  my $response = $self->{ua}->post( $url, $params );
  return $self->handle_response($response);
}

sub handle_response {
  my ( $self, $response ) = @_;
  my $content = $response->content();
  if ( $response->is_error ) {
	my $error_message = '';
	# HTML response.
	if ( $content =~ m/<h1>([^<]+)<\/h1>/ ) {
	  $error_message = $1;
	}
	#  XML response.
	elsif ( $content =~ m/<description>([^<]+)<\/description>/ ) {
	  $error_message = $1;
	}
	croak 'Could not run job: ' . $response->code .
	  ' ' . $response->message . '  ' . $error_message;
  }
  return $content;

}

1;
