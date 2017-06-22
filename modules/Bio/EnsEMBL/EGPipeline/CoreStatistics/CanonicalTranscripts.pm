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

Bio::EnsEMBL::EGPipeline::CoreStatistics::CanonicalTranscripts

=head1 DESCRIPTION

Calculate canonical transcripts, using the Ensembl script:
ensembl/misc-scripts/canonical_transcripts/set_canonical_transcripts.pl.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::CoreStatistics::CanonicalTranscripts;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Production::Pipeline::Common::Base);

sub run {
  my ($self) = @_;
  my $script = $self->param('script');
  my $out_dir = $self->param('out_dir');
  my $dba = $self->get_DBAdaptor;
  my $dbc = $dba->dbc;
  my $species_id = $dba->species_id;
  
  my $options = '';
  if ($out_dir) {
    $self->throw("Output directory '$out_dir' does not exist.") unless -e $out_dir;
    $options = "--verbose &> $out_dir/".$dbc->dbname."_$species_id.can_transcript.out";
  }
  
  my $command = sprintf(
    "perl $script ".
      "--dbhost %s ".
      "--dbport %d ".
      "--dbuser %s ".
      "--dbpass %s ".
      "--dbname %s ".
      "--coord_system toplevel ".
      "--write ".
      "$options",
    $dbc->host, $dbc->port, $dbc->username, $dbc->password, $dbc->dbname
  );
  
 	unless (system($command) == 0) {
    $self->throw("Failed to execute script: '$command'.");
	}
  
  my $sth = $dbc->prepare(
    "SELECT COUNT(*) FROM gene WHERE canonical_transcript_id IS NULL"
  );
	$sth->execute();
	my $count = ( $sth->fetchrow_array() )[0];
	if ($count != 0) {
    $self->throw("Canonical transcripts not specified for $count genes.");
	}
}

1;
