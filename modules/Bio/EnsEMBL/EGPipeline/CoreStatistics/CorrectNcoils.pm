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

Bio::EnsEMBL::EGPipeline::CoreStatistics::CorrectNcoils

=head1 DESCRIPTION

The ncoils program sometimes creates features that end one base beyond the
end of the translation. This seems to be a bug in the ncoils binary, rather
than anything we're doing wrong. Although a fix was attempted in the protein
features pipeline, it's actually much easier to fix here, once the translation
attributes have been calculated.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::CoreStatistics::CorrectNcoils;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Production::Pipeline::Common::Base);

sub run {
  my ($self) = @_;
  
  my $dba = $self->get_DBAdaptor;
  my $dbh = $dba->dbc()->db_handle();
  
  my $correct_ncoils_sql =
  'UPDATE '.
    'translation_attrib INNER JOIN '.
    'attrib_type USING (attrib_type_id) INNER JOIN '.
    'protein_feature USING (translation_id) INNER JOIN '.
    'analysis using (analysis_id) '.
  'SET seq_end = seq_end-1 '.
  'WHERE '.
    'code = "NumResidues" AND '.
    'logic_name = "ncoils" AND '.
    'seq_end = value+1; ';
  
  $dbh->do($correct_ncoils_sql) or $self->throw($dbh->errstr);
}

1;
