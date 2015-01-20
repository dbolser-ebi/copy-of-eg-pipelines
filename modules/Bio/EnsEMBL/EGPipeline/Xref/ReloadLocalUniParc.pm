
=head1 LICENSE

Copyright [1999-2015] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::EGPipeline::Xref::ReloadLocalUniParc

=head1 DESCRIPTION

Add UniParc checksums to local database.

=head1 Author

James Allen

=cut

use strict;
use warnings;

package Bio::EnsEMBL::EGPipeline::Xref::ReloadLocalUniParc;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub param_defaults {
  return {
    'ftp_file' => 'upidump.lis',
    'ftp_dir'  => '/ebi/ftp/pub/contrib/uniparc',
    'tmp_dir'  => '/tmp',
  };
}

sub run {
  my ($self) = @_;
  
  my $ftp_file    = $self->param_required('ftp_file');
  my $ftp_dir     = $self->param_required('ftp_dir');
  my $tmp_dir     = $self->param_required('tmp_dir');
  my $uniparc_db  = $self->param_required('uniparc_db');
  my $uniparc_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$uniparc_db);
  my $uniparc_dbh = $uniparc_dba->dbc->db_handle();
  
  `cp $ftp_dir/$ftp_file $tmp_dir/$ftp_file`;
  
  my @load_sql = (
    "DROP INDEX md5_idx ON protein;",
    "TRUNCATE TABLE protein;",
    "LOAD DATA LOCAL INFILE '$tmp_dir/$ftp_file' INTO TABLE protein FIELDS TERMINATED BY ' ';",
    "CREATE INDEX md5_idx ON protein(md5);",
  );
  
  foreach my $load_sql (@load_sql) {
    $uniparc_dbh->do($load_sql) or $self->throw("Failed to execute: $load_sql");
  }
  
}

1;
