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

Bio::EnsEMBL::EGPipeline::Xref::XrefLoader

=head1 DESCRIPTION

Base class for all other loaders in Bio::EnsEMBL::EGPipeline::Xref, providing some common methods

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Xref::XrefLoader;
use Log::Log4perl qw/:easy/
use base Bio::EnsEMBL::EGPipeline::BaseLoader;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis;
use Digest::MD5;

=head1 CONSTRUCTOR
=head2 new
  Example    : $info = Bio::EnsEMBL::EGPipeline::Xref::XrefLoader->new();
  Description: Creates a new loader object
  Returntype : Bio::EnsEMBL::EGPipeline::Xref::XrefLoader
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut
sub new {
  my ($proto, @args) = @_;
  my $self = $proto->SUPER::new(@args);
  return $self;
}

=head2 get_translation_uniprot
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Get mapping of UniProt to translation
  Returntype : Hash of translation IDs to UniProt accessions
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub get_translation_uniprot {
  my ($self, $dba) = @_;
  # get hash of xrefs by translation ID
  my $translation_accs = {};
  my $dbea             = $dba->get_DBEntryAdaptor();
  $dba->dbc()->sql_helper()->execute_no_return(
	-SQL => q/
		select tl.translation_id,unix.xref_id
		from
	translation tl
	join transcript tr using (transcript_id)
	join seq_region sr using (seq_region_id)
	join coord_system cs using (coord_system_id)
	join object_xref uniox on (uniox.ensembl_object_type='Translation' and uniox.ensembl_id=tl.translation_id)
	join xref unix using (xref_id) 
	join external_db unie using (external_db_id) 
        where
	cs.species_id=? and 
	unie.db_name in ('Uniprot\/SWISSPROT','Uniprot\/SPTREMBL')
	/,
	-CALLBACK => sub {
	  my ($tid, $xid) = @{$_[0]};
	  $translation_accs->{$tid} = $dbea->fetch_by_dbID($xid);
	  return;
	},
	-PARAMS => [$dba->species_id()]);
  return $translation_accs;
} ## end sub get_translation_uniprot

1;
