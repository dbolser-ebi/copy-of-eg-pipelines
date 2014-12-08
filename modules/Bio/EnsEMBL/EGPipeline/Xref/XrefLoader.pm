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
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis;
use Digest::MD5;
=head1 CONSTRUCTOR
=head2 new
  Arg [-NAME]  : 
       string - human readable version of the name of the genome
  Arg [-SPECIES]    : 
       string - computable version of the name of the genome (lower case, no spaces)
  Arg [-DBNAME] : 
       string - name of the core database in which the genome can be found
  Arg [-SPECIES_ID]  : 
       int - identifier of the species within the core database for this genome
  Arg [-TAXONOMY_ID] :
        string - NCBI taxonomy identifier
  Arg [-ASSEMBLY_NAME] :
        string - name of the assembly
  Arg [-ASSEMBLY_ID] :
        string - INSDC assembly accession
  Arg [-ASSEMBLY_LEVEL] :
        string - highest assembly level (chromosome, supercontig etc.)
  Arg [-GENEBUILD]:
        string - identifier for genebuild
  Arg [-DIVISION]:
        string - name of Ensembl Genomes division (e.g. EnsemblBacteria, EnsemblPlants)
  Arg [-STRAIN]:
        string - name of strain to which genome belongs
  Arg [-SEROTYPE]:
        string - name of serotype to which genome belongs

  Example    : $info = Bio::EnsEMBL::Utils::MetaData::GenomeInfo->new(...);
  Description: Creates a new info object
  Returntype : Bio::EnsEMBL::Utils::MetaData::GenomeInfo
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub new {
  my ($proto, @args) = @_;
  my $class = ref($proto) || $proto;
  my $self = bless({}, $class);
  $self->{logger} = get_logger();
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
  if(!defined $self->{logger}) {
  	  $self->{logger} = get_logger();
  }
  return $self->{logger};
}

=head2 get_analysis
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Arg        : logic name
  Description: Get specified analysis object
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : internal
  Status     : Stable
=cut
sub get_analysis {
  my ($self, $dba, $logic_name) = @_;
  
  my $aa = $dba->get_AnalysisAdaptor();
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  if (!defined $analysis) {
    $analysis = Bio::EnsEMBL::Analysis->new(
      -logic_name => $logic_name,
    );
    $aa->store($analysis);
    $analysis = $aa->fetch_by_logic_name($logic_name);
    if (!defined $analysis) {
      $self->logger()->warn("Analysis $logic_name could not be added to the core database.");
    }
  }
  return $analysis;
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
