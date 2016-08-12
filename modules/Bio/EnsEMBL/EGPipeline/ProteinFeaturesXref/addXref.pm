package Bio::EnsEMBL::EGPipeline::ProteinFeaturesXref::addXref;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  return {
   'external_db' => 'Uniprot/SPTREMBL',
  };
}

sub run {
  my ($self) = @_;
  my $logic_name = $self->param_required('logic_name');
  
  my $dba  = $self->core_dba();
  my $pfa  = $dba->get_adaptor('ProteinFeature');
  my $aa   = $dba->get_adaptor('Analysis');
  my $dbea = $dba->get_adaptor('DBEntry');
  
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  my @ProteinFeaturesArray = @{ $pfa->fetch_all_by_logic_name($logic_name) };
  
  foreach my $feature (@ProteinFeaturesArray){
    $self->add_xref($dbea, $feature, $analysis);
  }
}

sub add_xref {
  my ($self, $dbea, $feature, $analysis) = @_;
  my $external_db = $self->param_required('external_db');
  my $hit_name = $feature->hseqname();
  
  my $xref = Bio::EnsEMBL::DBEntry->new
    (
      -dbname      => $external_db,
      -primary_id  => $hit_name,
      -display_id  => $hit_name,
    );
  $xref->analysis($analysis);
  $dbea->store($xref, $feature->translation_id(), 'Translation');
}

1;
