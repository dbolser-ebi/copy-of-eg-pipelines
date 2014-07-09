=pod 

=head1 NAME

Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::GeneNamesProjection

=cut

=head1 DESCRIPTION

=head1 MAINTAINER

$Author: ckong $

=cut
package Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::GeneNamesProjection;

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor;
use Bio::EnsEMBL::Utils::SqlHelper;
use base ('Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::ProjectionBase');

sub param_defaults {
    return {
          
	   };
}

my ($flag_store_projections, $flag_backup);
my ($to_species, $from_species, $compara, $release);
my ($method_link_type, $homology_types_allowed, $percent_id_filter);
my ($log_file, $output_dir, $data);
my ($geneName_source, $geneDesc_rules, $taxon_filter);
my ($mlssa, $ha, $ma, $gdba);

sub fetch_input {
    my ($self) = @_;

    $flag_store_projections = $self->param('flag_store_projections');
    $flag_backup            = $self->param('flag_backup');

    $to_species             = $self->param('species');
    $from_species           = $self->param('from_species');
    $compara                = $self->param('compara');
    $release                = $self->param('release');
    $self->throw('to_species is obligatory parameter')   unless (defined $to_species);
    $self->throw('from_species is obligatory parameter') unless (defined $from_species);
    $self->throw('compara is obligatory parameter')      unless (defined $compara);
    $self->throw('release is obligatory parameter')       unless (defined $release);
    
    $method_link_type       = $self->param('method_link_type');
    $homology_types_allowed = $self->param('homology_types_allowed ');
    $percent_id_filter      = $self->param('percent_id_filter');
    $log_file               = $self->param('output_dir');
    $output_dir             = $self->param('output_dir');
    $self->throw('method_link_type is obligatory parameter')       unless (defined $method_link_type);
    $self->throw('homology_types_allowed is obligatory parameter') unless (defined $homology_types_allowed);
    $self->throw('percent_id_filter is obligatory parameter')      unless (defined $percent_id_filter);
    $self->throw('log_file is obligatory parameter')               unless (defined $log_file);
    $self->throw('output_dir is obligatory parameter')             unless (defined $output_dir);

    $geneName_source        = $self->param('geneName_source');
    $geneDesc_rules         = $self->param('geneDesc_rules');
    $taxon_filter           = $self->param('taxon_filter');
    $self->throw('geneName_source is obligatory parameter') unless (defined $geneName_source);
    $self->throw('geneDesc_rules is obligatory parameter') unless (defined $geneDesc_rules);
    $self->throw('taxon_filter is obligatory parameter')    unless (defined $taxon_filter);

return;
}

sub run {
    my ($self) = @_;

    # Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
    # Bio::EnsEMBL::Registry->disconnect_all();

    # Get taxon ancestry of the target species
    my $to_latin_species   = ucfirst(Bio::EnsEMBL::Registry->get_alias($to_species));
    my $meta_container     = Bio::EnsEMBL::Registry->get_adaptor($to_latin_species,'core','MetaContainer');
    my ($to_taxon_id)      = @{ $meta_container->list_value_by_key('species.taxonomy_id')};
    my ($ancestors,$names) = $self->get_taxon_ancestry($to_taxon_id);  

    # Exit projection if 'taxon_filter' is not found in the $ancestor list
    if (!grep (/$taxon_filter/, @$names)){
        die("$taxon_filter is not found in the ancestor list of $to_species\n")
    };

    # Creating adaptors
    my $from_ga   = Bio::EnsEMBL::Registry->get_adaptor($from_species, 'core', 'Gene');
    my $to_ga     = Bio::EnsEMBL::Registry->get_adaptor($to_species  , 'core', 'Gene');
    my $to_ta     = Bio::EnsEMBL::Registry->get_adaptor($to_species  , 'core', 'Transcript');
    my $to_dbea   = Bio::EnsEMBL::Registry->get_adaptor($to_species  , 'core', 'DBEntry');
    die("Problem getting DBadaptor(s) - check database connection details\n") if (!$from_ga || !$to_ga || !$to_ta || !$to_dbea);

    $mlssa = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'MethodLinkSpeciesSet'); 
    $ha    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Homology'); 
    $ma    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Member');   
    $gdba  = Bio::EnsEMBL::Registry->get_adaptor($compara, "compara", 'GenomeDB'); 
    die "Can't connect to Compara database specified by $compara - check command-line and registry file settings" if (!$mlssa || !$ha || !$ma ||!$gdba);

    $self->check_directory($log_file);
    $log_file  = $log_file."/".$from_species."-".$to_species."_GeneNamesProjection_logs.txt";

    # Write projection info metadata
    open $data,">","$log_file" or die $!;
    print $data "\n\tProjection log :\n";
    print $data "\t\trelease             :$release\n";
    print $data "\t\tfrom_db             :".$from_ga->dbc()->dbname()."\n";
    print $data "\t\tfrom_species_common :$from_species\n";
    print $data "\t\tto_db               :".$to_ga->dbc()->dbname()."\n";
    print $data "\t\tto_species_common   :$to_species\n";
   
    # Backup tables that will be updated
    $self->backup($to_ga) if($flag_backup==1);
    #backup($to_ga,$to_species)if($flag_backup==1);
	    
    # Build Compara GenomeDB objects
    my $from_GenomeDB = $gdba->fetch_by_registry_name($from_species);
    my $to_GenomeDB   = $gdba->fetch_by_registry_name($to_species);
    my $mlss          = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$from_GenomeDB, $to_GenomeDB]);
    my $mlss_id       = $mlss->dbID();
	    
    # Get homologies from compara - comes back as a hash of arrays
    print $data "\n\tRetrieving homologies of method link type $method_link_type for mlss_id $mlss_id \n";
    my $homologies    = $self->fetch_homologies($ha, $mlss, $from_species, $data, $gdba, $homology_types_allowed, $percent_id_filter);

    print $data "\n\tProjecting Gene Names & descriptions from $from_species to $to_species\n\n";

    my $total_genes   = scalar(keys %$homologies);

    foreach my $from_stable_id (keys %$homologies) {
       my $from_gene  = $from_ga->fetch_by_stable_id($from_stable_id);
       next if (!$from_gene);
       my @to_genes   = @{$homologies->{$from_stable_id}};
    
       foreach my $to_stable_id (@to_genes) {
          my $to_gene  = $to_ga->fetch_by_stable_id($to_stable_id);
          next if (!$to_gene);
          project_genenames($to_ga, $to_dbea, $from_gene, $to_gene);
       }
    }
    close($data);

return;
}

sub write_output {
    my ($self) = @_;

    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
    Bio::EnsEMBL::Registry->disconnect_all();

}

######################
## internal methods
######################
sub project_genenames {
    my ($to_geneAdaptor, $to_dbea, $from_gene, $to_gene,$ensemblObj_type) = @_;

    # Project when 'source gene' has display_xref and 'target gene' has NO display_xref
    if(defined $from_gene->display_xref() 
       && !defined $to_gene->display_xref()){

       my $from_gene_dbname     = $from_gene->display_xref->dbname();
       my $from_gene_display_id = $from_gene->display_xref->display_id();         

       # Get all DBEntries for 'source gene' base on the dbname of display_xref  
       foreach my $dbEntry (@{$from_gene->get_all_DBEntries($from_gene_dbname)}) { 

          if($dbEntry->display_id=~/$from_gene_display_id/  
               && $flag_store_projections==1 
               && grep (/$from_gene_dbname/, @$geneName_source))
          {
             print $data "\t\tProject from:".$from_gene->stable_id()."\t";
             print $data "to:".$to_gene->stable_id()."\t";
             print $data "GeneName:".$from_gene->display_xref->display_id()."\t";
             print $data "DB:".$from_gene->display_xref->dbname()."\n";

	     # Adding projection source information 
             $dbEntry->info_type("PROJECTION");
             $dbEntry->info_text("projected from $from_species,".$from_gene->stable_id());

             $to_dbea->store($dbEntry,$to_gene->dbID(), 'Gene', 1);
             $to_gene->display_xref($dbEntry);
             $to_geneAdaptor->update($to_gene);
          }
      }
   } 

   # Project gene_description to target_gene
   my $gene_desc = $from_gene->description();

   if(defined $from_gene->description() 
       && !defined $to_gene->description()
       && $flag_store_projections==1 
       && !grep (/$gene_desc/, @$geneDesc_rules)) 
   {
       $gene_desc    = $gene_desc."(projected from $from_species,".$from_gene->stable_id().")";

       print $data "\t\tProject from: ".$from_gene->stable_id()."\t";
       print $data "to: ".$to_gene->stable_id()."\t";
       print $data "Gene Description: $gene_desc\n";
 
       $to_gene->description($gene_desc);
       $to_geneAdaptor->update($to_gene);
   }    
}


1;
