#!/usr/local/ensembl/bin/perl  -w

use strict;
use SGD;
use SGDConf;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host => $SGD_DBHOST,
   -user => $SGD_DBUSER,
   -dbname => $SGD_DBNAME,
   -pass  => $SGD_DBPASS,
   -port  => $SGD_DBPORT,
  );


# adding assembly type to meta table


my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($SGD_LOGIC_NAME);
my $nc_analysis = $analysis_adaptor->fetch_by_logic_name('ncRNA');
my $operon_analysis = $analysis_adaptor->fetch_by_logic_name('operon');

my $sa = $db->get_SliceAdaptor;
my $ga = $db->get_GeneAdaptor;
my @chromosomes = @{$sa->fetch_all('chromosome')};
my %chromosomes_hash;
foreach my $chr (@chromosomes){
  $chromosomes_hash{$chr->seq_region_name} = $chr;
}
my ($genes, $operons,$xrefs) = &parse_gff(
				   $SGD_GFF_FILE,
				   \%chromosomes_hash,
				   $analysis,
				   $nc_analysis,
				   $operon_analysis,
				  );

write_genes($genes, $db,$xrefs);
write_simple_features($operons,$db);


my @genes = @{$ga->fetch_all};
open(TRANSLATE, "+>>".$SGD_NON_TRANSLATE) or die "couldn't open ".$SGD_NON_TRANSLATE." $!";
TRANSLATION: foreach my $gene(@genes){
  my $translation = &translation_check($gene);
  if($translation){
    next TRANSLATION;
  }else{
    print TRANSLATE $gene->stable_id." doesn't translate\n";
    next TRANSLATION;
  }
}
close(TRANSLATE);

