use Test::More;

use Bio::EnsEMBL::EGPipeline::RNAFeatures::GeneFilter;

#test 1
use_ok(Bio::EnsEMBL::Gene);

#test 2
use_ok(Bio::EnsEMBL::Exon);

#test 3
use_ok(Bio::EnsEMBL::Transcript);

#test 4
use_ok(Bio::EnsEMBL::Registry);

#test 5
use_ok(Bio::EnsEMBL::DBSQL::DBAdaptor);

#test 6
use_ok(Bio::EnsEMBL::DBEntry);

#test 7
can_ok(Bio::EnsEMBL::EGPipeline::RNAFeatures::GeneFilter , param_defaults);

#test 8
can_ok(Bio::EnsEMBL::EGPipeline::RNAFeatures::GeneFilter , run);

#test 9
can_ok(Bio::EnsEMBL::EGPipeline::RNAFeatures::GeneFilter , make_gene);

#test10
my $exon = Bio::EnsEMBL::Exon->new;
isa_ok( $exon, 'Bio::EnsEMBL::Exon' );

#test11
my $transcript = Bio::EnsEMBL::Transcript->new;
isa_ok( $transcript, 'Bio::EnsEMBL::Transcript' );

#test12
my $gene = Bio::EnsEMBL::Gene->new;
isa_ok( $gene, 'Bio::EnsEMBL::Gene' );

#test12
my $gene_xref = Bio::EnsEMBL::DBEntry->new;
isa_ok( $gene_xref, 'Bio::EnsEMBL::DBEntry' );


done_testing();
