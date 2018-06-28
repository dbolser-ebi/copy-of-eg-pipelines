# Copyright [1999-2014] EMBL-European Bioinformatics Institute
# and Wellcome Trust Sanger Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Bio::EnsEMBL::Utils::CliHelper;
use Carp;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Seq;

#### Functions definition ###

#sub trim { my $s = shift; $s =~ s/^_+|_+$//g; return $s };

sub find_translation { # REVERTING THe ORDER OF THE SEARCHES DID'T SHOWED SIGNIFICANT TIME IMPROVEMENT
  my ( $dba, $acc, $locus, $gene_name ) = @_;
  my $translation;
  my $identifier_type;
  my $dbentry_adaptor = $dba->get_adaptor("DBEntry");
  my @transcripts_ids = $dbentry_adaptor->list_transcript_ids_by_extids($acc);
  my @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($acc);

  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($locus);
    $identifier_type = 'locus';
  }
  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($gene_name);
    $identifier_type = 'name';
  }
  
  my $translation_adaptor = $dba->get_adaptor("Translation");
  my $transcript_adaptor  = $dba->get_adaptor("Transcript");
  my $gene_adaptor        = $dba->get_adaptor("Gene");
  my $transcript;
  my $gene;

  if ( scalar(@transcripts_ids) >= 1 ) {
    my $transcript_id = $transcripts_ids[0];
    $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    $translation = $translation_adaptor->fetch_by_Transcript($transcript);
    $identifier_type = 'accession';
  } elsif ( scalar(@gene_ids) >= 1 ) {
    $gene = $gene_adaptor->fetch_by_dbID( $gene_ids[0] );
    my @transcripts = @{ $transcript_adaptor->fetch_all_by_Gene($gene) };
    $transcript = $transcripts[0];
    $translation = $translation_adaptor->fetch_by_Transcript($transcript);
  }

  return $translation, $identifier_type;
} ## end sub find_translation

my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();

sub get_uniprot_seq {
  my ($acc) = @_;
  my $seq = $search->get('http://www.uniprot.org/uniprot/'.$acc.'.fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

sub get_uniparc_seq {
  my ($acc) = @_;
  my $seq = $search->get('http://www.uniprot.org/uniparc/?query=' . $acc . '&columns=sequence&format=fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

#### Load parameters and data ####

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

#### get the basic options for connecting to a database server
            #
my $optsd = [ @{ $cli_helper->get_dba_opts() }, 
              @{ $cli_helper->get_dba_opts('tax') } ,
              "division:s", 
              "uniprot_file:s", 
              "verbose", 
              "results_file:s", 
              "database_list_file:s", 
              "former_results_file:s" , 
              "blast_db_dir:s"
            ];

my $opts = $cli_helper->process_args( $optsd, \&pod2usage );

$opts->{results_file} ||= $opts->{uniprot_file}.".out";

Bio::EnsEMBL::Registry->load_registry_from_db( 
                         -USER => $opts->{user},
                         -PASS => $opts->{pass},
                         -HOST => $opts->{host},
                         -PORT => $opts->{port} );

print "Registry loaded\n";

my $lookup =
  Bio::EnsEMBL::LookUp::LocalLookUp->new( -SKIP_CONTIGS => 1,
                                          -NO_CACHE     => 1 ); # Used to check all leafs sub_specie/strains under a taxonomy_id (specie)

my ($tax_dba_details) = @{ $cli_helper->get_dba_args_for_opts( $opts, 1, 'tax' ) };
$tax_dba_details->{-GROUP}   = 'taxonomy';

my $tax_adaptor =
  Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor->new(
           Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$tax_dba_details} ) );

$lookup->taxonomy_adaptor($tax_adaptor);

print "Taxonomy lookup sevice loaded\n";

#### Load the PHI-base file

my @phibase_data;#contains each phi_base entry correctly loaded form the input phi file 
my $total_phi = 0;
my %mapped_annotations;

open my $uniprot, "<", $opts->{uniprot_file} or croak "Could not open ".$opts->{uniprot_file};

while( my $line = <$uniprot>) {
      chomp $line;
      $line =~ s/\r//g; 
      my @fields = split( ",", $line, -1 );

      if ( scalar(@fields) == 15) {
        push @phibase_data, \@fields;
        $total_phi++;
      } 
}

# 15 - Direct map ? Yes or No
# 16 - Sequence match ? Yes or No
# 17 - Database name where is matched gene
# 18 - Gene stable id
# 19 - Translation stable id
# 20 - Percentage identity

print "Total number of annotations: " . $total_phi . "\n";
print "PHI-base csv file loaded\n";

if ( $opts->{former_results_file} ) { 
  open my $former_results, "<", $opts->{former_results_file} or croak "Could not open ".$opts->{former_results_file};

  while(my $line =<$former_results>) {
    chomp $line;
    my @fields = split( ",", $line, -1 );
    $mapped_annotations{$fields[1]} = [@fields[15,16,17,18,19,20]];
  }
  print "Former results loaded\n";
} else {
  print "Former results not available\n";
}

#### Process each line of the PHI-base file

my $count_phi_zero_dbs = 0;
my $count_phi_many_dbs = 0;
my $uniprots = {};

open my $output, ">", $opts->{results_file} or croak "Could not open ".$opts->{results_file};

foreach my $annotation ( @phibase_data ) {

  my $annotn_tax_id = ${$annotation}[6];
  my $annotn_phi_entry = ${$annotation}[1];#PHI-base entry
  my $annotn_pathogene_specie = ${$annotation}[8];
  my $annotn_acc = ${$annotation}[2];#Uniprot protein id
  my $annotn_gene_name = ${$annotation}[3];#gene name
  my $annotn_locus = ${$annotation}[4];#gene_id

  my $result_line = join(',', @{$annotation});

  # print "\nANNOTATION:\n" . $annotn_phi_entry . "\nuniprot accession:" . $annotn_acc . "\nphi_base _gene_name:" . $annotn_gene_name .  "\ngene_locus_id:"
  #  . $annotn_locus . "\npathogen_species_taxon_id:" . $annotn_tax_id . "\npathogen_taxon_name_ncbi:" . $annotn_pathogene_specie . "\n******\n";

  # Get all available genomes for a given taxon id

  print "\nProcessing $annotn_phi_entry :\n";
  
  if ( $annotn_tax_id eq '' || ${$annotation}[0] < 0) {
    print "\tTax id is empty";
    next;
  }

  my $branch_dbas = $lookup->get_all_by_taxon_branch($annotn_tax_id);
  my $nb_taxon_branch_dbs = scalar(@{$branch_dbas});

  if ($nb_taxon_branch_dbs == 0) {
    print "\tNo dbs for this specie.\n";
    $count_phi_zero_dbs++;
    next;
  } elsif ( $nb_taxon_branch_dbs > 30 ) {
      print "\tToo many dbs for this specie.\n";
      $count_phi_many_dbs++;
      next;
  }

  if ( !exists($mapped_annotations{$annotn_phi_entry}) &&
       ($nb_taxon_branch_dbs <= 10 && $nb_taxon_branch_dbs > 0) )#$nb Of SubSpecies between 1 and 10
  {
    for my $brch_db (@{$branch_dbas}) {
      my $mc = $brch_db->get_MetaContainer();
      my $branch_species = $mc->single_value_by_key('species.production_name');
      my $division =  $mc->single_value_by_key('species.division');
      print $brch_db->dbc()->dbname() . "\n";

      if ( $division eq $opts->{division} ) {
        print "\tFound genome for species " . $branch_species . " in " . $brch_db->dbc()->dbname() . "/" . $brch_db->species_id ."\n";
        ( my $translation , my $id_type ) = find_translation( $brch_db, $annotn_acc, $annotn_locus, $annotn_gene_name );

        if ($translation) {
          print "\t\tFound gene on ensembl: \n";
          my $gene_adaptor        = $brch_db->get_adaptor("Gene");
          my $gene = $gene_adaptor->fetch_by_translation_stable_id($translation->stable_id);
          
          print "\t\t" . $gene->stable_id . "\t" . $translation->stable_id . "\n";
          print $output $result_line . ',1,' . $id_type . ',0,' . $branch_species . "," . $gene->stable_id . ',' . $translation->stable_id . ',' . 100 . "\n";

          last;

        } else {
          print "\t\tDidn't find direct match on Ensembl \n";
          eval {
            my $up = get_uniprot_seq($annotn_acc);
            $uniprots->{$annotn_acc} = $up;
            $id_type = 'uniprot';
          };
          if($@ || $uniprots->{$annotn_acc}->{seq} eq '') {
            print "Could not find entry : " . $annotn_acc;
            my $up = get_uniparc_seq($annotn_acc);
            $uniprots->{$annotn_acc} = $up;
            $id_type = 'uniparc';
          }

          print "\t\tGot uniprot/uniparc accession: \n";
          print Dumper($uniprots->{$annotn_acc});
          my @translations;

          if ($uniprots->{$annotn_acc}->{seq} ne '') {
            print "Proceeding with BLAST search \n";
            my $branch_species = $mc->single_value_by_key('species.production_name');
            my $blast_db_dir = $opts->{blast_db_dir};

            # Create BLAST DB for this specie
            print "creating peptide file if it doesn exist already in $blast_db_dir\n";
            my $blast_database = $blast_db_dir . "/" .$branch_species . ".pep.all";
            print $blast_database . "\n";

            if ( ! -f $blast_database ) { # Returns undef if the file doesn't exist, or if the file exists but it's not plain file (e.g. if it's a directory).
              open (my $db_file, ">", $blast_database) or croak "Could not open ".$blast_database;;
              print "getting all genes for $branch_species \n";
              my $translation_adaptor=$brch_db->get_TranslationAdaptor();
              my @translations = @{$translation_adaptor->fetch_all()}; 

              foreach my $translation (@translations) {
                my $stable_id = $translation->stable_id();
                my $sequence = $translation->seq();
                print $db_file ">" . $stable_id . "\n";
                print $db_file $sequence . "\n";
              }
              close $db_file;
            }
            my $fac = Bio::Tools::Run::StandAloneBlastPlus->new( # THIS -db_dir should not be hardcoded, make it a parameter instead
               -db_name => $branch_species,
               -db_dir => $blast_db_dir ,
               -db_data => $blast_database,
               -create => 1
            );
            $fac->make_db();

            print "BLAST database created for branch_species in $blast_database\n";

            my $query   = Bio::Seq->new( -display_id => $annotn_acc, -seq =>  $uniprots->{$annotn_acc}->{seq} );
            my $results = $fac->blastp( -query => $query, -method_args => [ -num_alignments => 1 ]);               
            my $query_length = length($uniprots->{$annotn_acc}->{seq});
            my $hit = $results->next_hit();
            my $hit_length = $hit->length();
            my $translation_stable_id = $hit->name();
            my $evalue = $hit->hsp->evalue();;              
            my $identity = $hit->hsp->percent_identity(); 
        
            print Dumper($evalue);               
            print Dumper($identity);

            if ( !$results ) {
              next;
            }
          
            my $query_length_covered = $hit_length / $query_length ; 

            print $query_length_covered . "\n";
            
            if ( $query_length_covered > 0.8) {
              print $output $result_line . ',1,' . $id_type . ',1,' . $branch_species . "," . 
                 "GENE_NAME" . "," . $translation_stable_id . "," . $evalue . ",". $identity . "\n" ;
            }
          }
        }
      }
    }
  }
}

close $output;

print "\n";
print "Total number of annotations with no species equivalent in Ensembl: " . $count_phi_zero_dbs . "\n";
print "Total number of annotations with too many species equivalent in Ensembl: " . $count_phi_many_dbs . "\n";

print "Execution time: $duration s\n";






