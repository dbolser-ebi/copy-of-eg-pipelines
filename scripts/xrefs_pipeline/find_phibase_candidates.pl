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
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::EGPipeline::Xref::Needle;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = [ @{ $cli_helper->get_dba_opts() }, "uniprot_file:s", "verbose", "results_file:s" ];

my $opts = $cli_helper->process_args( $optsd, \&pod2usage );
if ( $opts->{verbose} ) {
  Log::Log4perl->easy_init($DEBUG);
}
else {
  Log::Log4perl->easy_init($INFO);
}
$opts->{results_file} ||= $opts->{uniprot_file}.".out";

my $logger = get_logger();

my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
my $needle = Bio::EnsEMBL::EGPipeline::Xref::Needle->new();

$logger->info("Reading from ".$opts->{uniprot_file});
open my $uniprot, "<", $opts->{uniprot_file} or croak "Could not open ".$opts->{uniprot_file};

my $uniprots = {};
while(<$uniprot>) {
      chomp;
      my ($phi,$acc) = split;
      # fetch uniprot seq
      $logger->info("Fetching UniProt entry ".$acc);
      eval {
	  my $up = get_uniprot_seq($acc);
	  $uniprots->{$acc} = $up;
	  $up->{phi} = $phi; 
      };
      if($@) {
	  $logger->warn("Could not find entry ".$acc);
      }
}
close $uniprot;


open my $out, ">", $opts->{results_file} or croak "Could not open results file ".$opts->{results_file};
$logger->info("Connecting to core database(s)");
for my $core_dba_details (@{$cli_helper->get_dba_args_for_opts($opts)}) {
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$core_dba_details});
  my $ta = $dba->get_TranslationAdaptor();
  my $mc = $dba->get_MetaContainer();
  (my $div = lc($mc->get_division())) =~ s/ensembl//;
  my $ass = $mc->single_value_by_key('assembly.default');
  my $species = $mc->single_value_by_key('species.url');
  my $collection = "ensemblgenomes/$div/$species.$ass.pep.all";
  $logger->info("Searching $collection");

  my %input = map {$_ => $uniprots->{$_}->{seq} } keys %$uniprots;
  my $results = $search->search(\%input,$collection,'protein','blastp',{exp      => '1e-5', scores=>5, alignments=>5});

# get best hits in each case to allow for multi family
  my $best_hits = {};
  my $ts = {};
  my $needle_input = {};
  while(my($query,$qres) = each %$results) {
      my $qStr = $uniprots->{$query}->{seq};
      my $len = length($qStr);
      my $max_len = 1.1*$len; 
      my $min_len = 0.9*$len; 
      for my $qr (sort { $a->{expectation} <=> $b->{expectation}} grep {$_->{length}<=$max_len && $_->{length}>=$min_len} values %$qres) {
          push @{$best_hits->{$query}}, $qr->{id};
          my $t = $ta->fetch_by_stable_id($qr->{id});
          $ts->{$qr->{id}} = $t;
          $needle_input->{$query."-".$qr->{id}} = {asequence=>">".$query."\n".$qStr,bsequence=>">".$qr->{id}."\n".$t->seq()};
      }
  }
  my $needle_results = $needle->align($needle_input);
  while(my($query,$hits) = each %$best_hits) {
      my $outstr = $uniprots->{$query}->{phi}. " - annotated with UniProt:$query (".$uniprots->{$query}->{des}.")\n";
      my $i = 0;
      for my $hit (@{$hits}) {
          my $t = $ts->{$hit};
          my $g = $t->transcript()->get_Gene();
          my $gene_str = "Gene ".$g->stable_id();
          my $nom =  $g->external_name();
          $gene_str.=" ($nom)" if(defined $nom);
          $gene_str.=" ".$g->description() if(defined $nom);
          my $n = $needle_results->{$query.'-'.$hit};
	  if($n->{id} >= 90) {
	      $outstr.= "Hit ".++$i.": ".$n->{id}."% identity (".$n->{pos}."/".$n->{len}.") to $hit ($gene_str)";
	      $outstr.= $n->{aln}."\n";
	  }
      }
      print $out $outstr if($i>0);
  }
}
close $out;

sub get_uniprot_seq {
    my ($acc) = @_;
    my $seq = $search->get('http://www.uniprot.org/uniprot/'.$acc.'.fasta');    
    $seq =~ s/^>\S+\s+([^\n]+)\n//;
    my $des = $1;
    $seq =~ tr/\n//d;
    return {seq=>$seq,des=>$des};
}
