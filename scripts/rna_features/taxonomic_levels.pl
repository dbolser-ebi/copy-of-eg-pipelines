#!/usr/bin/env/perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;

# This script uses a file that maps Rfam accessions to taxon IDs, in order to
# determine the numbers of sequences assigned to an Rfam accession at a
# set of taxonomic levels. The path to the 'rfam2taxonomy' file is mandatory.
# 
# If the 'rfam2taxonomy' file does not exist it will be generated, and a
# further two files become mandatory, the 'full_region' and 'rfamseq' files.
# The full_region_file should be downloaded and unzipped from: 
#  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full_region.gz
# The rfamseq_file should be downloaded and unzipped from:
#  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/rfamseq.txt.gz
# Once the 'rfam2taxonomy' file has been created it should be saved (as ensgen)
# so that others don't need to recreate it:
#  /nfs/panda/ensemblgenomes/external/Rfam/<rfam_release>/rfam2taxonomy.txt
# 
# With the 'rfam2taxonomy' file it is then possible to use the mappings to
# tot up the numbers of sequences for each Rfam accession at a given set of
# taxonomic levels. By default, these levels correspond to the EG divisions,
# but any set of levels can be specified, using the 'level' parameter as many
# times as necessary.
# 
# The counts are sent to STDOUT; if using the default EG levels, please save
# the output file (as ensgen) so that others don't need to recreate it:
#  /nfs/panda/ensemblgenomes/external/Rfam/<rfam_release>/taxonomic_levels.txt
#

my ($host, $port, $user, $pass, $dbname, @levels, $root, $count_seqs, 
    $rfam2taxonomy_file, $full_region_file, $rfamseq_file, $tmp_dir, );

GetOptions(
  "host=s", \$host,
  "P|port=i", \$port,
  "user=s", \$user,
  "p|pass=s", \$pass,
  "dbname=s", \$dbname,
  "levels=s", \@levels,
  "root=s", \$root,
  "count_seqs", \$count_seqs,
  "rfam2taxonomy_file=s", \$rfam2taxonomy_file,
  "full_region_file=s", \$full_region_file,
  "rfamseq_file=s", \$rfamseq_file,
  "tmp_dir=s", \$tmp_dir,
);

die "-rfam2taxonomy_file is required" unless $rfam2taxonomy_file;

$host   = 'mysql-eg-pan-prod.ebi.ac.uk' unless $host;
$port   = '4276'                        unless $port;
$user   = 'ensro'                       unless $user;
$pass   = ''                            unless $pass;
$dbname = 'ncbi_taxonomy'               unless $dbname;

my $tax_dba = Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
  -host   => $host,
  -port   => $port,
  -user   => $user,
  -pass   => $pass,
  -dbname => $dbname,
);

my $nta = $tax_dba->get_TaxonomyNodeAdaptor();

my @eg_levels;
unless (scalar @levels) {
  @levels = qw(
    Archaea
    Bacteria
    Eukaryota
      Fungi
      Metazoa
        Chordata
      Viridiplantae
  );
  
  @eg_levels = qw(
    EnsemblBacteria
    EnsemblFungi
    EnsemblMetazoa
    EnsemblPlants
    EnsemblProtists
    Ensembl
  );
}

$root = 'cellular organisms' unless $root;

$count_seqs = 0 unless $count_seqs;
my $count_type = $count_seqs ? 'Sequences' : 'Taxa';

if (-e $rfam2taxonomy_file) {
  print STDERR "Mapping file '$rfam2taxonomy_file' exists, so won't be regenerated.\n";
} else {
  $tmp_dir = '/tmp'   unless $tmp_dir;
  make_path($tmp_dir) unless -e $tmp_dir;
  
  my $rfam2rfamseq_file = catdir($tmp_dir, 'rfam2rfamseq.txt');
  my $rfamseq2taxonomy_file = catdir($tmp_dir, 'rfamseq2taxonomy.txt');
  
  die '-full_region_file is required and must exist' unless $full_region_file && -e $full_region_file;
  die '-rfamseq_file is required and must exist' unless $rfamseq_file && -e $rfamseq_file;

  # The input files are huge, so trim down on the command line.
  my $cmd_1 = "cut -f1-2 $full_region_file | sort -k 2,2 -k 1,1 -u > $rfam2rfamseq_file";
  my $cmd_2 = "cut -f1,4 $rfamseq_file | sort -u > $rfamseq2taxonomy_file";
  my $cmd_3 = "join -1 2 -2 1 -t \$'\\t' $rfam2rfamseq_file $rfamseq2taxonomy_file > $rfam2taxonomy_file";

  system($cmd_1) == 0 || die "Failed to run: $cmd_1";
  system($cmd_2) == 0 || die "Failed to run: $cmd_2";
  system($cmd_3) == 0 || die "Failed to run: $cmd_3";
  
  unlink $rfam2rfamseq_file;
  unlink $rfamseq2taxonomy_file;
}

my $rfam2taxonomy_path = path($rfam2taxonomy_file);
my $rfam2taxonomy = $rfam2taxonomy_path->slurp;

my %rfam2taxonomy;
my %taxa;
foreach my $line (split(/\n/, $rfam2taxonomy)) {
  my (undef, $rfam_acc, $taxon_id) = split(/\s/, $line);
  $rfam2taxonomy{$rfam_acc}{$taxon_id}++;
  $taxa{$taxon_id} = undef;
}
my @rfam_acc = sort keys %rfam2taxonomy;

foreach my $taxon_id (keys %taxa) {
  my $node = $nta->fetch_by_taxon_id($taxon_id);
  
  if (! defined($node)) {
    # Because our NCBI taxonomy database is more recent than that used to
    # make the Rfam alignments, there are occasionally deleted or merged
    # taxon IDs. So just skip over them.
    foreach my $rfam_acc (@rfam_acc) {
      delete $rfam2taxonomy{$rfam_acc}{$taxon_id};
    }
    
  } elsif (! has_ancestor($node, $root)) {
    # Exclude species above the given root. (With the default settings,
    # this excludes viral and metagenomic sequences.)
    foreach my $rfam_acc (@rfam_acc) {
      delete $rfam2taxonomy{$rfam_acc}{$taxon_id};
    }
  }
}

my @columns = ('Rfam_acc', $count_type, 'LCA', @levels, @eg_levels);
print join("\t", @columns)."\n";

foreach my $rfam_acc (@rfam_acc) {
  print STDERR "Processing $rfam_acc\n";
  
  my $lca;
  my $total = 0;
  my %levels;
  foreach my $level (@levels) {
    $levels{$level} = 0;
  }
  
  foreach my $taxon_id (sort keys %{$rfam2taxonomy{$rfam_acc}}) {
    my $subtotal = $rfam2taxonomy{$rfam_acc}{$taxon_id};
    
    my $node = $nta->fetch_by_taxon_id($taxon_id);
    
    # Keep track of the last common ancestor of the species we've seen.
    if (defined $lca) {
      if (! has_ancestor($node, $lca->name)) {
        if ($node->name ne $lca->name) {
          $lca = $nta->fetch_common_ancestor($lca, $node);
        }
      }
    } else {
      $lca = $node;
    }
    
    foreach my $level (@levels) {
      if (has_ancestor($node, $level)) {
        if ($count_seqs) {
          $levels{$level} += $subtotal;
        } else {
          $levels{$level}++;
        }
      }
    }
    
    if ($count_seqs) {
      $total += $subtotal;
    } else {
      $total++;
    }
  }
  
  my %eg_levels;
  if (scalar @eg_levels) {
    $eg_levels{'EnsemblBacteria'} = $levels{'Archaea'}
                                  + $levels{'Bacteria'};
    $eg_levels{'EnsemblFungi'}    = $levels{'Fungi'};
    $eg_levels{'EnsemblMetazoa'}  = $levels{'Metazoa'}
                                  - $levels{'Chordata'};
    $eg_levels{'EnsemblPlants'}   = $levels{'Viridiplantae'};
    $eg_levels{'EnsemblProtists'} = $levels{'Eukaryota'}
                                  - $levels{'Fungi'}
                                  - $levels{'Metazoa'}
                                  - $levels{'Viridiplantae'};
    $eg_levels{'Ensembl'}         = $levels{'Chordata'};
  }
  
  my $lca_name;
  if (defined $lca) {
    $lca_name = $lca->name;
  } else {
    $lca_name = "beyond root ($root)";
  }
  
  my @level_counts    = map { $levels{$_} } @levels;
  my @eg_level_counts = map { $eg_levels{$_} } @eg_levels;
  
  my @row = (
    $rfam_acc,
    $total,
    $lca_name,
    @level_counts,
    @eg_level_counts,
  );
  
  print join("\t", @row)."\n";
  
  delete $rfam2taxonomy{$rfam_acc};
}

sub has_ancestor {
  my ($node, $name) = @_;
  
  $node = $node->parent;
  while ($node) {
    return 1 if ($node->name eq $name);
    $node = $node->parent;
  }
  return 0;
}
