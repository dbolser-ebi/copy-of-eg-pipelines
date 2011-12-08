#!/sw/arch/bin/perl -w

=head usage

cd /nfs/panda/ensemblgenomes/production/ncgenes_pipelines 

default_db_details="--host=mysql-eg-devel-3.ebi.ac.uk --port=4208 --user=ensrw --password=[replace_me]"
db=ashbya_gossypii_core_11_64_1

perl scripts/generate_ncrna_stable_ids_ashbya_gossypii.pl \
    ${default_db_details} \
    --dbname ${db} \
    --start EFAGO00000000000
=cut

use strict;
use DBI;
use Getopt::Long;
my $port = 3306; 
my ($host, $dbname, $user, $pass, @types, $start);

GetOptions('dbuser|user=s'       => \$user,
	   'dbpass|password=s'   => \$pass,
	   'dbhost|host=s'       => \$host,
	   'dbport|port=i'       => \$port,
	   'dbname=s'            => \$dbname,
	   'start=s'             => \$start,   # Pattern - USE ENS000001 or ENS for human, ENSMUS00001 or ENSMUS for mouse etc 
                                               # don't add G/T/E/P for specific types !!!
	   'types=s'             => \@types,
	   'help'                => sub { usage(); exit(0); },
            );

if (@types == 0) {
    @types = ('gene','transcript','translation','exon');
}

print STDERR "Types: @types\n";

if (!$user || !$host || !$dbname ) {
  usage();
  exit(1);
}

my $dbi = DBI->connect( "DBI:mysql:host=$host:port=$port;database=$dbname", $user, $pass,
			{'RaiseError' => 1}) || die "Can't connect to database\n";

foreach my $type (@types) {
  my $table = $type . "_stable_id";
  my $sth;

  # create and insert new stable_id;
  # get starting stable ID, either specified or current max
  my $new_stable_id=$start;

  # Get the list of objects which don't have a stable_id yet

  my $sql = "SELECT $type.${type}_id FROM $type WHERE stable_id like \"%RNA%\";";
  
  print STDERR "fetching the list of objects without a stable_id for type, $type\n";
  print STDERR "sql: $sql\n";

  $sth = $dbi->prepare($sql);
  $sth->execute();

  while (my ($object_id) = $sth->fetchrow_array()) { 

      print STDERR "processing object, $object_id...\n";

      # scalar localtime gives: Thu May 21 15:06:18 2009
      # should be: 2006-06-26 10:44:34.0

      use POSIX qw(strftime);
      my  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
      $mon++; 
      if ($mon < 10) {
	  $mon = "0" . $mon;
      }
      my $created_date = strftime "%Y-$mon-%e %H:%M:%S", localtime;
      my $modified_date = strftime "%Y-$mon-%e %H:%M:%S", localtime;

      ($new_stable_id,my $nis) = @{increment_stable_id($new_stable_id,$type)};
    
      my $sql = qq~UPDATE ${type} SET stable_id = "$nis", version = 1, created_date = NOW(), modified_date = NOW() WHERE ${type}_id = $object_id;~;
      
      $dbi->do($sql);
  }
  
  $sth->finish();

}


# --------------------------------------------------------------------------------

sub increment_stable_id {

  my $stable_id = shift;
  my $type = shift ;   
  
  my ($prefix, $suffix);  

  # check stable_id format ...  
  if ( $stable_id =~m/([a-zA-Z]+)([0-9]+)/ ){
    #($prefix,$suffix) = $stable_id =~ /([a-zA-Z]+)([0-9]+)/;
    $prefix = $1;
    $suffix = $2;
  } elsif( $stable_id =~m/([a-zA-Z]+)/){
    $prefix = $stable_id ; 
    $suffix = 0;
  }else { 
      print STDERR "stable_id: $stable_id\n";
      die "unrecongnized stable_id format - should match ([a-zA-Z]+)([0-9]+) or ([a-zA-Z]+) !!\n"; 
  } 
  my $new_sid; 

  if ($type=~m/gene/){ 
     $new_sid=$prefix."G"; 
  } elsif ($type=~m/transcript/){ 
     $new_sid=$prefix."T";
  } elsif ($type=~m/translation/){ 
     $new_sid=$prefix."P";
  } elsif ($type=~m/exon/){ 
     $new_sid=$prefix."E";
  }

  print STDERR "suffix: $suffix\n";

  my $new_stable_id = sprintf "%s%011d", $new_sid , $suffix+1 ;   

  my $old = sprintf "%s%011d", $prefix, $suffix+1 ;    
  #return [$old, $new_stable_id] ; 
  my @tmp;  
  
  push @tmp, $old ; 
  push @tmp, $new_stable_id ;  
 
  return \@tmp ; 

}

# -------------------------------------------------------------------------------- 
sub get_max_stable_id_from_gene_archive { 
  my ($dbi, $type) = @_;   

  # try to get from relevant archive
  my $sth = $dbi->prepare("SELECT MAX($type) FROM gene_archive");
  $sth->execute();  

  my $rs ;  
  if (my @row = $sth->fetchrow_array ) { 
    $rs = $row[0];
  }  
  if (length($rs)> 0 ) { 
   return $rs ; 
  } else { 
   print STDERR "no entry for $type found in gene_archive table - returning undef\n" ;  
   return undef ;
  }  

} 



sub usage {

  print << "EOF";

  USAGE :  

  generate_stable_ids.pl -dbuser|user {user} 
                         -dbpass|pass {password} 
                         -dbhost|host {host}
                         -dbport|port {port} 
                         -dbname {database} 
                         -types {gene,exon,transcript,translation} 
                         -start {first stable ID}

  Argument to -types is a comma-separated list of types of stable IDs to be produced.

  If the -types argument is ommitted, stable IDs are generated for all types (gene,transcript,translation,exon).

  Assigns stable IDs to objects that currently have none. Starting stable ID is found by incrementing the highest 
  current stable ID for that type *or* by using -start argument. The stable_ids are written to STDOUT. If no -start
  option is used the script tries to find the max. given stable_id for each object by looking up the <OBJ>_stable_id
  tables in the database and the gene_archive table ( only  for gene,translation and transcript, not for exon-stable-ids!)


  Note : 

  -start option requires to not submit an initial stable_id without any Gene/Transcript/Exon/Translation ending, like 
   ENSMUS000001 ( not ENSMUSG0001 than you end up with stable-ids like ENSMUSGG001 ENSMUSGT0001...) !

  Again,
  the parameter to -start should be the stable ID you wish to start from without the Gene/Transcript/Translation/Exon identifier

  Examples :  

    - to generate only Exon-stable-ids starting with 223 for Mouse, use 

                    -start ENSMUS222 -types exon 


    - to generate Exon-and Gene-stable-ids starting with 223 for Mouse, use 

                    -start ENSMUS222 -types exon,gene

    - to generate exon,transcript,translation and gene-stable-ids for Human starting, which all start with ID 666 use 

                    -start ENS665     

    - to generate a whole new set of stable_ids ( exon,transcript, translation, gene ) starting with 1 for an organism with 
      prefix ENSNEW you can use one of the following options : 
             
                    -start ENSNEW0          <or> 
                    -start ENSNEW0          <or> 
                    -start ENSNEW00000000
             


  Produces SQL which can be run against the target database.

EOF

}

# --------------------------------------------------------------------------------
