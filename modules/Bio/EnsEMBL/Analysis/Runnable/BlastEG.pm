=head1 LICENSE

# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::BlastEG

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastEG->
  new(
      -query    => $slice,
      -program  => 'blastn',
      -database => 'uniprot_sprot_invertebrates',
      -options  => 'hitdist=40',
      -parser   => $parser,
      -filter   => undef,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

  This module is a wrapper for running blast. It defaults to run
  ncbi-blast, but can run wu-blast. It requires a parser object,
  e.g. FilterBPlite.

=cut

package Bio::EnsEMBL::Analysis::Runnable::BlastEG;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($type, $database, $parser) =
    rearrange(['TYPE', 'DATABASE', 'PARSER'], @args);
  
  $type = 'ncbi' unless defined $type;
  throw("DATABASE option is mandatory") unless defined $database;
  throw("PARSER option is mandatory") unless defined $parser;
  
  $self->type($type);
  $self->database($database);
  $self->parser($parser);
  
  return $self;
}

sub type {
  my $self = shift;
  $self->{'type'} = shift if(@_);
  return $self->{'type'};
}

sub database {
  my $self = shift;
  $self->{'database'} = shift if(@_);
  return $self->{'database'};
}

sub parser {
  my $self = shift;
  $self->{'parser'} = shift if(@_);
  return $self->{'parser'};
}

sub run_analysis {
  my ($self) = @_;
  my $command     = $self->program;
  my $database    = $self->database;
  my $queryfile   = $self->queryfile;
  my $resultsfile = $self->resultsfile;

  if ($self->type eq 'wu') {
    if (! -e $ENV{BLASTMAT} && ! -e $ENV{WUBLASTMAT}) {  
      throw("Environment variable \$BLASTMAT is not set !!! ");
    }
    $command .= " $database $queryfile -o $resultsfile ";
  } else {
    $command .= " -db $database -query $queryfile -out $resultsfile ";
  }
  $command .= $self->options;
    
  print "Running blast $command\n\n";
  
  system($command) == 0 or throw("FAILED to run ".$command);
}

sub parse_results {
  my ($self, $resultsfile) = @_;
  $resultsfile = $self->resultsfile unless defined $resultsfile;
  
  # Remove NCBI footer as it confuses the Ensembl parser.
  {
    local $/ = undef;
    open(FILE, $resultsfile) or throw("FAILED to open $resultsfile: $!");
    my $results = <FILE>;
    close FILE;
    
    $results =~ s/^\s*Lambda\s+K\s+H\s+a\s+alpha.+\z//ms;
    
    open(FILE, ">$resultsfile") or throw("FAILED to open $resultsfile: $!");
    print FILE $results;
    close FILE;
  }
  
  $self->output($self->parser->parse_files([$resultsfile]));
}

1;
