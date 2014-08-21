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

=cut


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::ProteinSimilarity::CreateBlastDB

=head1 DESCRIPTION

Convert a fasta file into a BLAST database.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::ProteinSimilarity::CreateBlastDB;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Basename qw(fileparse);

sub param_defaults {
  return {
    'blast_type'    => 'ncbi',
    'blast_db'      => undef,
    'blast_db_type' => 'prot',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $makeblastdb_exe = $self->param_required('makeblastdb_exe');
  my $db_fasta_file = $self->param_required('db_fasta_file');
  
  if (!-e $makeblastdb_exe) {
    $self->throw("makeblastdb executable '$makeblastdb_exe' does not exist.");
  }
  if (!-e $db_fasta_file) {
    $self->throw("Fasta file '$db_fasta_file' does not exist.");
  }
  
  if (!$self->param_is_defined('blast_db')) {
    my ($name, $dirs, undef) = fileparse($db_fasta_file, qr/\.[^.]*/);
    $self->param('blast_db', "$dirs$name.blastdb");
  }
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({'blast_db' => $self->param('blast_db')}, 1);
}

sub run {
  my ($self) = @_;
  my $makeblastdb_exe = $self->param('makeblastdb_exe');
  my $db_fasta_file   = $self->param('db_fasta_file');
  my $blast_db        = $self->param('blast_db');
  my $blast_type      = $self->param('blast_type');
  my $blast_db_type   = $self->param('blast_db_type');
  
  #  #formatdb -i fly.fasta -p F -n fly for blast_type=wu dna, -p T for proteins
  my $cmd;
  if ($blast_type eq 'wu') {
    $cmd = "$makeblastdb_exe -i $db_fasta_file -n $blast_db";
    if ($blast_db_type eq 'prot') {
      $cmd .= " -p T";
    } else {
      $cmd .= " -p F";
    }
  } else {
    $cmd = "$makeblastdb_exe -in $db_fasta_file -out $blast_db -dbtype $blast_db_type";
  }
  
  my $out = `$cmd 2>&1`;
  if ($out =~ /error/mi) {
    $self->throw("Error when executing $cmd:\n$out");
  } else {
    if ($blast_type eq 'wu') {
      if (!(-e "$blast_db.nhr" && -e "$blast_db.nin" && -e "$blast_db.nsq")) {
        $self->throw("Error when executing $cmd:\n$out");
      }
    } else {
      if (!(-e "$blast_db.phr" && -e "$blast_db.pin" && -e "$blast_db.psq")) {
        $self->throw("Error when executing $cmd:\n$out");
      }
    }
  }
  
  # To keep the Ensembl Blast runnable happy, need to do this...
  `touch $blast_db` unless -e $blast_db;
}

1;
