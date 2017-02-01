# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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


# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::GenBlast
#
# Copyright (c) 2009 WormBase
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::HiveGenBlast

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::HiveGenBlast->
  new(
      -input_id => 'file_name',
      -db => $db,
      -analysis => $analysis,
     );
  $blast->fetch_input;
  $blast->run;
  $blast->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable GenBlast which wraps the program GenBlast

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the dna_align_feature  tables 

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::HiveGenBlast;

use strict;
use warnings;


use Bio::EnsEMBL::Analysis::Runnable::WormBase::GenBlast;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

#use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisRun');
#use base ('Bio::EnsEMBL::Hive::Process');
use base qw(Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base);

sub run {
  my ($self)=@_;
  my $runnable = $self->fetch_runnable();
  $self->dbc && $self->dbc->disconnect_if_idle();
  $self->dbc->reconnect_when_lost(1);
  $runnable->run();
  
  $self->write_output($runnable->output);
}

=head2 fetch_runnable

  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut

sub fetch_runnable {
  my ($self) = @_;
  my %genome_slices;

  my $dba = $self->get_DBAdaptor('core');
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name('genblast');


  my $genome_file = $analysis->db_file;
  open(my $fh, $genome_file) or throw("Could not open $genome_file for reading");
  while(<$fh>) {
    /^\>(\S+)/ and do {
      my $seq_name = $1;
      my $slice = $dba->get_SliceAdaptor->fetch_by_region('toplevel', $seq_name);
      if (not defined $slice) {
        throw("Could not extract slice for $seq_name from database");
      }
      $genome_slices{$seq_name} = $slice;
    }
  }

  # flag
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::WormBase::GenBlast->new
    (
     -query => $self->param('query'),
     -program => $analysis->program_file,
     -analysis => $analysis,
     -database => $analysis->db_file,
     -refslices => \%genome_slices,
     -options => $analysis->parameters,
    );
  return $runnable;
}


=head2 write_output

  Function  : writes the prediction transcripts back to the database
  after validation
  Returntype: none

=cut



sub write_output{
  my ($self,$output) = @_;

  my $dba = $self->get_DBAdaptor('core');
  my $adaptor = $dba->get_PredictionTranscriptAdaptor;
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name('genblast');

  foreach my $pt(@$output){
    $pt->analysis($analysis);
    $pt->slice($self->query) if(!$pt->slice);
    $adaptor->store($pt);
  }
}

1;
