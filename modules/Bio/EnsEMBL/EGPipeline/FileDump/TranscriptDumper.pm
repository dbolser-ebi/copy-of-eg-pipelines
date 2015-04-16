=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::FileDump::TranscriptDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Bio::EnsEMBL::Utils::IO::FASTASerializer;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'data_type' => 'transcripts',
    'file_type' => 'fa',
  };
}

sub run {
  my ($self) = @_;
  my $species   = $self->param_required('species');
  my $db_type   = $self->param_required('db_type');
  my $out_fh    = $self->param_required('out_fh');
  my $data_type = $self->param_required('data_type');
  
  my $reg = 'Bio::EnsEMBL::Registry';
  
  my $ta = $reg->get_adaptor($species, $db_type, 'Transcript');
  
  my $transcripts = $ta->fetch_all();
  
  my $serialiser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new( $out_fh );

  for my $t ( @$transcripts ){

#easiest to make the change to the header directly in the underlying Bio::PrimarySeqI object

      my $obj = undef ;

      if( $data_type eq 'transcripts' ){ $obj = $t->seq() }
      elsif( $data_type eq 'peptides' ){ $obj = $t->translate }
      else{ $obj = $self->throw( "unsupported -data_type option") }

#NB some transcripts have no translations, so you need to allow for that

      if ( $obj ){

	  my $id = undef;
	  my $start = undef;
	  my $end = undef;

	  if( $data_type eq 'transcripts' ){ 
	      $id = $t->stable_id();
	      $start = $t->seq_region_start();
	      $end = $t->seq_region_end();
	  }
	  elsif( $data_type eq 'peptides' ){ 
	      $id = $t->translation->stable_id() ;
	      $start = $t->translation->genomic_start();
	      $end = $t->translation->genomic_end();
	  }

	  $obj->display_id( rename_header( $id, $t, $start, $end ) );
	  $serialiser->print_Seq( $obj );
      }
  }

  close($out_fh);
 
}

sub rename_header{
#generate VectorBase format fasta header

    my $id = shift @_;
    my $t = shift @_;
    my $start = shift @_;
    my $end = shift @_;

    my $descr = $t->get_Gene->description()  ? $t->get_Gene->description() : 'no description' ;
    $descr =~ s/\s\[Source.+$//;


    my $header = join('|', 
		      "$id $descr" , 
		      $t->biotype() , 
		      join(":" , $t->seq_region_name() , $start ) . '-' . join(':' , $end , $t->strand() ),
		      'gene:' . $t->get_Gene->stable_id(),
	); 

    return $header;
}

1;
