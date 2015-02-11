
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

=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::Rfam::RfamLoader

=head1 DESCRIPTION

Module for adding genes from matches to INSDC sequences calculated by Rfam

=head1 Author

Dan Staines

=cut

package Bio::EnsEMBL::EGPipeline::Rfam::RfamLoader;
use Log::Log4perl qw/:easy/;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use base Bio::EnsEMBL::EGPipeline::BaseLoader;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;

=head1 CONSTRUCTOR
=head2 new
  Arg [-RFAM_DBA]  : 
       adaptor for UniProt Oracle database (e.g. SWPREAD)

  Example    : $ldr = Bio::EnsEMBL::EGPipeline::Rfam::RfamLoader->new(...);
  Description: Creates a new loader object
  Returntype : Bio::EnsEMBL::EGPipeline::Rfam::RfamLoader
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub new {
  my ($proto, @args) = @_;
  my $self = $proto->SUPER::new(@args);
  ( $self->{rfam_dba} )
	= rearrange(
				 [ 'RFAM_DBA' ],
				 @args );
  return $self;
}

=head1 METHODS
=head2 add_genes
  Arg        : Bio::EnsEMBL::DBSQL::DBAdaptor for core database to write to
  Description: Add genes to supplied core
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut
sub add_genes {
    my ( $self, $dba ) = @_;
    my $analysis = $self->get_analysis( $dba, 'rfam_genes' );
    # find all top-level slices   
    for my $slice (@{$dba->get_SliceAdaptor->fetch_all('toplevel')}) {
        $self->logger->info("Processing ".$slice->seq_region_name());
        $self->add_genes_to_slice($dba,$analysis,$slice);
    }
    return;
}

sub add_genes_to_slice {
    my ($self,$dba,$analysis,$slice) = @_;
    my $ga = $dba->get_GeneAdaptor();
    my $dba = $dba->get_DBEntryAdaptor();
    # for each slice, find its INSDC accession
    my $insdc = $slice->get_all_synonyms("EMBL");
    my $n = 0;
    if(defined $insdc && scalar(@$insdc)>0) {
        # for a given accession, query the RFAM database for matches
        # for each match, create a gene-transcript pair
        for my $acc (@{$insdc}) {
            $self->logger->info("Finding matches for ".$acc->name());
            $self->{rfam_dba}->dbc()->sql_helper()->execute_no_return(
                -SQL=>q/select rfam_acc,rfam_id,description,family.type,
                seq_start,seq_end,case when seq_start>seq_end then -1 else 1 end as seq_strand 
                from full_region join family using (rfam_acc) 
                where rfamseq_acc=? 
                and family.type like 'Gene%' and is_significant=1
                /,
                -PARAMS=>[$acc->name()],
                -CALLBACK => sub {
                    my ( $acc, $id, $des, $types, $start, $end, $strand ) = @{ $_[0] };
                    if($strand == -1) {
                        my $tmp = $start;
                        $start = $end;
                        $end = $tmp;
                    }
                    my ($g,$type) = split $types, '; ?';
                    $self->logger->debug("Adding $acc/$id to ".$slice->seq_region_name());
                    my $gene = Bio::EnsEMBL::Gene->new(-VERSION       => 1,
                                                       -SLICE         => $slice,
                                                       -BIOTYPE       => $type,
                                                       -SOURCE        => 'Rfam',
                                                       -DESCRIPTION   => $des,
                                                       -ANALYSIS      => $analysis);

                    my $xref = Bio::EnsEMBL::DBEntry->new(-DBNAME     => 'RFAM_GENE',
                                                          -PRIMARY_ID => $id,
                                                          -DISPLAY_ID => $id,
                                                          -ANALYSIS   => $analysis);
                    $dba->store($xref);
                    $gene->display_xref($xref);
                    $gene->add_DBEntry(Bio::EnsEMBL::DBEntry->new(-DBNAME     => 'RFAM',
                                                                  -PRIMARY_ID => $acc,
                                                                  -DISPLAY_ID => $acc,
                                                                  -DESCRIPTION=> $des,
                                                                  -ANALYSIS   => $analysis));

                    
                    my $transcript = Bio::EnsEMBL::Transcript->new(-VERSION       => 1,
                                                                   -BIOTYPE       => $type,
                                                                   -SLICE         => $slice,
                                                                   -ANALYSIS      => $analysis);

                    $transcript->add_Exon(Bio::EnsEMBL::Exon->new(-VERSION       => 1,
                                                                  -SLICE         => $slice,
                                                                  -START         => $start,
                                                                  -END           => $end,
                                                                  -STRAND        => $strand,
                                                                  -PHASE         => 0,
                                                                  -END_PHASE     => 0));
                    $self->logger()->debug("Adding transcript to $acc/$id");
                    $gene->add_Transcript($transcript);
                    $self->logger()->debug("Storing $acc/$id");
                    $ga->store($gene);
                    $n++;
                    return;
                }                
                );
        }
    }
    $self->logger()->info("Added $n genes to ".$slice->seq_region_name());
    return;
}

1;
