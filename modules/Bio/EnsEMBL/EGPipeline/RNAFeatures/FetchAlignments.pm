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

Bio::EnsEMBL::EGPipeline::RNAFeatures::FetchAlignments

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::FetchAlignments;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my $self = shift @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'db_type' => 'core',
  };
}

sub run {
  my ($self) = @_;
  my $species       = $self->param_required('species');
  my $run_cmscan    = $self->param_required('run_cmscan');
  my $run_trnascan  = $self->param_required('run_trnascan');
  my $alignment_dir = $self->param_required('alignment_dir');
  
  my $dba = $self->get_DBAdaptor($self->param_required('db_type'));
  my $dbh = $dba->dbc->db_handle();
  
  if ($run_cmscan) {
    my $rfam_logic_name   = $self->param_required('rfam_logic_name');
    my $cmscan_cm_file    = $self->param_required('cmscan_cm_file');
    my $cmscan_logic_name = $self->param_required('cmscan_logic_name');
    
    my $logic_name;
    if (! exists $$cmscan_cm_file{$species} && ! exists $$cmscan_cm_file{'all'}) {
      $logic_name = $rfam_logic_name;
    } elsif (exists $$cmscan_logic_name{$species}) {
      $logic_name = $$cmscan_logic_name{$species};
    } elsif (exists $$cmscan_logic_name{'all'}) {
      $logic_name = $$cmscan_logic_name{'all'};
    } else {
      $logic_name = 'cmscan_custom';
    }
    
    my $cmscan_file = "$alignment_dir/cmscan.txt";
    open(my $fh, '>', $cmscan_file) or $self->throw("Failed to open $cmscan_file: $!");
    $self->fetch_alignments($dbh, $logic_name, 'evalue', $fh);
    close($fh);
  }
  
  if ($run_trnascan) {    
    my $logic_name = $self->param_required('trnascan_logic_name');
    
    my $trnascan_file = "$alignment_dir/trnascan.txt";
    open(my $fh, '>', $trnascan_file) or $self->throw("Failed to open $trnascan_file: $!");
    $self->fetch_alignments($dbh, $logic_name, 'score', $fh);
    close($fh);
  }
}

sub fetch_alignments {
  my ($self, $dbh, $logic_name, $value_column, $fh) = @_;
  
  my $sql = "
    SELECT
      $value_column,
      meta_value AS species,
      LEFT(
        SUBSTRING(
          external_data,
          INSTR(
            external_data,
            \"'Biotype' => \"
          )+14
        ),
        INSTR(
          SUBSTRING(
            external_data,
            INSTR(
              external_data,
            \"'Biotype' => \"
            )+14
          ),
          \"'\"
        )-1
      ) AS biotype,
      hit_name
    FROM meta, analysis
    INNER JOIN dna_align_feature USING (analysis_id)
    WHERE logic_name = ?
    AND meta_key = 'species.display_name'
    ORDER BY species, biotype, hit_name
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  while (my @row = $sth->fetchrow_array) {
    print $fh join("\t", @row), "\n";
  }
}

1;
