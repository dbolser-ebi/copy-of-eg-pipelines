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

package Bio::EnsEMBL::EGPipeline::FileDump::XrefDumper;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'gene_centric' => 1,
    'data_type'    => 'xrefs',
    'file_type'    => 'tsv',
    'logic_name'   => [],
    'external_db'  => [],
  };
}

sub run {
  my ($self) = @_;
  my $species       = $self->param_required('species');
  my $db_type       = $self->param_required('db_type');
  my $out_file      = $self->param_required('out_file');
  my $logic_names   = $self->param_required('logic_name');
  my $external_dbs  = $self->param_required('external_db');
  
  # The API code accepts '%' as a pattern-matching character, so this
  # is an easy way to get all xrefs, if no external_dbs are given.
  if (scalar @$external_dbs == 0) {
    push @$external_dbs, '%';
  }
  
  my $reg = 'Bio::EnsEMBL::Registry';
  
  my $sa = $reg->get_adaptor($species, $db_type, 'Slice');
  my $slices = $sa->fetch_all('toplevel');
  
  open(my $out_fh, '>', $out_file) or $self->throw("Cannot open file $out_file: $!");
  $self->print_header($out_fh);
  
  my $ga = $reg->get_adaptor($species, $db_type, 'Gene');
  
  my $genes = [];
  foreach my $slice (@$slices) {
    $self->fetch_genes($ga, $logic_names, $slice, $genes);
  }
  $self->print_feature_list($out_fh, $genes, $external_dbs);
  
  close($out_fh);
  my $out_files = $self->param('out_files');
  
  $self->param('out_files', $out_files);
}

sub fetch_genes {
  my ($self, $adaptor, $logic_names, $slice, $genes) = @_;
  
  if (scalar(@$logic_names) == 0) {
    push @$genes, @{$adaptor->fetch_all_by_Slice($slice)};
  } else {
    foreach my $logic_name (@$logic_names) {
      push @$genes, @{$adaptor->fetch_all_by_Slice($slice, $logic_name)};
    }
  }
}

sub print_header {
  my ($self, $fh) = @_;
  
  say $fh join("\t",
    qw(
      gene_stable_id transcript_stable_id protein_stable_id xref
      db_name info_type source_identity xref_identity linkage_type
    )
  );
}

sub print_feature_list {
  my ($self, $fh, $genes, $external_dbs) = @_;
  
  my %xrefs;
  
  foreach my $gene (@$genes) {
    foreach my $external_db (@$external_dbs) {
      my $xref_list = $gene->get_all_DBEntries($external_db);
      push @{$xrefs{$gene->stable_id}{'-'}{'-'}}, @$xref_list;
    }
    
    my $transcripts = $gene->get_all_Transcripts;
    foreach my $transcript (@$transcripts) {
      my $tn_id = defined($transcript->translation) ? $transcript->translation->stable_id : '-';
      
      foreach my $external_db (@$external_dbs) {
        my $xref_list = $transcript->get_all_DBLinks($external_db);
        push @{$xrefs{$gene->stable_id}{$transcript->stable_id}{$tn_id}}, @$xref_list;
      }
    }
  }
        
  foreach my $g_id (sort keys %xrefs) {
    foreach my $tt_id (sort keys %{$xrefs{$g_id}}) {
      foreach my $tn_id (sort keys %{$xrefs{$g_id}{$tt_id}}) {
        foreach my $xref (sort {$a->primary_id cmp $b->primary_id} @{$xrefs{$g_id}{$tt_id}{$tn_id}}) {
          my $xref_id   = $xref->primary_id;
          my $xref_db   = $xref->dbname;
          my $info_type = $xref->info_type;
          
          my $src_identity  ='-';
          my $xref_identity ='-';
          if ($xref->isa('Bio::EnsEMBL::IdentityXref')) { 
            $src_identity  = $xref->ensembl_identity; 
            $xref_identity = $xref->xref_identity; 
          }

          my $linkage_type = '-';
          if ($xref->isa('Bio::EnsEMBL::OntologyXref')) { 
            $linkage_type = join(' ', @{$xref->get_all_linkage_types});
          }

          say $fh join("\t",
            $g_id, $tt_id, $tn_id, $xref_id, $xref_db, $info_type,
            $src_identity, $xref_identity, $linkage_type
          );
        }
      }
    }
  }
}

1;
