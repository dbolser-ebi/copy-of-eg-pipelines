=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::LoadGFF3::LoadPatch;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::SeqIO;
use Time::Piece;

sub run {
  my ($self) = @_;
  my $seq_file = $self->param_required('seq_file');
  
  my $dba = $self->get_DBAdaptor('core');
  my $csa = $dba->get_CoordSystemAdaptor();
  my $sa  = $dba->get_SliceAdaptor();
  my $aa  = $dba->get_AttributeAdaptor();
  
  my $contig_cs   = $csa->fetch_sequence_level();
  my $scaffold_cs = $csa->fetch_by_rank($contig_cs->rank() - 1);
  
  open(F, $seq_file);
  my $seq_io = Bio::SeqIO->new(
    -fh     => \*F,
    -format => 'fasta',
  );
  
  while (my $seq_obj = $seq_io->next_seq) {
    my $contig_name   = $seq_obj->id;
    my $scaffold_name = "VB_PATCH_$contig_name"; # Maybe keep the same? Want the alignment on toplevel, really.
    my $length        = $seq_obj->length;
    
    my $contig_seq_region_id =
      $self->fetch_or_create_slice($sa, $contig_cs, $contig_name, $length, \$seq_obj->seq);
    
    my $scaffold_seq_region_id =
      $self->fetch_or_create_slice($sa, $scaffold_cs, $scaffold_name, $length);
    
    $self->add_slice_mapping($dba, $scaffold_seq_region_id, $contig_seq_region_id, $length);
    
    $self->add_attributes($aa, $scaffold_seq_region_id);
  }
  
  close(F);
}

sub fetch_or_create_slice {
  my ($self, $sa, $coordinate_system, $name, $length, $sequence) = @_;
  
  my $seq_region_id;
  my $slice = $sa->fetch_by_region($coordinate_system->name, $name);
  if (defined $slice) {
    $seq_region_id = $slice->get_seq_region_id;
  } else {
    $slice = $self->create_slice($coordinate_system, $name, $length);
    $seq_region_id = $sa->store($slice, $sequence);
  }
  
  return $seq_region_id;
}

sub create_slice {
  my ($self, $coordinate_system, $name, $length) = @_;
  
  my $slice = Bio::EnsEMBL::Slice->new(
    -seq_region_name   => $name,
    -start             => 1,
    -end               => $length,
    -seq_region_length => $length,
    -strand            => 1,
    -coord_system      => $coordinate_system,
  );
  
  return $slice;
}

sub add_slice_mapping {
  my ($self, $dba, $asm_seq_region_id, $cmp_seq_region_id, $length) = @_;
  
  my $sql =
    "INSERT IGNORE INTO assembly (
      asm_seq_region_id,
      asm_start,
      asm_end,
      cmp_seq_region_id,
      cmp_start,
      cmp_end,
      ori
    ) VALUES (?, ?, ?, ?, ?, ?, ?)";
  my $dbh = $dba->dbc->db_handle();
  my $sth = $dbh->prepare($sql);
  
  $sth->execute($asm_seq_region_id, 1, $length, $cmp_seq_region_id, 1, $length, 1)
    or throw("Failed to insert assembly mapping with $sql: ".$sth->errstr);
}

sub add_attributes {
  my ($self, $aa, $seq_region_id) = @_;
  my @attributes;
  
  my $time = localtime;
  my $timestamp = $time->strftime('%Y-%m-%d %H-%M-%S');
  
  push @attributes, $self->create_attribute('toplevel',    1);
  push @attributes, $self->create_attribute('non_ref',     1);
  push @attributes, $self->create_attribute('patch_novel', $timestamp);
  
  $aa->store_on_Object($seq_region_id, \@attributes, 'seq_region');
}

sub create_attribute {
  my ($self, $key, $value) = @_;
  
  my $attribute = Bio::EnsEMBL::Attribute->new(
    -CODE  => $key,
    -VALUE => $value
  );
  
  return $attribute;
}

1;
