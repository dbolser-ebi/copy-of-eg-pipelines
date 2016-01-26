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

package Bio::EnsEMBL::EGPipeline::FileDump::PairwiseAlignmentFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Registry;
use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'skip_dumps'       => [],
    'files_per_subdir' => 500,
    'mlss_ids'         => [],
    'method_links'     => ['LASTZ_NET', 'TRANSLATED_BLAT_NET'],
  };
}

sub write_output {
  my ($self) = @_;
  my $dump_types       = $self->param_required('dump_types');
  my $skip_dumps       = $self->param_required('skip_dumps');
  my $results_dir      = $self->param_required('results_dir');
  my $files_per_subdir = $self->param_required('files_per_subdir');
  my $mlss_ids         = $self->param_required('mlss_ids');
  my $method_links     = $self->param_required('method_links');
  
  my %skip_dumps = map { $_ => 1 } @$skip_dumps;
  
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Multi', 'compara');
  
  my $mlssa = $dba->get_adaptor("MethodLinkSpeciesSet");
  
  my $gdba = $dba->get_adaptor("GenomeDB");
  my $gdbs = $gdba->fetch_all_current;
  my %assemblies;
  foreach my $gdb (@$gdbs) {
    $assemblies{$gdb->name} = $gdb->assembly;
  }
  
  my @mlss;
  if (scalar @$mlss_ids) {
    @mlss = @{ $mlssa->fetch_all_by_dbID_list($mlss_ids) };
  } else {
    foreach my $method_link (@$method_links) {
      push @mlss, @{ $mlssa->fetch_all_by_method_link_type($method_link) };
    }
  }
  
  foreach my $flow (keys %$dump_types) {
    foreach my $dump_type (@{$$dump_types{$flow}}) {
      if (!exists $skip_dumps{$dump_type}) {        
        foreach my $mlss (@mlss) {
          my $mlss_id       = $mlss->dbID;
          my ($ref, $label) = $self->generate_label($mlss, \%assemblies);
          my $out_dir       = catdir($results_dir, $dump_type, $label);
          path($out_dir)->mkpath;
          
          my %output_ids = (
            'mlss_id'     => $mlss_id,
            'ref_species' => $ref,
            'out_dir'     => $out_dir,
          );
          
          $self->dataflow_output_id(\%output_ids, $flow);
        }
      }
    }
  }
}

sub generate_label {
  my ($self, $mlss, $assemblies) = @_;
  
  my $ref     = $mlss->get_tagvalue('reference_species');
  my $non_ref = $mlss->get_tagvalue('non_reference_species');
  
  my $ref_species = $ref;
  
  my $ref_assembly     = $$assemblies{$ref};
  my $non_ref_assembly = $$assemblies{$non_ref};
  
  $ref     =~ s/_/-/g;
  $non_ref =~ s/_/-/g;
  
  $ref     = ucfirst($ref)    ."-$ref_assembly";
  $non_ref = ucfirst($non_ref)."-$non_ref_assembly";
  
  my $ml = $mlss->method;
  my $ml_name = lc($ml->type);
  $ml_name =~ s/_net$//;
  $ml_name =~ s/translated_blat/tblat/;
  
  return ($ref_species, "$ref\_$non_ref\_$ml_name");
}

1;
