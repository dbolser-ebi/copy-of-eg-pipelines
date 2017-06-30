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
    'compara'          => 'multi',
    'skip_dumps'       => [],
    'files_per_subdir' => 500,
    'mlss_ids'         => [],
    'method_links'     => ['LASTZ_NET', 'TRANSLATED_BLAT_NET'],
  };
}

sub write_output {
  my ($self) = @_;
  my $compara          = $self->param_required('compara');
  my $dump_types       = $self->param_required('dump_types');
  my $compara_dumps    = $self->param_required('compara_dumps');
  my $skip_dumps       = $self->param_required('skip_dumps');
  my $results_dir      = $self->param_required('results_dir');
  my $files_per_subdir = $self->param_required('files_per_subdir');
  my $release_date     = $self->param_required('release_date');
  my $ensembl_release  = $self->param_required('ensembl_release');
  my $mlss_ids         = $self->param_required('mlss_ids');
  my $method_links     = $self->param_required('method_links');
  
  my %skip_dumps = map { $_ => 1 } @$skip_dumps;
  
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($compara, 'compara');
  
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
  
  
  # Check first_release of mlss, if < current release, don't export.
  
  foreach my $flow (keys %$dump_types) {
    foreach my $dump_type (@{$$dump_types{$flow}}) {
      if (!exists $skip_dumps{$dump_type}) {
        my $prefix = $$compara_dumps{$dump_type}."_$release_date";
  
        foreach my $mlss (@mlss) {
          my $mlss_id = $mlss->dbID;
          
          if ($mlss->first_release < $ensembl_release) {
            $self->warning("WGA with MLSS ID $mlss_id should already exist, so dumps will not be generated.");
          } else {
            my ($ref, $sub_dir) = $self->sub_dir($prefix, $mlss, \%assemblies);
            my $aln_dir         = catdir($results_dir, $sub_dir);
            path($aln_dir)->mkpath;
            
            my %output_ids = (
              'mlss_id'     => $mlss_id,
              'ref_species' => $ref,
              'aln_dir'     => $aln_dir,
              'out_dir'     => $results_dir,
              'sub_dir'     => $sub_dir,
            );
            
            $self->dataflow_output_id(\%output_ids, $flow);
          }
        }
      }
    }
  }
}

sub sub_dir {
  my ($self, $prefix, $mlss, $assemblies) = @_;
  
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
  
  return ($ref_species, "$prefix\_$ref\_$non_ref\_$ml_name");
}

1;
