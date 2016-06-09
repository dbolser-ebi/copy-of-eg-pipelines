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

package Bio::EnsEMBL::EGPipeline::FileDump::VEPDumper;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use File::Copy;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'data_type' => 'vep',
    'file_type' => 'vep',
  };
}

sub run {
  my ($self) = @_;
  my $perl_lib   = $self->param_required('perl_lib');
  my $vep_script = $self->param_required('vep_script');
  my $vep_params = $self->param_required('vep_params');
  my $species    = $self->param_required('species');
  my $dir        = $self->param_required('results_dir');
  
  my $mc       = $self->core_dba->get_MetaContainer();
  my $assembly = $mc->single_value_by_key('assembly.default');
  my $host     = $mc->dbc->host();
  my $port     = $mc->dbc->port();
  
  $vep_params .= " --species $species";
  $vep_params .= " --assembly $assembly";
  $vep_params .= " --dir $dir";
  $vep_params .= " --host $host";
  $vep_params .= " --port $port";
  $vep_params .= " --user ensro";
  
	my $vep_cmd = "perl -I $perl_lib $vep_script $vep_params";
  
  my $finished = 0;
  
  open CMD, "$vep_cmd 2>&1 |" or $self->throw("Failed to run command: $vep_cmd");
  my @buffer;
  while(<CMD>) {
    $finished = 1 if /Finished/;
    push @buffer, $_;
    shift @buffer if scalar @buffer > 5;
  }
  close CMD;
  
  $self->throw("Error running VEP command: $vep_cmd\n".join("", @buffer)."\n") unless $finished;
}

sub write_output {
  my ($self) = @_;
  my $species  = $self->param_required('species');
  my $dir      = $self->param_required('results_dir');
  my $out_file = $self->param_required('out_file');
  my $vep_dir  = "$dir/$species";
  
  # Reinstate the upper case letter (for alternative strains)
  # at the end of species names.
  my $lc_species = lc($species);
  if ($species ne $lc_species) {
    if (-e "$dir/$lc_species") {
      move("$dir/$lc_species", "$dir/$species");
    }
  }
  
  $self->dataflow_output_id({
    out_file => $out_file,
    vep_dir  => $vep_dir,
  }, 1);
}

1;
