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

package Bio::EnsEMBL::EGPipeline::FileDump::ValidateVEP;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use Path::Tiny qw(path);
use TAP::Harness;

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'data_type' => 'QC-REPORT',
    'file_type' => 'txt',
  };
}

sub run {
  my ($self) = @_;
  my $vep_hc_script = $self->param_required('vep_hc_script');
  my $vep_hc_params = $self->param_required('vep_hc_params');
  my $species       = $self->param_required('species');
  my $dir           = $self->param_required('results_dir');
  my $out_file      = $self->param_required('out_file');
  
  my $dbc  = $self->core_dbc();
  my $host = $dbc->host();
  my $port = $dbc->port();
  
  $vep_hc_params .= " --species $species";
  $vep_hc_params .= " --dir $dir";
  $vep_hc_params .= " --host $host";
  $vep_hc_params .= " --port $port";
  $vep_hc_params .= " --user ensro";
  
  $self->warning("Running command $vep_hc_script $vep_hc_params");
  
  my @test_args = split(/\s+/, $vep_hc_params);
  
  my $fh = path($out_file)->filehandle('>');
  
  my $harness = TAP::Harness->new({
    'test_args' => \@test_args,
    'stdout'    => $fh,
  });
  $harness->runtests($vep_hc_script);
  
  close $fh;
  
  my $data = path($out_file)->slurp;
  if ($data !~ /Result:\s+PASS/m) {
    die "Validation failed, details in: $out_file";
  }
}

1;
