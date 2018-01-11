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

use base ('Bio::EnsEMBL::VEP::Pipeline::DumpVEP::QCDump');

use File::Path qw(make_path rmtree);

# The QC module makes assumptions about what directories have been
# created by earlier modules, which won't be true for VB.
sub dump_dir {
  my $self = shift;
  return $self->param('pipeline_dir').$self->param('dir_suffix');
}

# Override bits of the qc function that assume you have a tarball of the
# dumps which you need to unpack to test; we have them sitting in a
# directory already, so just need to point there. (The function makes
# some assumptions about job input_ids which are violated in the case
# of VB, so we can't just let it do the tarballing and get on with it...)
sub qc {
  my ($self, $mod) = @_;

  my $type      = $self->param('type');
  my $has_var   = $self->param('variation');
  my $has_reg   = $self->param('regulation');
  my $converted = $mod && $mod =~ /tabix/;
  my $species   = $self->required_param('species');
  my $assembly  = $self->required_param('assembly');
  
  my $method_name = ($type eq 'core' ? '' : $type.'_').'species_suffix';
  my $source_dir  = $self->dump_dir.'/'.$self->$method_name;
  
  die("ERROR: Expected to find $source_dir\n") unless -d $source_dir;
  die("ERROR: Expected to find $source_dir/info.txt") unless -e "$source_dir/info.txt";
  
  my $qc_dir = $self->dump_dir."/qc/$species/$assembly";
  unless(-d $qc_dir) {
    make_path($qc_dir) or die "Could not make directory $qc_dir";
  }
  
  my $config_obj = Bio::EnsEMBL::VEP::Config->new({
    dir => $source_dir,
    offline => 1,
    species => $species,
    assembly => $assembly,
    check_existing => 1,
    regulatory => 1,
  });
  my $cache_dir_obj = Bio::EnsEMBL::VEP::CacheDir->new({dir => $source_dir, config => $config_obj});

  $self->check_info($cache_dir_obj, $converted);

  $self->check_annotation_sources($cache_dir_obj, $converted);

  $self->check_dirs($cache_dir_obj, $converted);

  $self->run_test_set($self->dump_dir) if $has_var && $type ne 'refseq';

  rmtree($qc_dir);
}

1;
