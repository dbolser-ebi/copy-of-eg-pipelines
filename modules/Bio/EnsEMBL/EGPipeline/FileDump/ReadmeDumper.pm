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

package Bio::EnsEMBL::EGPipeline::FileDump::ReadmeDumper;

use strict;
use warnings;
no  warnings 'redefine';
use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);

sub run {
  my ($self) = @_;
  
  my $out_fh          = $self->param('out_fh');
  my $readme_template = $self->param('readme_template');
  
  my ($division, undef) = $self->get_division();
  $readme_template =~ s/DIVISION/$division/gm;
  
  print $out_fh $readme_template;
  
  close($out_fh);
}

1;
