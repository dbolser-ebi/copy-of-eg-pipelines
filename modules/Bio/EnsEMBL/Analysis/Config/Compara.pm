# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Config::Compara

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::Compara;
    use Bio::EnsEMBL::Analysis::Config::Compara qw();

=head1 DESCRIPTION

Compara pipeline configuration.

It imports and sets a number of global variables, mainly paths, into the
calling package. Without arguments all the variables are set, and with a
list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Config> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Config> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Analysis::Config::Compara;

use strict;
use vars qw(%Config);

my $installation = '/nfs/panda/ensemblgenomes/production/compara/binaries';
my $java_home = '/sw/arch/pkg/jdk1.5';

%Config = (
    PYTHON    => 'python',
    JAVA      => "${java_home}/bin/java",
    EXONERATE => "${installation}/exonerate",
    SEMPHY    => "${installation}/semphy",
    ORTHEUS   => "${installation}/ortheus/Ortheus.py",
    TREEBEST  => "${installation}/treebest",
);

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else all
    my @vars = @_ ? @_ : keys(%Config);
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if (defined $Config{ $_ }) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Config{ $_ };
	} else {
	    die "Error: Config: $_ not known\n";
	}
    }
}

1;
