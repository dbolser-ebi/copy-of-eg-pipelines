package Bio::EnsEMBL::Pipeline::Config::General;

use strict;
use vars qw(%Config);
%Config = (
  BIN_DIR  => '/usr/local/ensembl/bin',
  DATA_DIR => '/usr/local/ensembl/data',
  LIB_DIR  => '/usr/local/ensembl/lib',
  PIPELINE_WORK_DIR   => '/tmp',
#  PIPELINE_INPUT_DIR => '/nfs/acari/jb16/projects/Anasp/chunks',
#  PIPELINE_TARGET_DIR => '/data/blastdb/Ensembl/Large/Gallus_gallus/WASHUC1/softmasked_dusted/Gallus_gallus.WASHUC1.softmasked.fa',
  SLICE_INPUT_ID_REGEX => '(\S+)\.(\d+)-(\d+):?([^:]*)',
  PIPELINE_REPEAT_MASKING => ['RepeatMask'],	
  SNAP_MASKING => [],
  MAX_JOB_TIME => 86400,
  KILLED_INPUT_IDS => '',
  RENAME_ON_RETRY => 1,
);



  sub import {
    my ($callpack) = caller(0);
    my $pack = shift;
    my @vars = @_ ? @_ : keys(%Config);
    return unless @vars;
    eval "package $callpack; use vars qw(".
      join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;
    foreach (@vars) {
      if (defined $Config{ $_ }) {
        no strict 'refs';
	*{"${callpack}::$_"} = \$Config{ $_ };
      }else {
	die "Error: Config: $_ not known\n";
      }
    }
  }
  1;
