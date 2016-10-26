=head1 LICENSE

Copyright [1999-2014] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

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


=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::PipeConfig::GenBlast_conf

=head1 DESCRIPTION

Align amino acids sequences against a set of core databases.

=head1 Author

XX

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::GenBlast_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::Version 2.4;

use File::Spec::Functions qw(catdir);

sub default_options {
  my $self = shift @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'genblastg_protein_'.$self->o('ensembl_release'),
    
    species      => [],
    antispecies  => [],
    division     => [],
    run_all      => 0,
    meta_filters => {},
    max_hive_capacity => 150,
    
    source_species      => [],
    source_antispecies  => [],
    source_division     => [],
    source_run_all      => 0,
    source_meta_filters => {},
    blast_threads       => 1,

  };
}

sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  
  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );
  
  return $options;
}

sub hive_meta_table {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_analyses {
  my $self = shift @_;
  
  return [
    {
      -logic_name      => 'InitialisePipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -input_ids       => [{}],
      -flow_into       => {
                            '1' => ['SpeciesFactory'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count => 1,
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            variation_flow  => 0,
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['GenBlastFactory'],
                          },
      -meadow_type     => 'LOCAL',
    },

    { # read in the protein chunks and fan out the analysis
      -logic_name      => 'GenBlastFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::GenBlastFactory',
      -max_retry_count => 1,
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['GenBlast'],
                          },
    },

    { # -query => protein_chunk / -database => genome_file |db_file is the respective species genome file
      -logic_name      => 'GenBlast',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::HiveGenBlast',
      -can_be_empty    => 1,
      -parameters      => { 
                             logic_name => 'genblast',
                             queryfile =>  '#query#',
                             save_object_type => 'PredictionTranscript',
                          },
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 3,
      -rc_name         => '8GB_threads',
      -flow_into       => {'-1' => 'GenBlast_HighMem'},
    },

    { # -query => protein_chunk / -database => genome_file |db_file is the respective species genome file
      -logic_name      => 'GenBlast_HighMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::HiveGenBlast',
      -can_be_empty    => 1,
      -parameters      => { 
                             logic_name => 'genblast',
                             queryfile =>  '#query#',
                             save_object_type => 'PredictionTranscript',
                          },
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 3,
      -rc_name         => '16GB_threads',
      -flow_into       => {'-1' => 'GenBlast_HigherMem'},

    },
    { # -query => protein_chunk / -database => genome_file |db_file is the respective species genome file
      -logic_name      => 'GenBlast_HigherMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::HiveGenBlast',
      -can_be_empty    => 1,
      -parameters      => { 
                             logic_name => 'genblast',
                             queryfile =>  '#query#',
                             save_object_type => 'PredictionTranscript',
                          },
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 3,
      -rc_name         => '32GB_threads',
    },

  ];
}

sub resource_classes {
  my ($self) = @_;
  
  my $blast_threads = $self->o('blast_threads');

  return {
    %{$self->SUPER::resource_classes},
    '8GB_threads' => {'LSF' => '-q production-rh7 -n ' . ($blast_threads + 1) . ' -M 8000 -R "rusage[mem=8000,tmp=8000]"'},
    '16GB_threads' => {'LSF' => '-q production-rh7 -n ' . ($blast_threads + 1) . ' -M 16000 -R "rusage[mem=16000,tmp=16000]"'},
    '32GB_threads' => {'LSF' => '-q production-rh7 -n ' . ($blast_threads + 1) . ' -M 32000 -R "rusage[mem=32000,tmp=32000]"'},
  }
}

1;
