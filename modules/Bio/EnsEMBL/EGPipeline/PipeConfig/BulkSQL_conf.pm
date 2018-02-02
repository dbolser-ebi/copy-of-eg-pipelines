=head1 LICENSE

Copyright [2009-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::EGPipeline::PipeConfig::BulkSQL_conf

=head1 DESCRIPTION

Pipeline to run SQL commands across a set of databases.

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::BulkSQL_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::PipeConfig::Base_conf');

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'bulk_sql',
    
    species      => [],
    division     => [],
    antispecies  => [],
    run_all      => 0,
    meta_filters => {},
    
    sql_file => undef,
  };
}

sub pipeline_analyses {
  my ($self) = @_;
  
  return [
    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
      -max_retry_count => 0,
      -input_ids       => [ {} ],
      -parameters      => {
                           species      => $self->o('species'),
                           antispecies  => $self->o('antispecies'),
                           division     => $self->o('division'),
                           run_all      => $self->o('run_all'),
                           meta_filters => $self->o('meta_filters'),
                          },
      -flow_into       => {
                           '2' => 'SqlExecute',
                          },
      -rc_name         => 'default',
    },

    {
      -logic_name      => 'SqlExecute',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::SqlExecute',
      -parameters      => {
                            sql_file => $self->o('sql_file'),
                          },
      -max_retry_count => 1,
      -hive_capacity   => 10,
      -rc_name         => 'default',
    },

  ];
}

1;
