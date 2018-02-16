=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

Bio::EnsEMBL::EGPipeline::PipeConfig::GeneTreeHighlighting_conf

=head1 DESCRIPTION

Populate compara table with GO and InterPro terms,
to enable highlighting in the genome browser.

=head1 AUTHOR

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::GeneTreeHighlighting_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::Version 2.4;

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    pipeline_name => 'gene_tree_highlighting_' . $self->o('ensembl_release'),

    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},

    production_db_host => 'mysql-eg-pan-prod.ebi.ac.uk',
    production_db_port => 4276,
    production_db_user => 'ensro',
    production_db_name => 'ensembl_production',

    highlighting_capacity => 10,
    compara_division      => 'vb',
  }
}

sub beekeeper_extra_cmdline_options {
  my ($self) = @_;

  my $options = join(' ',
    $self->SUPER::beekeeper_extra_cmdline_options,
    "-reg_conf ".$self->o('registry')
  );

  return $options;
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name      => 'CreateExternalDB',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -input_ids       => [ {} ],
      -parameters      =>
                          {
                            cmd =>
                            'mysqldump '.
                              ' -h '.$self->o('production_db_host').
                              ' -P '.$self->o('production_db_port').
                              ' -u '.$self->o('production_db_user').
                              '    '.$self->o('production_db_name').
                              ' master_external_db '.
                            ' | '.
                            "sed -e 's/master_external_db/external_db/g' ".
                            ' | '.
                            'mysql '.
                              ' -h '.$self->o('compara_db_host').
                              ' -P '.$self->o('compara_db_port').
                              ' -u '.$self->o('compara_db_user').
                              ' -p'. $self->o('compara_db_pass').
                              '    '.$self->o('compara_db_name')
                          },
      -max_retry_count => 1,
      -flow_into       => ['SpeciesFactory'],
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            regulation_flow => 0,
                            variation_flow  => 0,
                          },
      -max_retry_count => 1,
      -flow_into       => {
                            '2->A' => ['HighlightGO'],
                            'A->2' => ['HighlightInterPro'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'HighlightGO',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::GeneTreeHighlight::HighlightGO',
      -hive_capacity   => $self->o('highlighting_capacity'),
      -parameters      => {
                            compara_division => $self->o('compara_division'),
                          },
      -max_retry_count => 1,
      -rc_name 	       => 'normal',
    },

    {
      -logic_name      => 'HighlightInterPro',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::GeneTreeHighlight::HighlightInterPro',
      -hive_capacity   => $self->o('highlighting_capacity'),
      -parameters      => {
                            compara_division => $self->o('compara_division'),
                          },
      -max_retry_count => 1,
      -rc_name 	       => 'normal',
    },
  ];
}

1;
