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

Bio::EnsEMBL::EGPipeline::PipeConfig::BlastProtein_conf

=head1 DESCRIPTION

Align amino acids sequences against a set of core databases.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::WormBlast2_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');

use Bio::EnsEMBL::Hive::Version 2.4;

use File::Spec::Functions qw(catdir);

sub default_options {
  my $self = shift @_;
  return {
    %{$self->SUPER::default_options},
    
    pipeline_name => 'blast_protein_'.$self->o('ensembl_release'),
    
    species      => [],
    antispecies  => [],
    division     => [],
    run_all      => 0,
    meta_filters => {},
    
    # Parameters for dumping and splitting Fasta query files.
    max_seq_length          => 250000,
    max_seq_length_per_file => $self->o('max_seq_length'),
    max_seqs_per_file       => 500,
    max_files_per_directory => 500,
    max_dirs_per_directory  => $self->o('max_files_per_directory'),
    
    max_hive_capacity => 100,
    
    # Source for proteomes can be one of: 'file', 'database', 'uniprot'
    proteome_source    => 'file',
    source_label       => 'Protein alignments',
    is_canonical       => 1,
    
    # This parameter should only be specified if proteome_source = 'file'
    db_fasta_file => undef,
    
    # Dump the proteomes of other species and align them.
    source_species      => [],
    source_antispecies  => [],
    source_division     => [],
    source_run_all      => 0,
    source_meta_filters => {},
    
    # Default is to use ncbi-blast; wu-blast is possible, for backwards
    # compatability, but it is now unsupported and orders of magnitude
    # slower than ncbi-blast. Setting BLASTMAT is irrelevant for ncbi-blast,
    # but the Ensembl code requires it to be set as an environment variable,
    # pointing to a directory that exists.
    blast_type       => 'ncbi',
    blast_dir        => '/nfs/panda/ensemblgenomes/external/ncbi-blast-2+',
    makeblastdb_exe  => catdir($self->o('blast_dir'), 'bin/makeblastdb'),
    blastp_exe       => catdir($self->o('blast_dir'), 'bin/blastp'),
    blastx_exe       => catdir($self->o('blast_dir'), 'bin/blastx'),
    blast_matrix     => undef,
    blast_threads    => 4,
    blast_parameters => '-evalue 0.0001 -seg no -soft_masking true -word_size 2 -num_alignments 500000 -num_descriptions 500000 -lcase_masking -num_threads '.$self->o('blast_threads'),
#   blast_parameters => '-evalue 0.0001 -num_alignments 400000 -num_descriptions 400000 -lcase_masking -seg yes -num_threads '.$self->o('blast_threads'),
    
    # For wu-blast, set the following  parameters instead of the above.
    # blast_type       => 'wu',
    # blast_dir        => '/nfs/panda/ensemblgenomes/external/wublast',
    # makeblastdb_exe  => catdir($self->o('blast_dir'), 'wu-formatdb'),
    # blastp_exe       => catdir($self->o('blast_dir'), 'blastp'),
    # blastx_exe       => catdir($self->o('blast_dir'), 'blastx'),
    # blast_matrix     => catdir($self->o('blast_dir'), 'matrix'),
    # blast_parameters => '-W 3 -B 100000 -V 100000 -hspmax=0 -lcmask -wordmask=seg',
    
    # For parsing the output.
    output_regex     => '^\s*([\w\-\.]+)',
    pvalue_threshold => 0.0001,
    filter_prune     => 1,
    filter_min_score => 100,
    
    # By default, do both blastp and blastx.
    # Specify blast_db if you don't want it in the directory alongside db_file.
    blast_db => undef,
    blastp   => 1,
    blastx   => 1,
    logic_name_prefix => 'worm',
    
    analyses =>
    [
      {
        'logic_name'    => 'gadflyp',
        'display_label' => $self->o('source_label'),
        'description'   => 'FlyBase peptide sequences aligned to the proteome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'gadfly',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/gadfly92.pep',
      },
      
      {                 
        'logic_name'    => 'gadflyx',
        'display_label' => $self->o('source_label'),
        'description'   => 'FlyBase peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'gadfly',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/gadfly92.pep',
      },
      {
        'logic_name'    => 'yeastp',
        'display_label' => $self->o('source_label'),
        'description'   => 'SGD peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'yeast',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/yeast63.pep',
      },
      
      {                 
        'logic_name'    => 'yeastx',
        'display_label' => $self->o('source_label'),
        'description'   => 'yeast peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'yeast',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/yeast63.pep',
      },
      {
        'logic_name'    => 'slimswissprotp',
        'display_label' => $self->o('source_label'),
        'description'   => 'UniProt/SW peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'slimswissprot',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/slimswissprot20169.pep',
      },
      
      {                 
        'logic_name'    => 'slimswissprotx',
        'display_label' => $self->o('source_label'),
        'description'   => 'UniProt/SW peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'slimswissprot',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/slimswissprot20169.pep',
      },
      {
        'logic_name'    => 'ipi_humanp',
        'display_label' => $self->o('source_label'),
        'description'   => 'Human IPI peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'ipi_human',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/ipi_human_07_20.pep',
      },
      
      {                 
        'logic_name'    => 'ipi_humanx',
        'display_label' => $self->o('source_label'),
        'description'   => 'Human IPI peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'ipi_human',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/ipi_human_07_20.pep',
      },
      {
        'logic_name'    => 'wormpepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.elegans peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'wormpep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/wormpep257.pep',
      },
      {                 
        'logic_name'    => 'wormpepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.elegans peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'wormpep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/wormpep257.pep',
      },
      {
        'logic_name'    => 'remapepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.remanei peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'remapep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/remapep256.pep',
      },
      {                 
        'logic_name'    => 'remapepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.remanei peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'remapep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/remapep256.pep',
      },
      {
        'logic_name'    => 'ppapepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'P.pacificus peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'ppapep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/ppapep250.pep',
      },
      {                 
        'logic_name'    => 'ppapepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'P.pacificus peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'ppapep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/ppapep250.pep',
      },
      {
        'logic_name'    => 'jappepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.japonica peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'jappep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/jappep250.pep',
      },
      {                 
        'logic_name'    => 'jappepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.japonica peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'jappep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/jappep250.pep',
      },
      {
        'logic_name'    => 'brepepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.brenneri peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brepep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/brepep250.pep',
      },
      {                 
        'logic_name'    => 'brepepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.brenneri peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brepep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/brepep250.pep',
      },
      {
        'logic_name'    => 'brigpepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.briggsae peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brigpep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/brigpep257.pep',
      },
      {                 
        'logic_name'    => 'brigpepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'C.briggsae peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brigpep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/brigpep257.pep',
      },
      {
        'logic_name'    => 'brugpepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'B.malayi peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brugpep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/brugpep257.pep',
      },
      {                 
        'logic_name'    => 'brugpepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'B.malayi peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brigpep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/brugpep257.pep',
      },
      {
        'logic_name'    => 'ovolpepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'O.volvulus peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'ovolpep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/ovolpep253.pep',
      },
      {                 
        'logic_name'    => 'ovolpepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'O.volvulus peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'brigpep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/ovolpep253.pep',
      },
      {
        'logic_name'    => 'srapepp',
        'display_label' => $self->o('source_label'),
        'description'   => 'S.ratti peptide sequences aligned to the proteome with <em>blastp</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'srapep',
        'db_version'    => '1',
        'program'       => 'blastp',
        'program_file'  => $self->o('blastp_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/srapep255.pep',
      },
      {                 
        'logic_name'    => 'srapepx',
        'display_label' => $self->o('source_label'),
        'description'   => 'S.ratti peptide sequences aligned to the genome with <em>blastx</em>.',
        'displayable'   => 1,
        'web_data'      => '{"type" => "protein"}',
        'db'            => 'srapep',
        'db_version'    => '1',
        'program'       => 'blastx',
        'program_file'  => $self->o('blastx_exe'),
        'parameters'    => $self->o('blast_parameters'),
        'module'        => 'Bio::EnsEMBL::Analysis::Runnable::BlastEG',
        'linked_tables' => ['protein_align_feature'],
        'db_type'       => 'core',
        'db_file'       => '/nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline/BlastDB/srapep255.pep',
      },
      
    ],

    # Remove existing *_align_features; if => 0 then existing analyses
    # and their features will remain, with the logic_name suffixed by '_bkp'.
    delete_existing => 1,

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 0,

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

sub pipeline_create_commands {
  my ($self) = @_;
  
  return [
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('pipeline_dir'),
    $self->db_cmd("CREATE TABLE split_genome (species varchar(100) NOT NULL, split_file varchar(255) NOT NULL)"),
    $self->db_cmd("CREATE TABLE split_proteome (species varchar(100) NOT NULL, split_file varchar(255) NOT NULL)"),
  ];
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
                            '1->A' => ['SpeciesFactory'],
                            'A->1' => ['TargetSpeciesFactory'],
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
                            '2->A' => ['TargetDatabase'], # <= Semaphore
                            'A->2' => ['GenBlastFactory'],
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
      -meadow_type     => 'LOCAL',
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
      -rc_name         => '8GB',
    },


    {
      -logic_name      => 'TargetDatabase',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::TargetDatabase',
      -max_retry_count => 1,
      -parameters      => {
                            blastp => $self->o('blastp'),
                            blastx => $self->o('blastx'),
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['DumpProteome'],
                            '3' => ['DumpGenome'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name        => 'DumpProteome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpProteome',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              proteome_dir => catdir($self->o('pipeline_dir'), '#species#', 'proteome'),
                              use_dbID     => 1,
                              allow_stop_codons => 1,
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SplitProteome'],
    },

    {
      -logic_name        => 'SplitProteome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              fasta_file              => '#proteome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => [ '?table_name=split_proteome' ],
                            }
    },

    {
      -logic_name        => 'DumpGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpGenome',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              genome_dir => catdir($self->o('pipeline_dir'), '#species#', 'genome'),
                            },
      -rc_name           => 'normal',
      -flow_into         => ['SplitGenome'],
    },

    {
      -logic_name        => 'SplitGenome',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::FastaSplit',
      -can_be_empty      => 1,
      -analysis_capacity => 5,
      -parameters        => {
                              fasta_file              => '#genome_file#',
                              max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                              max_seqs_per_file       => $self->o('max_seqs_per_file'),
                              max_files_per_directory => $self->o('max_files_per_directory'),
                              max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                            },
      -rc_name           => 'normal',
      -flow_into         => {
                              '2' => [ '?table_name=split_genome' ],
                            }
    },

    {
      -logic_name      => 'TargetSpeciesFactory',
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
                            '2' => ['AnalysisFactory'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'AnalysisFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::AnalysisFactory',
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            analyses     => $self->o('analyses'),
                            blastp       => $self->o('blastp'),
                            blastx       => $self->o('blastx'),
                            logic_name_prefix => $self->o('logic_name_prefix'),
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['AnalysisSetupCoreP'],
                            '3' => ['AnalysisSetupCoreX'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'AnalysisSetupCoreP',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::AnalysisSetup',
      -can_be_empty    => 1,
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            db_backup_required => 0,
                            db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'core_bkp.sql.gz'),
                            delete_existing    => $self->o('delete_existing'),
                            production_lookup  => $self->o('production_lookup'),
                            production_db      => $self->o('production_db'),
                          },
      -meadow_type     => 'LOCAL',
      -flow_into       => {
                            '2' => ['FetchProteomeFiles'],
                          },
    },
    {
      -logic_name      => 'AnalysisSetupCoreX',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::AnalysisSetup',
      -can_be_empty    => 1,
      -max_retry_count => 0,
      -batch_size      => 10,
      -parameters      => {
                            db_backup_required => 0,
                            db_backup_file     => catdir($self->o('pipeline_dir'), '#species#', 'core_bkp.sql.gz'),
                            delete_existing    => $self->o('delete_existing'),
                            production_lookup  => $self->o('production_lookup'),
                            production_db      => $self->o('production_db'),
                          },
      -meadow_type     => 'LOCAL',
      -flow_into       => {
                            '2' => ['FetchGenomeFiles'],
                          },
    },

    {
      -logic_name      => 'FetchProteomeFiles',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchSplitFiles',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            seq_type => 'proteome',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                           '2' => ['BlastPFactory'],
                          },
    },

    {
      -logic_name      => 'FetchGenomeFiles',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::FetchSplitFiles',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            seq_type => 'genome',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BlastXFactory'],
                          },
    },

    {
      -logic_name      => 'BlastPFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastFactory',
      -can_be_empty    => 1,
      -max_retry_count => 1,
      -parameters      => {
                            max_seq_length => $self->o('max_seq_length'),
                            queryfile      => '#split_file#',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BlastP'],
                          },
    },

    {
      -logic_name      => 'BlastXFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::BlastFactory',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            max_seq_length => $self->o('max_seq_length'),
                            queryfile      => '#split_file#',
                          },
      -rc_name         => 'normal',
      -flow_into       => {
                            '2' => ['BlastX'],
                          },
    },

    {
      -logic_name      => 'BlastP',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'pep',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '8GB_threads',
      -flow_into       => {
                            '-1' => ['BlastP_HighMem'],
                          },
    },

    {
      -logic_name      => 'BlastP_HighMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'pep',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '16GB_threads',
      -flow_into       => {
                            '-1' => ['BlastP_HigherMem'],
                          },
    },

    {
      -logic_name      => 'BlastP_HigherMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'pep',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                          },
      -rc_name         => '32GB_threads',
    },

    {
      -logic_name      => 'BlastX',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'dna',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '8GB_threads',
      -flow_into       => {
                            '-1' => ['BlastX_HighMem'],
                          },
    },

    {
      -logic_name      => 'BlastX_HighMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'dna',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                            escape_branch    => -1,
                          },
      -rc_name         => '16GB_threads',
      -flow_into       => {
                            '-1' => ['BlastX_HigherMem'],
                          },
    },

    {
      -logic_name      => 'BlastX_HigherMem',
      -module          => 'Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::Blast',
      -can_be_empty    => 1,
      -hive_capacity   => $self->o('max_hive_capacity'),
      -max_retry_count => 1,
      -parameters      => {
                            db_type          => 'core',
                            blast_type       => $self->o('blast_type'),
                            blast_matrix     => $self->o('blast_matrix'),
                            output_regex     => $self->o('output_regex'),
                            query_type       => 'dna',
                            database_type    => 'pep',
                            pvalue_threshold => $self->o('pvalue_threshold'),
                            filter_prune     => $self->o('filter_prune'),
                            filter_min_score => $self->o('filter_min_score'),
                          },
      -rc_name         => '32GB_threads',
    },

  ];
}

sub resource_classes {
  my ($self) = @_;
  
  my $blast_threads = $self->o('blast_threads');

  return {
    %{$self->SUPER::resource_classes},
    '8GB' => {'LSF' => '-q production-rh7 -M 8000 -R "rusage[mem=8000,tmp=8000]"'},
    '8GB_threads' => {'LSF' => '-q production-rh7 -n ' . ($blast_threads + 1) . ' -M 8000 -R "rusage[mem=8000,tmp=8000]"'},
    '16GB_threads' => {'LSF' => '-q production-rh7 -n ' . ($blast_threads + 1) . ' -M 16000 -R "rusage[mem=16000,tmp=16000]"'},
    '32GB_threads' => {'LSF' => '-q production-rh7 -n ' . ($blast_threads + 1) . ' -M 32000 -R "rusage[mem=32000,tmp=32000]"'},
  }
}

1;
