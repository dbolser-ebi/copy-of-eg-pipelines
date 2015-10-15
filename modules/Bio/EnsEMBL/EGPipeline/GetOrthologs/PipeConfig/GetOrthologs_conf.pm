package Bio::EnsEMBL::EGPipeline::GetOrthologs::PipeConfig::GetOrthologs_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');
use Bio::EnsEMBL::Hive::Version 2.3;
use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;

    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },      

		'registry'  	   => '',
        'pipeline_name'    => $self->o('hive_dbname'),       
        'output_dir'       => '/nfs/ftp/pub/databases/ensembl/projections/'.$self->o('ENV', 'USER').'/workspace/'.$self->o('pipeline_name'),     
		'method_link_type' => 'ENSEMBL_ORTHOLOGUES',

        # hive_capacity values for analysis
	    'getOrthologs_capacity'  => '50',

	 	'species_config' => 
		{ 
	 	  '1'=>{
	 	  		# compara database to get orthologs from
	 	  		#  'plants', 'protists', 'fungi', 'metazoa', 'multi' 
	 	  		'compara' => '',
	 	  		# source species  
	 	  		'source'  => '', # 'schizosaccharomyces_pombe'  	  		
				# target species
				'target'  => undef, # ['aspergillus_nidulans','puccinia_graminis']
				}, 

#	 	  '2'=>{
#	 	  		'compara'  => '',
#	 	  		'source'   => '',  	
#	 	  		'target'   => undef,  	
#	 	  	    },  		
    	},

       'pipeline_db' => {  
		     -host   => $self->o('hive_host'),
        	 -port   => $self->o('hive_port'),
        	 -user   => $self->o('hive_user'),
        	 -pass   => $self->o('hive_password'),
	         -dbname => $self->o('hive_dbname'),
        	 -driver => 'mysql',
      	},
		
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('output_dir'),
    ];
}

# Ensures output parameters gets propagated implicitly
sub hive_meta_table {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

# override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  return 
      ' -reg_conf ' . $self->o('registry'),
  ;
}

sub pipeline_analyses {
    my ($self) = @_;
 
    return [
    {  -logic_name    => 'backbone_fire_GetOrthologs',
       -module        => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -input_ids     => [ {} ] , 
       -flow_into 	  => { '1' => ['SourceFactory'], }
    },   
 
    {  -logic_name    => 'SourceFactory',
       -module        => 'Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::SourceFactory',
       -parameters    => { 'species_config'  => $self->o('species_config'), }, 
       -flow_into     => { '2' => ['MLSSJobFactory'], },          
       -rc_name       => 'default',
    },    
 
    {  -logic_name    => 'MLSSJobFactory',
       -module        => 'Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::MLSSJobFactory',
       -parameters    => { 'method_link_type' => $self->o('method_link_type'), },
       -flow_into     => { '2' => ['GetOrthologs'], },
       -rc_name       => 'default',
    },
  
    {  -logic_name    => 'GetOrthologs',
       -module        => 'Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::GetOrthologs',
       -parameters    => {	'output_dir'             => $self->o('output_dir'),
							'method_link_type'       => $self->o('method_link_type'),
    	 				 },
       -batch_size    =>  1,
       -rc_name       => 'default',
	   -hive_capacity => $self->o('getOrthologs_capacity'), 
	   -flow_into     => { '-1' => 'GetOrthologs_16GB', }, 
	 },
	 
    {  -logic_name    => 'GetOrthologs_16GB',
       -module        => 'Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::GetOrthologs',
       -parameters    => {	'output_dir'             => $self->o('output_dir'),
							'method_link_type'       => $self->o('method_link_type'),
    	 				 },
       -batch_size    =>  1,
       -rc_name       => '16Gb_mem',
	   -hive_capacity => $self->o('getOrthologs_capacity'), 
	   -flow_into     => { '-1' => 'GetOrthologs_32GB', }, 
	 },

    {  -logic_name    => 'GetOrthologs_32GB',
       -module        => 'Bio::EnsEMBL::EGPipeline::GetOrthologs::RunnableDB::GetOrthologs',
       -parameters    => {	'output_dir'             => $self->o('output_dir'),
							'method_link_type'       => $self->o('method_link_type'),
    	 				 },
       -batch_size    =>  1,
       -rc_name       => '32Gb_mem',
	   -hive_capacity => $self->o('getOrthologs_capacity'), 
	 },
	 	 
  ];
}

1;
