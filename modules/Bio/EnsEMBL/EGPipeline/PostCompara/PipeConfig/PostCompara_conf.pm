package Bio::EnsEMBL::EGPipeline::PostCompara::PipeConfig::PostCompara_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;

    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

        release  		=> software_version(),
	    registry        => [],
        compara         => 'plants',
        pipeline_name   => 'PostCompara_'.$self->o('compara').'_'.$self->o('release'),
        email           => $self->o('ENV', 'USER').'@ebi.ac.uk', 
        output_dir      => '/nfs/nobackup2/ensemblgenomes/'.$self->o('ENV', 'USER').'/workspace/'.$self->o('pipeline_name'),     

	## Flags controlling pipeline to run
	    # on '1' by default, set to '0' if want to skip this analysis
    	flag_GO         => '0',     
    	flag_GeneNames  => '1',     
    	
	## GeneName/Description Projection 
		# source species 
		gn_from_species  => 'arabidopsis_thaliana',

		# target species
	    gn_species       => ['vitis_vinifera'],
	    gn_antispecies   => [],
        gn_division 	 => [], # EnsemblMetazoa, EnsemblProtists, EnsemblFungi, EnsemblPlants
	    gn_run_all       => 0,

        taxon_filter     => 'eudicotyledons', # i.e Liliopsida,eudicotyledons
		geneName_source  => ['UniProtKB/Swiss-Prot', 'TAIR_SYMBOL'],
		#geneDesc_source  => ['UniProtKB/Swiss-Prot', 'TAIR_LOCUS', 'UniProtKB/TrEMBL'] ,

		gn_method_link_type       => 'ENSEMBL_ORTHOLOGUES',
		# only certain types of homology are considered
		gn_homology_types_allowed => ['ortholog_one2one'],
		
        # Percentage identify filter for the homology
        'gn_percent_id_filter'    => '10',

	## GO Projection  
		# source species 
		go_from_species  => 'arabidopsis_thaliana',,

		# target species
	    go_species       => ['vitis_vinifera'],
	    go_antispecies   => [],
        go_division 	 => [], # EnsemblMetazoa, EnsemblProtists, EnsemblFungi, EnsemblPlants
	    go_run_all       => 0,

		go_method_link_type       => 'ENSEMBL_ORTHOLOGUES',
		# only certain types of homology are considered
		go_homology_types_allowed => ['ortholog_one2one'],
		
        # Percentage identify filter for the homology
        'go_percent_id_filter'    => '10',

		# ensembl object type to attach GO projection, default 'Translation', options 'Transcript'
		ensemblObj_type           => 'Translation', 

        ## GOA webservice parameters
        goa_webservice   => 'http://www.ebi.ac.uk/QuickGO/',
		goa_params       => 'GValidate?service=taxon&action=getBlacklist&taxon=',

		## only certain types of homology are considered
		homology_types_allowed => ['ortholog_one2one','apparent_ortholog_one2one'],
		
		## only these evidence codes will be considered for GO term projection
		## See https://www.ebi.ac.uk/panda/jira/browse/EG-974
		evidence_codes         => ['IEA','IDA','IC','IGI','IMP','IPI','ISS','NAS','ND','RCA','TAS'],
		#  IC Inferred by curator
		#  IDA Inferred from direct assay
		#  IEA Inferred from electronic annotation
		#  IGI Inferred from genetic interaction
		#  IMP Inferred from mutant phenotype
		#  IPI Inferred from physical interaction
		#  ISS Inferred from sequence or structural similarity
		#  NAS Non-traceable author statement
		#  ND No biological data available
		#  RCA Reviewed computational analysis
		#  TAS Traceable author statement

	## For all pipelines
        ## flags
		'flag_store_projections' => '0', #  Off by default. Control the storing of projections into database. 
		'flag_backup'			 => '1', #  On by default. Dumping of table, backup to_species db. 
		
        'pipeline_db' => {  
     	   -host   => $self->o('host'),
           -port   => $self->o('port'),
           -user   => $self->o('user'),
           -pass   => $self->o('pass'),
           -dbname => $self->o('dbname'),
           -driver => 'mysql',
      },
		
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
    ];
}

# Ensures species output parameter gets propagated implicitly
sub hive_meta_table {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  return 
      ' -reg_conf ' . $self->o('registry')
  ;
}

sub pipeline_analyses {
    my ($self) = @_;
 
 	# Control which pipelines to run
  	my $pipeline_flow;
  
  	if ($self->o('flag_GO') && $self->o('flag_GeneNames')) {
    	$pipeline_flow  = ['GeneNamesProjectionFactory', 'GOProjectionFactory'];
  	} elsif ($self->o('flag_GO')) {
    	$pipeline_flow  = ['GOProjectionFactory'];
  	} elsif ($self->o('flag_GeneNames')) {
  	    $pipeline_flow  = ['GeneNamesProjectionFactory'];
  	}
 
    return [
      {  -logic_name    => 'backbone_fire_PostCompara',
         -module        => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
         -input_ids     => [ {} ], # Needed to create jobs
         -hive_capacity => -1,
		 -flow_into     => $pipeline_flow,
#         -flow_into 	=> {
#			'1' => ['GeneNamesProjectionFactory'],
#			'2' => ['GOProjectionFactory'],
#         },
      },

    {  -logic_name      => 'GeneNamesProjectionFactory',
        -module         => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
        -parameters     => {
                              species     => $self->o('gn_species'),
                              antispecies => $self->o('gn_antispecies'),
                              division    => $self->o('gn_division'),
                              run_all     => $self->o('gn_run_all'),
                            },
#       -input_ids       => [ {} ],
       -max_retry_count => 1,
       -rc_name         => 'default',
       -flow_into       => {
				             '2->A' => [ 'GeneNamesProjection' ],
				             'A->1' => [ 'NotifyUser' ],
                           },
    },
       
    {  -logic_name      => 'GOProjectionFactory',
        -module         => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
        -parameters     => {
                              species     => $self->o('go_species'),
                              antispecies => $self->o('go_antispecies'),
                              division    => $self->o('go_division'),
                              run_all     => $self->o('go_run_all'),
                            },
#       -input_ids       => [ {} ],
       -max_retry_count => 1,
       -rc_name         => 'default',
       -flow_into       => {
				             '2->A' => [ 'GOProjection' ],
				             'A->1' => [ 'NotifyUser' ],
                           },
    },


    {  -logic_name => 'GeneNamesProjection',
       -module     => 'Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::GeneNamesProjection',
       -parameters => {
			'geneName_source'		  => $self->o('geneName_source'),  
#			'geneDesc_source'		  => $self->o('geneDesc_source'),  
		    'taxon_filter'			  => $self->o('taxon_filter'),
		    'from_species'            => $self->o('gn_from_species'),
		    'compara'                 => $self->o('compara'),
   		    'release'                 => $self->o('release'),
   		    'method_link_type'        => $self->o('gn_method_link_type'),
   		    'homology_types_allowed ' => $self->o('gn_homology_types_allowed'),
            'percent_id_filter'       => $self->o('gn_percent_id_filter'),
            'output_dir'              => $self->o('output_dir'),
   	   },
       -rc_name       => 'default',
    },

    {  -logic_name => 'GOProjection',
         -module     => 'Bio::EnsEMBL::EGPipeline::PostCompara::RunnableDB::GOProjection',
         -parameters => {
		    'from_species'            => $self->o('go_from_species'),
		    #'division'                => $self->o('division'),
   		    'release'                 => $self->o('release'),
   		    'ensemblObj_type'		  => $self->o('ensemblObj_type'),
   		    'goa_webservice'          => $self->o('goa_webservice'),
   		    'goa_params'              => $self->o('goa_params'),
   		    'method_link_type'        => $self->o('go_method_link_type'),
   		    'homology_types_allowed ' => $self->o('go_homology_types_allowed'),
   		    'evidence_codes'		  => $self->o('evidence_codes'),
   		    'percent_id_filter'       => $self->o('go_percent_id_filter'),
            'output_dir'              => $self->o('output_dir'),
     	 },
#         -hive_capacity => $self->o('blastp_capacity'),
#         -batch_size    =>  50, 
          -rc_name       => 'default',
	  },


    {  -logic_name => 'NotifyUser',
       -module     => 'Bio::EnsEMBL::EGPipeline::ProjectGeneNames::RunnableDB::NotifyUser',
       -parameters => {
          	'email'      => $self->o('email'),
          	'subject'    => $self->o('pipeline_name').' has finished',
          	'output_dir' => $self->o('output_dir'),
       },
    }
  ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
            'flag_store_projections' => $self->o('flag_store_projections'),
       		'flag_backup'            => $self->o('flag_backup'),
    };
}

# override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
    my $self = shift;    
    return "-reg_conf ".$self->o("registry");
}

sub resource_classes {
    my $self = shift;
    return {
      'default'  	 => { 'LSF' => '-q production-rh6 -n 4 -M 4000 -R "rusage[mem=4000]"'},
      'mem'     	 => { 'LSF' => '-q production-rh6 -n 4 -M 12000 -R "rusage[mem=12000]"'},
      '2Gb_job'      => {'LSF' => '-q production-rh6 -C0 -M2000  -R"select[mem>2000]  rusage[mem=2000]"' },
      '24Gb_job'     => {'LSF' => '-q production-rh6 -C0 -M24000 -R"select[mem>24000] rusage[mem=24000]"' },
      '250Mb_job'    => {'LSF' => '-q production-rh6 -C0 -M250   -R"select[mem>250]   rusage[mem=250]"' },
      '500Mb_job'    => {'LSF' => '-q production-rh6 -C0 -M500   -R"select[mem>500]   rusage[mem=500]"' },
	  '1Gb_job'      => {'LSF' => '-q production-rh6 -C0 -M1000  -R"select[mem>1000]  rusage[mem=1000]"' },
	  '2Gb_job'      => {'LSF' => '-q production-rh6 -C0 -M2000  -R"select[mem>2000]  rusage[mem=2000]"' },
	  '8Gb_job'      => {'LSF' => '-q production-rh6 -C0 -M8000  -R"select[mem>8000]  rusage[mem=8000]"' },
	  '24Gb_job'     => {'LSF' => '-q production-rh6 -C0 -M24000 -R"select[mem>24000] rusage[mem=24000]"' },
	  'msa'          => {'LSF' => '-q production-rh6 -W 24:00' },
	  'msa_himem'    => {'LSF' => '-q production-rh6 -M 32768 -R "rusage[mem=32768]" -W 24:00' },
	  'urgent_hcluster'      => {'LSF' => '-q production-rh6 -C0 -M8000  -R"select[mem>8000]  rusage[mem=8000]"' },
    }
}


1;
