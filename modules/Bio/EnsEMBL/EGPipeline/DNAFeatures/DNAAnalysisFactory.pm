=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::DNAAnalysisFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'no_dust'              => 0,
    'no_repeatmasker'      => 0,
    'no_trf'               => 0,
    'repeatmasker_library' => undef,
  };
}

sub run {
  my ($self) = @_;
  my $species = $self->param_required('species');
  my $dna_analyses = $self->param_required('dna_analyses');
  my $repeatmasker_library = $self->param_required('repeatmasker_library');
  my $logic_names = $self->param_required('logic_name');
  my $always_use_repbase = $self->param_required('always_use_repbase');
  my $pipeline_dir = $self->param_required('pipeline_dir');
  my $db_backup_file = $self->param('db_backup_file');
  
  # RepeatMasker requires a species name with a space.
  (my $species_rm = ucfirst($species)) =~ s/_/ /;
  
  my $filtered_analyses = [];
  foreach my $analysis (@{$dna_analyses}) {
    my $logic_name = $$analysis{'logic_name'};
    
    if ($logic_name eq 'dust') {
      push @$filtered_analyses, $analysis unless $self->param('no_dust');
      
    } elsif ($logic_name eq 'repeatmask') {
      unless ($self->param('no_repeatmasker')) {
        if ($always_use_repbase || (! exists $$repeatmasker_library{$species} && ! exists $$repeatmasker_library{'all'})) {
          if ($self->check_repeatmasker($analysis, $species_rm, $pipeline_dir)) {
            $$analysis{'parameters'} .= " -species \"$species_rm\"";
          }
          push @$filtered_analyses, $analysis;
        }
      }
      
    } elsif ($logic_name eq 'repeatmask_customlib') {
      unless ($self->param('no_repeatmasker')) {
        if (exists $$repeatmasker_library{$species} || exists $$repeatmasker_library{'all'}) {
          if (exists $$logic_names{$species}) {
            $$analysis{'logic_name'} = $$logic_names{$species};
          }
          my $library = $$repeatmasker_library{$species} || $$repeatmasker_library{'all'};
          $$analysis{'db_file'} = $library;
          $$analysis{'parameters'} .= ' -lib "'.$library.'"';
          push @$filtered_analyses, $analysis;
        }
      }
      
    } elsif ($logic_name eq 'trf') {
      push @$filtered_analyses, $analysis unless $self->param('no_trf');
    } else {
      $self->warning("Unrecognised analysis '$logic_name' in 'dna_analyses' parameter.");
    }
  }
  
  foreach my $analysis (@{$filtered_analyses}) {
    $$analysis{'species'} = $species;
    $$analysis{'db_backup_file'} = $db_backup_file if defined $db_backup_file;
  }
  
  $self->param('filtered_analyses', $filtered_analyses);
  
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('filtered_analyses'), 2);
  
}

# If RepeatMasker hasn't got a species in its taxonomy, it will fail.
sub check_repeatmasker {
  my ($self, $analysis, $species_rm, $pipeline_dir) = @_;
  
  my $program_file = $$analysis{'program_file'};
  my $parameters = $$analysis{'parameters'};
  $parameters .= " -species \"$species_rm\"";
  
  my $test_fasta_file = "$pipeline_dir/rm_test.fa";
  open (FASTA, '>$test_fasta_file') or die "Failed to create test file '$test_fasta_file'";
  print FASTA ">test\nATTATT\n";
  close (FASTA); 
  my $rm_out = `$program_file $parameters $test_fasta_file 2>&1`;
  unlink $test_fasta_file;
  
  if ($rm_out =~ /Species.*is not known to RepeatMasker/) {
    return 0;
  } else {
    return 1;
  }
}

1;

