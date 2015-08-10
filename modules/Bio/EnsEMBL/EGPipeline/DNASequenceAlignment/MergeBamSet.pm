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

package Bio::EnsEMBL::EGPipeline::DNASequenceAlignment::MergeBamSet;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::EGPipeline::Common::Aligner;
use Bio::EnsEMBL::ENA::SRA::BaseSraAdaptor qw(get_adaptor);
use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'samtools_dir' => '/nfs/panda/ensemblgenomes/external/samtools',
    'vcf'          => 0,
    'use_csi'      => 0,
    'clean_up'     => 1,
  };
}

sub fetch_input {
	my ($self) = @_;
  
	my $merge = $self->param_required('merge');
  # Invert the hash and gather old keys into a new value list.
  my %merge_bam;
  foreach my $bam_file (keys %$merge) {
    my $merge_id = $$merge{$bam_file};
    push @{$merge_bam{$merge_id}}, $bam_file;
  }
  
  $self->param('merge_bam', \%merge_bam);
}
  
sub run {
	my ($self) = @_;
  
  my $work_dir     = $self->param_required('work_directory');
  my $samtools_dir = $self->param_required('samtools_dir');
	my $vcf          = $self->param_required("vcf");
	my $use_csi      = $self->param_required("use_csi");
  my $clean_up     = $self->param_required('clean_up');
	my $merge_bam    = $self->param_required('merge_bam');
  
  my $aligner = Bio::EnsEMBL::EGPipeline::Common::Aligner->new(
    -samtools_dir => $samtools_dir,
  );
  
  my @output_files;
  
  foreach my $merge_id (keys %$merge_bam) {
    my $merged_bam_file = catdir($work_dir, "/$merge_id.bam");
    my @bam_files = @{$$merge_bam{$merge_id}};
    my $size = scalar @bam_files;
    
    if ($size == 1) {
      my ($bam_file) = $bam_files[0];
      rename $bam_file, $merged_bam_file;
    } else {
      $aligner->merge_bam(\@bam_files, $merged_bam_file);
      if ($clean_up) {
        map { unlink $_ } @bam_files;
			}
		}
    
    my $sorted_bam = $aligner->sort_bam($merged_bam_file);
    rename $sorted_bam, $merged_bam_file;
    $aligner->index_bam($merged_bam_file, $use_csi);
    if ($vcf) {
      $aligner->generate_vcf($merged_bam_file);
    }
    
    my $ini_file = $self->write_inifile($merge_id, $merged_bam_file);
    push @output_files,
      {
        merged_bam_file  => $merged_bam_file,
        ini_file => $ini_file,
      };
	}
  
  $self->param('output_files', \@output_files);
}

sub write_output {
  my ($self) = @_;
  
  foreach my $output (@{$self->param('output_files')}) {
    $self->dataflow_output_id($output, 1);
  }
}

sub write_inifile {
  my ($self, $merge_id, $merged_bam_file) = @_;
  
  my $mode        = $self->param_required('mode');
  my $work_dir    = $self->param_required('work_directory');
  my $merge_level = lc($self->param_required('merge_level'));
  
  my ($name, $caption, $description);
  
  if ($mode eq 'file') {
    ($name, $caption, $description) = $self->file_ini($merge_id, $merge_level);
  } elsif ($mode eq 'study') {
    ($name, $caption, $description) = $self->study_ini($merge_id, $merge_level);
  } else {
    $self->throw("Unrecognised mode of operation, '$mode'");
  }
  
  my $header = "$merge_id = rnaseq_align\n";
  my $body = 
    "[$merge_id]\n".
    "source_name = $name\n".
    "caption     = $caption\n".
    "description = $description\n".
    "source_url  = $merged_bam_file\n".
    "source_type = bam\n".
    "display     = off\n";
  
  my $ini_file = catdir($work_dir, "/$merge_id.ini");
  open (my $fh, '>', $ini_file) or die "Failed to open file '$ini_file'";
	print $fh "[ENSEMBL_INTERNAL_BAM_SOURCES]\n$header\n$body\n";
  close($fh);
  
  return $ini_file;
}

sub file_ini {
  my ($self, $merge_id, $merge_level) = @_;
  
  my $name;
  if ($merge_level eq 'file') {
    $name = $merge_id;
  } else {
    my $species      = $self->param('species');
    my $seq_file     = $self->param('seq_file');
    my $species_file = $self->param('species_file');
    if (exists $$species_file{$species}) {
      push @$seq_file, $$species_file{$species};
    }
    $name = join('; ', @$seq_file);
  }
  
  return ($name, '', '');
}

sub study_ini {
  my ($self, $merge_id, $merge_level) = @_;
  
  my ($name, $caption, $description);
  my $ena_link = "<a href='http://www.ebi.ac.uk/ena/data/view/$merge_id' target='_blank'>$merge_id</a>";
  
  my $species       = $self->param('species');
  my $study_list    = $self->param('study');
  my $species_study = $self->param('species_study');
  if (exists $$species_study{$species}) {
    push @$study_list, $$species_study{$species};
  }
  
  $caption = $merge_id;
  
  my $study_adaptor = get_adaptor('Study');
  my @study_details;
  foreach my $study_acc (@$study_list) {
    foreach my $study (@{$study_adaptor->get_by_accession($study_acc)}) {
      (my $study_link = $ena_link) =~ s/$merge_id/$study_acc/g;
      push @study_details, $study->title()." (Study  $study_link)";
    }
  }
  $description = join('; ', @study_details);
  
  if ($merge_level eq 'study') {
    my ($study) = @{$study_adaptor->get_by_accession($merge_id)};
    $name = $study->title();
    $description = "$name (Study $ena_link)";
    
  } elsif ($merge_level eq 'taxon') {
    my $taxon = $study_adaptor->fetch_by_taxon_id($merge_id);
    $name = $taxon->name();
    
  } elsif ($merge_level eq 'sample') {
    my $sample_adaptor = get_adaptor('Sample');
    my ($sample) = @{$sample_adaptor->get_by_accession($merge_id)};
    $name = $sample->title() || $sample->description();
    
  } elsif ($merge_level eq 'experiment') {
    my $experiment_adaptor = get_adaptor('Experiment');
    my ($experiment) = @{$experiment_adaptor->get_by_accession($merge_id)};
    $name = $experiment->title();
    
  } else {
    $name = join('; ', @$study_list);
    
  }
    
  return ($name, $caption, $description);
}

1;
