=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::EGPipeline::FileDump::WriteDrupalFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);

sub param_defaults {  
  return {
    'drupal_desc' => {
      'fasta_toplevel'    => '<STRAIN> strain genomic <SEQTYPE> sequences, <ASSEMBLY> assembly, softmasked using RepeatMasker, Dust, and TRF.',
      'fasta_seqlevel'    => '<STRAIN> strain genomic <SEQTYPE> sequences, <ASSEMBLY> assembly.',
      'agp_assembly'      => 'AGP (v2.0) file relating <MAPPING> for the <SPECIES> <STRAIN> strain, <ASSEMBLY> assembly.',
      'fasta_transcripts' => '<STRAIN> strain transcript sequences, <GENESET> geneset.',
      'fasta_peptides'    => '<STRAIN> strain peptide sequences, <GENESET> geneset.',
      'gtf_genes'         => '<STRAIN> strain <GENESET> geneset in GTF (v2.2) format.',
      'gff3_genes'        => '<STRAIN> strain <GENESET> geneset in GFF3 format.',
      'gff3_repeats'      => '<STRAIN> strain <ASSEMBLY> repeat features (RepeatMasker, Dust, TRF) in GFF3 format.',
    },
    
    'drupal_desc_exception' => {
      'fasta_toplevel' => {
        'Musca domestica' => '<STRAIN> strain genomic <TOPLEVEL> sequences, <ASSEMBLY> assembly, softmasked using WindowMasker, Dust, and TRF.',
      },
      'gff3_repeats' => {
        'Musca domestica' => '<STRAIN> strain <GENESET> repeat features (WindowMasker, Dust, TRF) in GFF3 format.',
      },
    },
    
    'gene_dumps' => [
      'fasta_transcripts',
      'fasta_peptides',
      'gtf_genes',
      'gff3_genes',
    ],
    
    'drupal_species' => {
      'Anopheles culicifacies' => 'Anopheles culicifacies A',
    },
    
    'staging_dir'  => 'sites/default/files/ftp/staging',
  };
}


sub run {
  my ($self) = @_;
  
  my $results_dir = $self->param_required('results_dir');
  my $drupal_file = $self->param_required('drupal_file');
  
  opendir(my $dh, $results_dir) || die "Failed to open '$results_dir': $!";
  my @file_names = grep { !/.md5/ && !/^\./ && -f "$results_dir/$_" } readdir($dh);
  closedir $dh;
  
  my @fields = (
    'GUID', 'File', 'Organism', 'File Type', 'File Format', 'Status',
    'Description', 'Latest Change', 'Tags', 'Release Date Start',
    'Release Date End', 'Version', 'Display Version', 'Previous Version',
    'Xgrid_enabled', 'Fasta Header Regex', 'Download Count',
    'Title', 'URL', 'Ensembl organism name'
  );
  
  my %data = map { $_ => [] } @fields;
  my $guid = 1;
  
  foreach my $file_name (sort @file_names) {
    $self->process_file($file_name, $guid++, \%data);
  }
  
  open(my $fh, '>', $drupal_file) || die "Failed to open '$drupal_file': $!";
  
  print $fh join(",", @fields)."\n";
  
  for (my $i=0; $i<($guid-1); $i++) {
    my @row;
    foreach my $column (@fields) {
      push @row, '"' . $data{$column}[$i] . '"';
    }
    print $fh join(",", @row)."\n";
  }
}

sub process_file {
  my ($self, $file_name, $guid, $data) = @_;
  
  my $staging_dir  = $self->param_required('staging_dir');
  my $release_date = $self->param_required('release_date');
  
  my ($species, $strain, $data_type, $assembly, $geneset, $dump_type) =
    $self->parse_filename($file_name);
  
  my $organism = $self->organism($species);
  my $file_type = $self->file_type($data_type);
  my $file_format = $self->file_format($dump_type);
  my $description = $self->description($dump_type, $data_type, $species, $strain, $assembly, $geneset);
  my $display_version = $self->display_version($dump_type, $assembly, $geneset);
  my $xgrid = $self->xgrid($dump_type);
  
  push $$data{'GUID'}, $guid;
  push $$data{'File'}, catdir($staging_dir, $file_name);
  push $$data{'Organism'}, $organism;
  push $$data{'File Type'},$file_type ;
  push $$data{'File Format'},$file_format ;
  push $$data{'Status'}, 'Current';
  push $$data{'Description'}, $description;
  push $$data{'Latest Change'}, '';
  push $$data{'Tags'}, '';
  push $$data{'Release Date Start'}, $release_date;
  push $$data{'Release Date End'}, $release_date;
  push $$data{'Version'}, '1';
  push $$data{'Display Version'}, $display_version;
  push $$data{'Previous Version'}, '';
  push $$data{'Xgrid_enabled'}, $xgrid;
  push $$data{'Fasta Header Regex'}, '(.*?)\s';
  push $$data{'Download Count'}, '0';
  push $$data{'Title'}, '';
  push $$data{'URL'}, '';
  push $$data{'Ensembl organism name'}, $species;
}

sub parse_filename {
  my ($self, $file_name) = @_;
  
  my ($species, $strain, $data_type, $assembly, $geneset_version, $file_type) =
    $file_name =~ /^(\w+\-\w+)\-([\w\-\.]+)_([A-Z0-9]+)_(\w+)(\.\d+)?\.(\w+)/;
  
  my $geneset = '';
  if ($geneset_version) {
    $geneset = "$assembly$geneset_version";
  }
  $species =~ s/\-/ /;
  
  my ($dump_type);
  if ($file_type eq 'agp') {
    $dump_type = 'agp_assembly';
    
  } elsif ($file_type eq 'gtf') {
    $dump_type = 'gtf_genes';
    
  } elsif ($file_type eq 'gff3') {
    if ($data_type eq 'BASEFEATURES') {
      $dump_type = 'gff3_genes';
    } elsif ($data_type eq 'REPEATFEATURES') {
      $dump_type = 'gff3_repeats';
    }
    
  } elsif ($file_type eq 'fa') {
    if ($data_type eq 'TRANSCRIPTS') {
      $dump_type = 'fasta_transcripts';
    } elsif ($data_type eq 'PEPTIDES') {
      $dump_type = 'fasta_peptides';
    } elsif ($data_type eq 'CHROMOSOMES') {
      $dump_type = 'fasta_toplevel';
    } elsif ($data_type eq 'SCAFFOLDS') {
      $dump_type = 'fasta_toplevel';
    } elsif ($data_type eq 'CONTIGS') {
      $dump_type = 'fasta_seqlevel';
    }
    
  }
  
  return ($species, $strain, $data_type, $assembly, $geneset, $dump_type);
}

sub organism {
  my ($self, $species) = @_;
  
  my $drupal_species = $self->param_required('drupal_species');
  my $organism = $species;
  if (exists $$drupal_species{$species}) {
    $organism = $$drupal_species{$species};
  }
  
  return $organism;
}

sub file_type {
  my ($self, $data_type) = @_;
  
  my $file_type;
  if ($data_type eq 'REPEATFEATURES') {
    $file_type = 'Repeat features';
  } elsif ($data_type =~ /(CONTIG|SCAFFOLD)2(SCAFFOLD|CHROMOSOME)/) {
    $file_type = ucfirst(lc($1)) . ' to ' . ucfirst(lc($2)) . ' mapping';
  } else {
    $file_type = ucfirst(lc($data_type));
  } 
  return $file_type;
}

sub file_format {
  my ($self, $dump_type) = @_;
  
  my ($file_format) = $dump_type =~ /^([a-z0-9]+)/;
  if ($file_format eq 'fasta') {
    $file_format = ucfirst($file_format);
  } else {
    $file_format = uc($file_format);
  }
  return $file_format;
}

sub description {
  my ($self, $dump_type, $data_type, $species, $strain, $assembly, $geneset) = @_;
  
  my $drupal_desc           = $self->param_required('drupal_desc');
  my $drupal_desc_exception = $self->param_required('drupal_desc_exception');
  
  my $description;
  if (exists $$drupal_desc_exception{$dump_type}{$species}) {
    $description = $$drupal_desc_exception{$dump_type}{$species};
  } else {
    $description = $$drupal_desc{$dump_type};
  }
  
  my $seqtype = lc($data_type);
  my ($from, $to) = $data_type =~ /(\w+)2(\w+)/;
  my $mapping = ($from && $to) ? lc($from).'s to '.lc($to).'s' : '';
  
  $description =~ s/<ASSEMBLY>/$assembly/;
  $description =~ s/<GENESET>/$geneset/;
  $description =~ s/<MAPPING>/$mapping/;
  $description =~ s/<SEQTYPE>/$seqtype/;
  $description =~ s/<SPECIES>/$species/;
  $description =~ s/<STRAIN>/$strain/;
  
  return $description;
}

sub display_version {
  my ($self, $dump_type, $assembly, $geneset) = @_;
  
  my $gene_dumps = $self->param_required('gene_dumps');
  my %gene_dumps = map { $_ => 1 } @$gene_dumps;
  
  my $display_version;
  if (exists $gene_dumps{$dump_type}) {
    $display_version = $geneset;
  } else {
    $display_version = $assembly;
  }
  return $display_version;
}

sub xgrid {
  my ($self, $dump_type) = @_;
  
  my $xgrid = 'No';
  if ($dump_type =~ /fasta/) {
    $xgrid = 'Yes';
  }
  return $xgrid;
}

1;
