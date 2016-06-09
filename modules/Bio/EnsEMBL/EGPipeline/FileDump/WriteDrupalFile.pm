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
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'compara_files' => 0,
  };
}

sub run {
  my ($self) = @_;
  my $new_files     = $self->param_required('new_files');
  my $missing_files = $self->param_required('missing_files');
  my $changed_files = $self->param_required('changed_files');
  
  my @fields = (
    'Title', 'File', 'Organism', 'File Type', 'File Format', 'Status',
    'Description', 'Latest Change', 'Tags', 'Release Date Start',
    'Release Date End', 'Version', 'Display Version', 'Previous Version',
    'Xgrid_enabled', 'Fasta Header Regex', 'Download Count',
    'URL', 'Ensembl organism name', 'md5'
  );
  
  my %new = map { $_ => [] } @fields;
  foreach my $file (sort @$new_files) {
    $self->process_new_file($file, \%new);
  }
  $self->drupal_bulk_update(\@fields, \%new);
  
  my %missing;
  foreach my $file (sort @$missing_files) {
    $self->process_missing_file($file, \%missing);
  }
  $self->drupal_manual_update(\%missing);
  
  $self->drupal_shell_cmds($new_files, $changed_files);
}

sub process_new_file {
  my ($self, $file, $data) = @_;
  
  my $staging_dir   = $self->param_required('staging_dir');
  my $release_date  = $self->param_required('release_date');
  my $compara_files = $self->param_required('compara_files');
  
  my ($species, $strain, $data_type, $assembly, $geneset, $dump_type);
  my ($organism, $file_type, $file_format, $description, $display_version, $xgrid);
  
  if ($compara_files) {
    $species = '';
    $organism = '';
    ($file_type, $file_format) = $self->parse_compara_filename($file);
    $description = $self->compara_description($file);
    $display_version = $release_date;
    $xgrid = 0;
    
  } else {
    ($species, $strain, $data_type, $assembly, $geneset, $dump_type) = $self->parse_filename($file);
    $organism = $self->organism($species);
    $file_type = $self->file_type($data_type);
    $file_format = $self->file_format($dump_type);
    $description = $self->description($dump_type, $data_type, $species, $strain, $assembly, $geneset);
    $display_version = $self->display_version($dump_type, $assembly, $geneset);
    $xgrid = $self->xgrid($dump_type);
  }
  
  push $$data{'Title'}, $file;
  push $$data{'File'}, catdir($staging_dir, $file);
  push $$data{'Organism'}, $organism;
  push $$data{'File Type'}, $file_type;
  push $$data{'File Format'}, $file_format;
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
  push $$data{'URL'}, '';
  push $$data{'Ensembl organism name'}, $species;
  push $$data{'md5'}, catdir($staging_dir, "$file.md5");
}

sub drupal_bulk_update {
  my ($self, $fields, $data) = @_;
  my $drupal_file = $self->param_required('drupal_file');
  
  open(my $fh, '>', $drupal_file) || die "Failed to open '$drupal_file': $!";
  
  print $fh join(",", @$fields)."\n";
  
  my $rows = scalar(@{$$data{'Title'}});
  
  for (my $i=0; $i<=$rows; $i++) {
    my @row;
    foreach my $column (@$fields) {
      push @row, '"' . $$data{$column}[$i] . '"';
    }
    print $fh join(",", @row)."\n";
  }
  
  close $fh;
}

sub process_missing_file {
  my ($self, $file, $data) = @_;
  
  push @{$$data{$file}}, ['Status', 'Archived'];
}

sub drupal_manual_update {
  my ($self, $data) = @_;
  my $manual_file = $self->param_required('manual_file');
  
  open(my $fh, '>', $manual_file) || die "Failed to open '$manual_file': $!";
  
  print $fh "Manual updates to the following pages are required:\n";
  
  foreach my $title (keys %$data) {
    print $fh "$title\n";
    foreach my $field_value (@{$$data{$title}}) {
      print $fh "\t".join(": ", @$field_value)."\n";
    }
  }
  
  close $fh;
}

sub drupal_shell_cmds {
  my ($self, $new_files, $changed_files) = @_;
  my $results_dir      = $self->param_required('results_dir');
  my $sh_ebi_file      = $self->param_required('sh_ebi_file');
  my $sh_nd_file       = $self->param_required('sh_nd_file');
  my $nd_login         = $self->param_required('nd_login');
  my $nd_downloads_dir = $self->param_required('nd_downloads_dir');
  my $nd_staging_dir   = $self->param_required('nd_staging_dir');
  my $release_date     = $self->param_required('release_date');
  
  open(my $ebi, '>', $sh_ebi_file) || die "Failed to open '$sh_ebi_file': $!";
  open(my $nd, '>', $sh_nd_file) || die "Failed to open '$sh_nd_file': $!";
  
  if (@$changed_files) {
    print $nd "# Archive changed files:\n";
    print $nd "cd $nd_downloads_dir\n";
    print $nd "mkdir archive/$release_date\n";
    foreach my $file (@$changed_files) {
      
      print $nd "mv $file archive/$release_date\n";
      print $nd "rm -f $file.md5\n";
    }
    
    my $changed_dir = "$results_dir/changed";
    mkdir $changed_dir;
    
    print $ebi "# Copy changed files to the downloads directory at ND:\n";
    print $ebi "cd $changed_dir\n";
    foreach my $file (@$changed_files) {
      print $ebi "ln -s ../$file";
      print $ebi "ln -s ../$file.md5";
    }
    print $ebi "scp $changed_dir/* $nd_login:$nd_downloads_dir\n";
  }
  
  if (@$new_files) {
    print $nd "# Clear the staging directory:\n";
    print $nd "rm -rf $nd_staging_dir/*\n";
    
    my $new_dir = "$results_dir/new";
    mkdir $new_dir;
    
    print $ebi "# Copy new files to the staging directory at ND:\n";
    print $ebi "cd $new_dir\n";
    foreach my $file (@$new_files) {
      print $ebi "ln -s ../$file";
      print $ebi "ln -s ../$file.md5";
    }
    print $ebi "scp $new_dir/* $nd_login:$nd_staging_dir\n";
  }
  
  close $ebi;
  close $nd;
}

sub parse_filename {
  my ($self, $file) = @_;
  
  my ($species, $strain, $data_type, $assembly, $geneset_version, $file_type) =
    $file =~ /^(\w+\-\w+)\-([\w\-\.]+)_([A-Z0-9]+)_(\w+)(\.\d+)?\.(\w+)/;
  
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

sub parse_compara_filename {
  my ($self, $file) = @_;
  
  my ($prefix) = $file =~ /^([^_]+)/;
  
  if ($prefix eq 'GENE-TREES-NEWICK') {
    return ('Gene trees', 'Newick (tar)');
  } elsif ($prefix eq 'GENE-ALIGN-TRANSCRIPTS') {
    return ('Gene alignments', 'Fasta (tar)');
  } elsif ($prefix eq 'GENE-ALIGN-PEPTIDES') {
    return ('Gene alignments', 'Fasta (tar)');
  } elsif ($prefix eq 'GENE-TREES-TRANSCRIPTS') {
    return ('Gene trees', 'XML (tar)');
  } elsif ($prefix eq 'GENE-TREES-PEPTIDES') {
    return ('Gene trees', 'XML (tar)');
  } elsif ($prefix eq 'HOMOLOGS') {
    return ('Homologs', 'XML (tar)');
  } elsif ($prefix eq 'WG-ALIGN') {
    return ('Whole genome alignments', 'MAF (tar)');
  }
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

sub compara_description {
  my ($self, $file) = @_;
  
  my $drupal_desc = $self->param_required('drupal_desc');
  my ($prefix) = $file =~ /^([^_]+)/;
  my $description = $$drupal_desc{$prefix};
  
  if ($prefix eq 'WG-ALIGN') {
    $file =~ /^$prefix\_[\w\-]+_(\w+)\-(\w+)\-(\w+)_(\w+)\-(\w+)\-([a-zA-Z0-9]+)/;
    my $species1 = "$1 $2 ($3)";
    my $species2 = "$4 $5 ($6)";
    
    $description =~ s/<SPECIES1>/$species1/;
    $description =~ s/<SPECIES2>/$species2/;
  }
  
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
  
  my $xgrid = '0';
  if ($dump_type =~ /fasta/) {
    $xgrid = '1';
  }
  return $xgrid;
}

1;
