=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::CMScan

=head1 SYNOPSIS

  my $cmscan = Bio::EnsEMBL::Analysis::Runnable::CMScan->
  new(
      -query => $slice,
      -program => 'cmscan',
     );
  $cmscan->run;
  my @features = @{$cmscan->output};

=head1 DESCRIPTION

CMScan uses the cmscan program in the Infernal 1.1 suite to annotate RNA genes 
on a given chunk of sequence. The source of the covariance models will probably 
be from Rfam (ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz), but
any set of covariance models will work, as long as they conform to the Infernal
specification (http://selab.janelia.org/software/infernal/Userguide.pdf).

=cut

package Bio::EnsEMBL::Analysis::Runnable::CMScan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

use List::Util qw(min max);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CMScan
  Arg [2]   : string, cmpress:     path to the 'cmpress' executable
  Arg [3]   : string, cm_file:     path to the file of covariance models
  Arg [4]   : int, cpu:            number of processors to use
  Arg [5]   : string, heuristics:  "slowest", "slower", "slow", "default", "faster", or "fastest"
  Arg [6]   : string, threshold:   E-value threshold
  Arg [7]   : string, db_name:     db_name from the external_db table
  Arg [8]   : int, recreate_index: overwrite existing indexes
  Function  : create a new  Bio::EnsEMBL::Analysis::Runnable::CMScan
  Returntype: Bio::EnsEMBL::Analysis::Runnable::CMScan
  Exceptions: 
  Example   : 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($cmpress, $cm_file, $cpu, $heuristics, $threshold, $db_name, $recreate_index) =
    rearrange(['CMPRESS', 'CM_FILE', 'CPU', 'HEURISTICS', 'THRESHOLD', 'DB_NAME', 'RECREATE_INDEX'], @args);
  
  $self->program('cmscan') if (!$self->program);
  
  ($cmpress = $self->program('cmscan')) =~ s/cmscan/cmpress/ unless defined $cmpress;
  $self->cmpress($cmpress);
  
  if (defined $cm_file) {
    $self->cm_file($cm_file);
  } else {
    throw('A file with covariance models (-cm_file) is required.');
  }
  
  # By default, cmscan will grab 1 more CPU than you ask for...
  $cpu = 1 unless defined $cpu;
  $self->cpu($cpu - 1);
  
  $heuristics = 'default' unless defined $heuristics;
  $self->heuristics($heuristics);
  
  $threshold = '0.001' unless defined $threshold;
  $self->threshold($threshold);
  
  $self->db_name($db_name);
  
  $self->recreate_index($recreate_index);
  
  return $self;
}

sub cmpress {
  my $self = shift;
  $self->{'cmpress'} = shift if (@_);
  return $self->{'cmpress'};
}

sub cm_file {
  my $self = shift;
  $self->{'cm_file'} = shift if (@_);
  return $self->{'cm_file'};
}

sub cpu {
  my $self = shift;
  $self->{'cpu'} = shift if (@_);
  return $self->{'cpu'};
}

sub heuristics {
  my $self = shift;
  $self->{'heuristics'} = shift if (@_);
  return $self->{'heuristics'};
}

sub threshold {
  my $self = shift;
  $self->{'threshold'} = shift if (@_);
  return $self->{'threshold'};
}

sub db_name {
  my $self = shift;
  $self->{'db_name'} = shift if (@_);
  return $self->{'db_name'};
}

sub recreate_index {
  my $self = shift;
  $self->{'recreate_index'} = shift if (@_);
  return $self->{'recreate_index'};
}

=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CMScan
  Arg [2]   : string, program name
  Function  : generates and executes the cmscan command
  Returntype: none
  Exceptions: throws if execution failed
  Example   : 

=cut

sub run_analysis {
  my ($self, $program) = @_;
  if (!$program) {
    $program = $self->program;
  }
  throw("$program is not executable") unless($program && -x $program);
  
  $self->prepare_index();
  
  my $options = $self->options;
  $options .= " --notextw ";
  
  $options .= ' --cpu '.$self->cpu.' ';
  $options .= ' --incE '.$self->threshold.' ';
  
  my $heuristics = $self->heuristics;
  if ($heuristics eq 'slowest') {
    $options .= ' --max ';
  } elsif ($heuristics eq 'slower') {
    $options .= ' --nohmm ';
  } elsif ($heuristics eq 'slow') {
    $options .= ' --mid ';
  } elsif ($heuristics eq 'faster') {
    $options .= ' --rfam ';
  } elsif ($heuristics eq 'fastest') {
    $options .= ' --hmmonly ';
  } else {
    $options .= ' --default ';
  }
  
  $self->options($options);
  
  my $cm_file  = $self->cm_file;
  my $in_file  = $self->queryfile;
  my $out_file = $self->queryfile.'.out';
  my $aln_file = $self->queryfile.'.aln';
  $self->resultsfile($out_file);
  
  my $cmd = "$program $options --tblout $out_file -o $aln_file $cm_file $in_file";
  system($cmd) == 0 or throw("Failed to run: $cmd");
  
  my $cat_cmd = "cat $aln_file >> $out_file";
  system($cat_cmd) == 0 or throw("Failed to run: $cat_cmd");
}

sub prepare_index {
  my ($self, $program) = @_;
  
  my $cm_file  = $self->cm_file;
  
  my @suffixes = qw( i1f i1i i1m i1p );
  my @indexes = map { "$cm_file.$_" } @suffixes;
  my $missing = 0;
  foreach my $index (@indexes) {
    $missing++ unless -e $index;
  }
  
  if ($missing || $self->recreate_index) {
    if (!$program) {
      $program = $self->cmpress;
    }
    throw("$program is not executable") unless($program && -x $program);
    
    my $cmd = "$program $cm_file";
    system($cmd) == 0 or throw("Failed to run: $cmd");
  }
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CMScan
  Arg [2]   : string, filename
  Function  : open and parse the results file into sequence alignments
  Returntype: none 
  Exceptions: throws on failure to open or close output file
  Example   : 

=cut

sub parse_results {
  my ($self, $results) = @_;
  if (!$results) {
    $results = $self->resultsfile;
  }
  
  # There is an easy-to-parse table output, but it lacks the alignment, which
  # is required to add gaps for the alignment feature (as a cigar line).
  # Unfortunately, until the next release of Infernal we also need to parse
  # the table output, which is prepended to the results, in order to get both
  # the Rfam name and accession...
  
  open(my $fh, $results) or throw("Failed to open results file '$results'");
  my $results_text = do { local $/; <$fh> };
  close($fh) or throw("Failed to close results file '$results'");
  
  $results_text =~ s/^#.*\n(^[^#]+)^#.*\n//m;
  my %names = $1 =~ /^(\S+)\s+(\S+)/gm;
  
  my @lines = split(/\n/, $results_text);
  
  my @results;
  my $current;
  
  foreach my $line (@lines) {
    if ($line =~ /^\S/) {
      push(@results, $current) if $current;
      if ($line =~ /^>>/) {
        $current = "$line\n";
      } else {
        $current = undef;
      }
    } elsif ($current) {
      $current .= "$line\n";
    }
  }
  push(@results, $current) if $current;
  
  my @features;
  
  foreach my $result (@results) {
    my $hit = $self->parse_result($result, \%names);
    
    if (defined $hit) {
      my $feature = Bio::EnsEMBL::DnaDnaAlignFeature->new(
        -slice            => $self->query,
        -external_db_name => $self->db_name,
        %$hit,
      );
      push (@features, $feature);
    }
  }
  
  $self->output(\@features);
}

sub parse_result {
  my ($self, $result, $names) = @_;
  
  my ($rna_name, $rna_desc) = $result =~ /^>>\s+(\S+)[ ]+(.*)/;
  my $rna_acc = $$names{$rna_name};
  my $biotype = $self->biotype($rna_name, $rna_desc);
  
  my ($stats) = $result =~ /^\s+\(\d+\)\s+(.+)/m;
  
  my
  ( $significant, $evalue, $score, $bias,
    $model, $model_start, $model_end, undef,
    $start, $end, $strand, undef, $accuracy, $trunc, $gc
  ) = split(/\s+/, $stats);
  
  if ($significant eq '!') {
    my ($structure) = $result =~ /^\s+(\S+)\s+CS$/m;
    my ($subject, $query) = $result =~ /CS\n\s+(.+)\n.+\n\s+(.+)/m;
    my ($subject_seq) = $subject =~ /^\s*\S+\s+\d+\s+(.+)\s+\d+\s*$/;
    my ($query_seq) = $query =~ /^\s*\S+\s+\d+\s+(.+)\s+\d+\s*$/;
    $subject_seq =~ s/\s+$//;
    $query_seq =~ s/\s+$//;
    
    my ($ungapped_structure, $cigar) =
      $self->parse_alignment($structure, $subject_seq, $query_seq);
    
    my %extra_data = (
      'Biotype' => $biotype,
      'Desc' => $rna_desc,
      'Trunc' => $trunc,
      'Accuracy' => $accuracy,
      'Bias' => $bias,
      'GC' => $gc,
      'Significant' => $significant,
      'Structure' => $ungapped_structure,
    );
    $extra_data{'Accession'} = $rna_acc if $rna_acc && $rna_acc ne '-';
    
    my %hit =
    ( 
      -start        => min($start, $end),
      -end          => max($start, $end),
      -strand       => $start > $end ? -1 : 1,
      -hstart       => $model_start,
      -hend         => $model_end,
      -hstrand      => $model_start > $model_end ? -1 : 1,
      -score        => $score,
      -p_value      => $evalue,
      -hseqname     => $rna_name,
      -cigar_string => $cigar,
      -extra_data   => \%extra_data,
    );
  
    return \%hit;
    
  } else {
    return;
    
  }
}

sub parse_alignment {
  my ($self, $structure, $subject_seq, $query_seq) = @_;
  
  my @subject_seq = split(//, $subject_seq);
  my @query_seq   = split(//, $query_seq);
  
  # To make parsing a bit easier, project the deletions, from the
  # subject sequence to the query sequence...
  for (my $i=0; $i<length($subject_seq); $i++) {
    if ($subject_seq[$i] eq '.') {
      $query_seq[$i] = '.';
    }
  }

  $subject_seq = join('', @subject_seq);
  $query_seq   = join('', @query_seq);
  
  # ... then expand long insertions and deletions.
  my @delete_lengths = $subject_seq =~ /\*\[\s*(\d*)\]\*/g;
  my @insert_lengths = $query_seq =~ /\*\[\s*(\d*)\]\*/g;
  
  my $start_q = $query_seq;
  my $start_struc = $structure;
  
  for (my $i=0; $i<scalar(@delete_lengths); $i++) {
    my $delete = $delete_lengths[$i] ? '.' x $delete_lengths[$i] : '';
    my $insert = $insert_lengths[$i] ? '-' x $insert_lengths[$i] : '';
    
    $query_seq =~ s/(\*\[\s*(\d+)\]\*)/$insert$delete/;
    my $sub_structure = '~' x length($1);
    $structure =~ s/$sub_structure/$insert$delete/;
  }
  
  @query_seq = split(//, $query_seq);
  my @structure = split(//, $structure);
  
  # Keep track of whether we're in a match, insert, or delete state; as long
  # as we stay in that state, increment the counter. When the state changes,
  # write out the current state, then reset the state and counter.
  my @ungapped_structure;
  my ($cigar, $state, $count) = ('', 'M', 0);
  
  for (my $i=0; $i<length($query_seq); $i++) {
    if ($query_seq[$i] eq '-') {
      if ($state eq 'I') {
        $count++;
      } else {
        $cigar .= $count if $count > 1;
        $cigar .= $state;
        $state  = 'I';
        $count  = 1;
      }
    } else {
      push @ungapped_structure, $structure[$i];
      
      if ($query_seq[$i] eq '.') {
        if ($state eq 'D') {
          $count++;
        } else {
          $cigar .= $count if $count > 1;
          $cigar .= $state;
          $state  = 'D';
          $count  = 1;
        }
      } else {
        if ($state eq 'M') {
          $count++;
        } else {
          $cigar .= $count if $count > 1;
          $cigar .= $state;
          $state  = 'M';
          $count  = 1;
        }
      }
    }
  }
  $cigar .= $count if $count > 1;
  $cigar .= $state;
  
  my $ungapped_structure = join("", @ungapped_structure);
  
  return ($ungapped_structure, $cigar);
}

sub biotype {
  my ($self, $rna_name, $rna_desc) = @_;
  
  my $biotype;
  
  # Note that the order in the following conditional is important, as
  # it goes from the specific to the generic (e.g. RNase_P is a ribozyme,
  # but it's sufficiently noteworthy to have its own biotype).
  
  if ($rna_name =~ /^class_I_RNA/) {
    $biotype = 'class_I_RNA';
    
  } elsif ($rna_name =~ /^class_I_RNA/) {
    $biotype = 'class_II_RNA';
    
  } elsif ($rna_name =~ /(^IRES_|_IRES$)/) {
    $biotype = 'misc_RNA';
    
  } elsif ($rna_name =~ /_SRP$/) {
    $biotype = 'SRP_RNA';
    
  } elsif ($rna_name =~ /_tmRNA$/) {
    $biotype = 'tmRNA';
    
  } elsif ($rna_name =~ /^tRNA/) {
    $biotype = 'tRNA';
    
  } elsif ($rna_name =~ /Y_RNA$/) {
    $biotype = 'Y_RNA';
    
  } elsif ($rna_desc =~ /antisense RNA/i) {
    $biotype = 'antisense_RNA';
    
  } elsif ($rna_desc =~ /catalytic intron/i) {
    $biotype = 'sense_intronic';
    
  } elsif ($rna_desc =~ /(conserved region|long ncRNA)/i) {
    $biotype = 'lncRNA'; # Need to add this biotype for transcripts
    
  } elsif ($rna_desc =~ /(exon|transcript) \d+$/i) {
    $biotype = 'lncRNA';
    
  } elsif ($rna_desc =~ /CRISPR/i) {
    $biotype = 'CRISPR';
    
  } elsif ($rna_desc =~ /microRNA/i) {
    $biotype = 'pre_miRNA';
    
  } elsif ($rna_desc =~ /ribosomal RNA/i) {
    $biotype = 'rRNA';
    
  } elsif ($rna_desc =~ /RNase[ _]*MRP/i) {
    $biotype = 'RNase_MRP_RNA';
    
  } elsif ($rna_desc =~ /RNase[ _]*P/i) {
    $biotype = 'RNase_P_RNA';
    
  } elsif ($rna_desc =~ /ribozyme/i) {
    $biotype = 'ribozyme';
    
  } elsif ($rna_desc =~ /small Cajal/i) {
    $biotype = 'scaRNA';
    
  } elsif ($rna_desc =~ /(small nucleolar RNA|snoRNA)/i) {
    $biotype = 'snoRNA';
    
  } elsif ($rna_desc =~ /(small nuclear RNA|spliceosomal RNA)/i) {
    $biotype = 'snRNA';
    
  } elsif ($rna_desc =~ /(sRNA|small RNA)/i) {
    $biotype = 'sRNA';
    
  } elsif ($rna_desc =~ /telomerase/i) {
    $biotype = 'telomerase_RNA';
    
  } elsif ($rna_desc =~ /^vault RNA/i) {
    $biotype = 'vaultRNA';
    
  } elsif ($rna_desc =~ /ncRNA/) {
    $biotype = 'ncRNA';
    
  } elsif ($rna_desc =~ /(element|frameshift|hairpin|insertion sequence|leader|motif|pseudoknot|riboswitch|signal|stem[ \-]loop)/i) {
    $biotype = 'misc_RNA';
    
  } else {
    $biotype = 'ncRNA';
  }
  
  return $biotype;
}

1;
