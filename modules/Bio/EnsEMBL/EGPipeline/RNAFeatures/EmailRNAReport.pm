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

Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::EmailRNAReport;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport');

use MIME::Lite;

sub fetch_input {
  my ($self) = @_;
  my $run_cmscan    = $self->param_required('run_cmscan');
  my $run_trnascan  = $self->param_required('run_trnascan');
  my $pipeline_dir  = $self->param_required('pipeline_dir');
  my $evalue_levels = $self->param_required('evalue_levels');
  
  my $text = "The RNA Features pipeline has completed.\n";
  
  if ($run_cmscan) {
    my @evalue_levels = keys %$evalue_levels;
    my $cmscan_file = "$pipeline_dir/cmscan.txt";
    $text .= "\nThe cmscan alignments are summarised in the table below, and in the attached plots.";
    $text .= $self->summarise($cmscan_file, 'E-value', \@evalue_levels);
    
  }
  
  if ($run_trnascan) {
    my $trnascan_file = "$pipeline_dir/trnascan.txt";
    $text .= "\nThe trnascan alignments are summarised in the table below.";
    $text .= $self->summarise($trnascan_file, 'Score', [0, 40]);
  }
  
  $self->param('text', $text);
}

sub run {
  my ($self) = @_;
  my $email         = $self->param_required('email');
  my $subject       = $self->param_required('subject');
  my $text          = $self->param_required('text');
  my $pipeline_dir  = $self->param_required('pipeline_dir');
  my $evalue_levels = $self->param_required('evalue_levels');
  
  my $msg = MIME::Lite->new(
    From    => $email,
    To      => $email,
    Subject => $subject,
    Type    => 'multipart/mixed',
  );
  
  $msg->attach(
    Type => 'TEXT',
    Data => $text,
  );
  
  foreach my $evalue (keys %$evalue_levels) {
    my $biotypes_file = "biotypes_$evalue.svg";
    my $biotypes_path = "$pipeline_dir/$biotypes_file";
    my $distinct_file = "distinct_$evalue.svg";
    my $distinct_path = "$pipeline_dir/$distinct_file";
    
    if (-e $biotypes_path) {
      $msg->attach(
        Type        => 'image/svg+xml',
        Path        => $biotypes_path,
        Filename    => $biotypes_file,
        Disposition => 'attachment',
      );
    }
    
    if (-e $distinct_path) {
      $msg->attach(
        Type        => 'image/svg+xml',
        Path        => $distinct_path,
        Filename    => $distinct_file,
        Disposition => 'attachment',
      );
    }
  }
  
  $msg->send;
}

sub summarise {
  my ($self, $file, $value_type, $thresholds) = @_;
  
  local $/;
  open(my $fh, '<', $file) or $self->throw("Failed to open $file: $!");
  my @rows = split(/\n/, <$fh>);
  close($fh);

  my %summary;
  foreach my $row (@rows) {
    my ($value, $species, $biotype, $name) = split(/\t/, $row);
    foreach my $threshold (@$thresholds) {
      if ($value_type eq 'E-value') {
        if ($value <= $threshold) {
          $summary{$species}{$threshold}{$biotype}{$name}++;
        }
      } else {
        if ($value >= $threshold) {
          $summary{$species}{$threshold}{$biotype}{$name}++;
        }
      }
    }
  }
  
  my @columns = ('Species', $value_type, '# Biotypes', '# Models', '# Alignments');
  my @results;
  foreach my $species (sort keys %summary) {
    foreach my $threshold (sort keys %{$summary{$species}}) {
      my ($biotypes, $models, $alignments) = (0, 0, 0);
      
      foreach my $biotype (keys %{$summary{$species}{$threshold}}) {
        $biotypes++;
        foreach my $name (keys %{$summary{$species}{$threshold}{$biotype}}) {
          $models++;
          $alignments += $summary{$species}{$threshold}{$biotype}{$name};
        }
      }
      
      push @results, [$species, $threshold, $biotypes, $models, $alignments];
    }
  }
  
  return $self->format_table('', \@columns, \@results);
}

1;
