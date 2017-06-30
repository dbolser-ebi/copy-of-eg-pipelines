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

Bio::EnsEMBL::EGPipeline::CoreStatistics::EmailSummary

=head1 DESCRIPTION

Collate job data then email it.
Adapted from Bio::EnsEMBL::Production::Pipeline::Production::EmailSummaryCore.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::CoreStatistics::EmailSummary;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Hive::RunnableDB::NotifyByEmail Bio::EnsEMBL::Production::Pipeline::Common::Base);

use Bio::EnsEMBL::Hive::Utils qw/destringify/;

sub fetch_input {
  my ($self) = @_;
  
  $self->assert_executable('sendmail');
  
  my $can_transcripts = $self->jobs('CanonicalTranscripts');
  my $ct_exons = $self->jobs('ConstitutiveExons');
  my $gene_count = $self->jobs('GeneCount');
  my $gene_gc = $self->jobs('GeneGC');
  my $genome_stats = $self->jobs('GenomeStats');
  my $meta_coords = $self->jobs('MetaCoords');
  my $meta_levels = $self->jobs('MetaLevels');
  my $pep_stats = $self->jobs('PepStats');
  my $coding_density = $self->jobs('CodingDensity');
  my $pseudogene_density = $self->jobs('PseudogeneDensity');
  my $short_non_coding_density = $self->jobs('ShortNonCodingDensity');
  my $long_non_coding_density = $self->jobs('LongNonCodingDensity');
  my $percent_gc = $self->jobs('PercentGC');
  my $percent_repeat = $self->jobs('PercentRepeat');
  my $snp_density = $self->jobs('SnpDensity');
  my $snp_count = $self->jobs('SnpCount');
  my $analyze_tables = $self->jobs('AnalyzeTables');

    
  my @args = (
    $can_transcripts->{successful_jobs},
    $can_transcripts->{failed_jobs},
    $ct_exons->{successful_jobs},
    $ct_exons->{failed_jobs},
    $gene_count->{successful_jobs},
    $gene_count->{failed_jobs},
    $gene_gc->{successful_jobs},
    $gene_gc->{failed_jobs},
    $genome_stats->{successful_jobs},
    $genome_stats->{failed_jobs},
    $meta_coords->{successful_jobs},
    $meta_coords->{failed_jobs},
    $meta_levels->{successful_jobs},
    $meta_levels->{failed_jobs},
    $pep_stats->{successful_jobs},
    $pep_stats->{failed_jobs},
    $coding_density->{successful_jobs},
    $coding_density->{failed_jobs},
    $pseudogene_density->{successful_jobs},
    $pseudogene_density->{failed_jobs},
    $short_non_coding_density->{successful_jobs},
    $short_non_coding_density->{failed_jobs},
    $long_non_coding_density->{successful_jobs},
    $long_non_coding_density->{failed_jobs},
    $percent_gc->{successful_jobs},
    $percent_gc->{failed_jobs},
    $percent_repeat->{successful_jobs},
    $percent_repeat->{failed_jobs},
    $snp_density->{successful_jobs},
    $snp_density->{failed_jobs},
    $snp_count->{successful_jobs},
    $snp_count->{failed_jobs},
    $analyze_tables->{successful_jobs},
    $analyze_tables->{failed_jobs},
    $self->failed(),
    $self->summary($can_transcripts),
    $self->summary($ct_exons),
    $self->summary($gene_count),
    $self->summary($gene_gc),
    $self->summary($genome_stats),
    $self->summary($pep_stats),
    $self->summary($meta_coords),
    $self->summary($meta_levels),
    $self->summary($coding_density),
    $self->summary($pseudogene_density),
    $self->summary($short_non_coding_density),
    $self->summary($long_non_coding_density),
    $self->summary($percent_gc),
    $self->summary($percent_repeat),
    $self->summary($snp_density),
    $self->summary($snp_count),
    $self->summary($analyze_tables),
  );
  
  my $msg = sprintf(<<'MSG', @args);
Your Production Pipeline has finished. We have:

  * %d species with canonical transcripts (%d failed)
  * %d species with constitutive exons (%d failed)
  * %d species with gene count (%d failed)
  * %d species with gene gc (%d failed)
  * %d species with genome_stats attributes (%d failed)
  * %d species with meta coords (%d failed)
  * %d species with meta levels (%d failed)
  * %d species with pep stats (%d failed)
  * %d species with coding density (%d failed)
  * %d species with pseudogene density (%d failed)
  * %d species with short noncoding density (%d failed)
  * %d species with long noncoding density (%d failed)
  * %d species with percent gc (%d failed)
  * %d species with percent repeat (%d failed)
  * %d species with snp density (%d failed)
  * %d species with snp count (%d failed)

  * %d species with analyzed tables (%d failed)

%s

===============================================================================

Full breakdown follows ...

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

%s

MSG
  $self->param('text', $msg);
  return;
}

sub jobs {
  my ($self, $logic_name) = @_;
  my $aa = $self->db->get_AnalysisAdaptor();
  my $aja = $self->db->get_AnalysisJobAdaptor();
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  my @jobs;
  if (!$analysis) {
    return {
      name => $logic_name,
      successful_jobs => 0,
      failed_jobs => 0,
      jobs => \@jobs,
    };
  }
  my $id = $analysis->dbID();
  @jobs = @{$aja->fetch_all_by_analysis_id($id)};
  $_->{input} = destringify($_->input_id()) for @jobs;
  @jobs = sort { $a->{input}->{species} cmp $b->{input}->{species} } @jobs;
  my %passed_species = map { $_->{input}->{species}, 1 } grep { $_->status() eq 'DONE' } @jobs;
  my %failed_species = map { $_->{input}->{species}, 1 } grep { $_->status() eq 'FAILED' } @jobs;
  return {
    analysis => $analysis,
    name => $logic_name,
    jobs => \@jobs,
    successful_jobs => scalar(keys %passed_species),
    failed_jobs => scalar(keys %failed_species),
  };
}


sub failed {
  my ($self) = @_;
  my $failed = $self->db()->get_AnalysisJobAdaptor()->fetch_all_by_analysis_id_status(undef, 'FAILED');
  if(! @{$failed}) {
    return 'No jobs failed. Congratulations!';
  }
  my $output = <<'MSG';
The following jobs have failed during this run. Please check your hive's error msg table for the following jobs:

MSG
  foreach my $job (@{$failed}) {
    my $analysis = $self->db()->get_AnalysisAdaptor()->fetch_by_dbID($job->analysis_id());
    my $line = sprintf(q{  * job_id=%d %s(%5d) input_id='%s'}, $job->dbID(), $analysis->logic_name(), $analysis->dbID(), $job->input_id());
    $output .= $line;
    $output .= "\n";
  }
  return $output;
}

my $sorter = sub {
  my $status_to_int = sub {
    my ($v) = @_;
    return ($v->status() eq 'FAILED') ? 0 : 1;
  };
  my $status_sort = $status_to_int->($a) <=> $status_to_int->($b);
  return $status_sort if $status_sort != 0;
  return $a->{input}->{species} cmp $b->{input}->{species};
};

sub summary {
  my ($self, $data) = @_;
  my $name = $data->{name};
  my $underline = '~'x(length($name));
  my $output = "$name\n$underline\n\n";
  my @jobs = @{$data->{jobs}};
  if(@jobs) {
    foreach my $job (sort $sorter @{$data->{jobs}}) {
      my $species = $job->{input}->{species};
      $output .= sprintf("  * %s - job_id=%d %s\n", $species, $job->dbID(), $job->status());
    }
  }
  else {
    $output .= "No jobs run for this analysis\n";
  }
  $output .= "\n";
  return $output;
}

1;
