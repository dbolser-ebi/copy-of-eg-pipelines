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

Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::AnalysisSetup

=head1 DESCRIPTION



=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::BlastAlignment::WormBase::AnalysisSetup;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::EGPipeline::Common::RunnableDB::AnalysisSetup);
use Bio::EnsEMBL::Analysis;

sub run {
  my $self = shift @_;
  my $species = $self->param_required('species');
  my $logic_name = $self->param_required('logic_name');
  
  my $dba = $self->get_DBAdaptor($self->param('db_type'));
  my $dbh = $dba->dbc->db_handle;
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  
  if (defined $analysis) {

    if ($self->skip($analysis,$aa,$dbh)){
      $dba->dbc->disconnect_if_idle(); 
      return; 
    }

    $self->index_if_needed($analysis);

    if ($self->param('delete_existing')) {
      my $analysis_id = $analysis->dbID;
      foreach my $table (@{$self->param('linked_tables')}) {
        my $sql = "DELETE FROM $table WHERE analysis_id = $analysis_id";
        my $sth = $dbh->prepare($sql) or throw("Failed to delete rows using '$sql': ".$dbh->errstr);
        $sth->execute or throw("Failed to delete rows using '$sql': ".$sth->errstr);
      }
    } else {
      my $logic_rename = $self->param_required('logic_rename');
      my $renamed_analysis = $aa->fetch_by_logic_name($logic_rename);
      if (defined $renamed_analysis) {
        $self->throw(
          "Cannot rename '$logic_name' to '$logic_rename' because '$logic_rename' already exists.\n".
          "Either provide a different 'logic_rename' parameter, or delete the '$logic_rename' analysis.\n"
        );
      } else {
        $analysis->logic_name($logic_rename);
        $aa->update($analysis);
        my $new_analysis = $self->create_analysis;
        $new_analysis->db_version((stat($new_analysis->db_file))[9]);
        $aa->store($new_analysis);
      }
    }
  } else {

    my $new_analysis = $self->create_analysis;
    $new_analysis->db_version((stat($new_analysis->db_file))[9]);
    $aa->store($new_analysis);
  }
 
  if ($self->param('production_lookup')) {
    $self->production_updates;
  }
 
  $dba->dbc->disconnect_if_idle(); 
}

# run the makeblastdb if needed
sub index_if_needed{
   my ($self,$analysis)=@_;
   unless (-e $analysis->db_file . '.phr'){
       my $cmd = $self->param_required('makeblastdb_exe') . ' -in '. $analysis->db_file . ' -dbtype prot';
       my $out = `$cmd 2>&1`;
       $self->throw("Error when executing $cmd:\n$out\n") if $out =~/error/mi;
   }
}

# compare timestamps to determine if the analysis needs to be run
sub skip {

  my ($self,$analysis,$aa,$dbh)=@_;
  my $analysis_timestamp = $analysis->db_version;
  my $blast_timestamp = (stat($analysis->db_file))[9];

  if ($blast_timestamp != $analysis_timestamp) {
    $analysis->db_version($blast_timestamp);
    $aa->update($analysis);
    return undef;
  } else {

    my $analysis_id = $analysis->dbID;
    my $count=0;

    foreach my $table (@{$self->param('linked_tables')}) { # technically there should be only one linked_table
      my $sql = "SELECT COUNT(*) FROM $table WHERE analysis_id = $analysis_id";
      my $sth = $dbh->prepare($sql) or throw("Failed to count rows using '$sql': ".$dbh->errstr);
      $sth->execute or throw("Failed to count rows using '$sql': ".$sth->errstr);
      my($c)=$sth->fetchrow_array;
      # print $self->param('species'),' : ',$analysis->logic_name,"/",$self->param('logic_name')," query:$sql count:$c\n";
      $count+=$c;
      $sth->finish();
    }
    return undef unless $count; # don't skip it if there are zero features
  }
  $self->{skipped_analysis} = $self->param('logic_name');
  return 1;
}

# send skippable analysis to #3 as dead end
sub write_output {
  my ($self) = @_;
  if ($self->{skipped_analysis} eq $self->param('logic_name')) {
     # print "skipping ",$self->param('logic_name')," ",$self->{skipped_analysis},"\n";
     $self->dataflow_output_id($self->{'skipped_analysis'},3);
  }else{
     # print "continuing ",$self->param('logic_name'),"\n";
     $self->dataflow_output_id({},2);
  }
}

1;
