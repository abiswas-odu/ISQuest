#!/usr/bin/perl

use strict;
use warnings;
use File::Spec;
use Bio::Perl;
use File::Basename;

my ($fileName1, $fileName2,$baseSeqFileName) = @ARGV;
my $numSequencesFile1 = numSequences($fileName1);
my $reportEveryFile1 = int($numSequencesFile1/100) || 1;
print "$numSequencesFile1 sequences to convert in $fileName1\n";

my $numSequencesFile2 = numSequences($fileName2);
my $reportEveryFile2=int($numSequencesFile2/100) || 1;
print "$numSequencesFile2 sequences to convert in $fileName2\n";

my $in1=Bio::SeqIO->new(-file=>$fileName1,-format=>"fastq");
my $in2=Bio::SeqIO->new(-file=>$fileName2,-format=>"fastq");

#my $baseSeqFileName = basename($fileName1);
#$baseSeqFileName = substr($baseSeqFileName,0,index($baseSeqFileName,'.'));
$baseSeqFileName = $baseSeqFileName . '.fasta';

my $seqOut=Bio::SeqIO->new(-file=>">$baseSeqFileName",-format=>"fasta");

my $seqCount=0;
my $percentDone=0;
while(my $seq=$in1->next_seq){
    my $fastaSeqID=$seq->id;
    if(defined $seq->desc){
    	my $fastaSeqID=$seq->id . '_' . $seq->desc;
    }
    $fastaSeqID =~ s/\s//g;
    my $fastaSeq=Bio::PrimarySeq->new(-seq => $seq->seq(), -display_id => $fastaSeqID);
    $seqOut->write_seq($fastaSeq);
	$seqCount++;
    if($seqCount%$reportEveryFile1 == 0){
      $percentDone++;
      print "$percentDone%..";
    }
}
print "Done with subfile $fileName1.\n";
$seqCount=0;
$percentDone=0;
while(my $seq=$in2->next_seq){
    my $fastaSeqID=$seq->id;
    if(defined $seq->desc){
       my $fastaSeqID=$seq->id . '_' . $seq->desc;
    }
    $fastaSeqID =~ s/\s//g;
    my $fastaSeq=Bio::PrimarySeq->new(-seq => $seq->seq(), -display_id => $fastaSeqID);
    $seqOut->write_seq($fastaSeq);
	$seqCount++;
    if($seqCount%$reportEveryFile2 == 0){
      $percentDone++;
      print "$percentDone%..";
    }
}
print "Done with subfile $fileName2.\n";


sub numSequences{
  my ($file) = @_;
  my $num = `wc -l $file | cut -f1 -d' '`;
  chomp($num);
  $num = $num / 4;
  return $num;
}
