#!/usr/bin/env perl

#Created by Sarah Schmedes and August Woerner
#Will run metaphlan2 on a group of desginated SRR files
#The methaphlan2 run will group all SRR files by SRS number
#to create one output per sample
#Date Released: 08/01/2017

use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

#select samples to run through metaphlan2
my $metafile = $OH_METAFILE;
my ($bodysite, $QCdatapath, $outdir);

#Set flags
GetOptions('bodysite=s' => \$bodysite,
	   'input=s' => \$QCdatapath,
	   'output=s' => \$outdir)
    or die " Must declar -bodysite flag and bodysite symbol to run script\n";

#collect all SRSs for a particular body site
my @SRSs = queryDB($metafile, $bodysite, "bodysiteID", "srsID");
my $SRSs;

foreach my $sample (@SRSs)
{
    #Collect all SRR numbers for each SRS sample, one at a time
    my ($metainput, @metainput, $filesin);
    
    my @sampleSRRs = queryDB($metafile, $sample, "srsID", "srrID");


    #create input variable for each SRR for metaphlan2 command
    $filesin='';
    foreach (@sampleSRRs)
    {
	$metainput = "$QCdatapath/$_\.1.qc.fastq.bz2 $QCdatapath/$_\.2.qc.fastq.bz2";
	$filesin .= ' ' . $metainput;
    }

#run metaphlan2 on all samples
#all metaphlan2/strainphlan scripts must be in the environmental path
    if (system("bash", "-c", "metaphlan2.py --input_type multifastq <(bzcat $filesin) $outdir/$sample\_profile.txt --bowtie2out $outdir/$sample\_bowtie2.txt --samout $outdir/$sample\.sam.bz2 --nproc 20"))
    {
	die "Metaphlan2 ERROR: $!";
    }
    
}









