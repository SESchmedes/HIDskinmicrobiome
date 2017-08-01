#!/usr/bin/env perl

#Authored by Sarah E. Schmedes
#Date Released: 08/01/2017
use strict;
use warnings;
use Genomix qw(:constants :meta); #imports genomic constants and meta information
use Getopt::Long;

print $SAMTOOLS , "\nIs the current samtools!\n";
my ($bodysite, $Ohdatapath, $qcOUTpath, $OhqcFastqs);
#Set flags
GetOptions('bodysite=s' => \$bodysite,
	   'fastq_input=s' => \$Ohdatapath,
	   'intermediate_path=s' => \$qcOUTpath, 
	   'qc_outpath=s' => \$OhqcFastqs)
    or die "Bodysite flag is required for run\n";

#Select files to run based on criteria in OhMeta
my $metafile = $OH_METAFILE;
my @SRRs = queryDB($metafile, $bodysite, "bodysiteID", "srrID"); #comment out when not using these samples
my $SRRs;

if (! @SRRs)
{
    die "Sample Array did not form!\n";
}

foreach $SRRs (@SRRs)
{
    #cutadapt parameters
    my $Nextera3a = "CTGTCTCTTATA"; #comment out if not using Nextera
    my $Nextera5a = "TATAAGAGACAG"; #comment out if not using Nextera
    my $qual = 20; #can change quality score threshold
    my $length = 50; #can change length filter cutoff

    print "Cutadapt running on $SRRs\n";

    if (system ("cutadapt -a $Nextera3a -A $Nextera3a -g $Nextera5a -G $Nextera5a -n 2 --trim-n -q $qual -m $length -o $qcOUTpath/$SRRs\.1.fastq -p $qcOUTpath/$SRRs\.2.fastq $Ohdatapath/$SRRs\_1.fastq.bz2 $Ohdatapath/$SRRs\_2.fastq.bz2 > $qcOUTpath/$bodysite\_cutadapt_report.txt"))
    {
	die "Cutadapt failed\nError is: $!\n";
    }

    if (! -f "$qcOUTpath/$SRRs\.1.fastq" || ! -f "$qcOUTpath/$SRRs\.2.fastq")
    {
	die "Cutadapt failed\nError is: $!\n";
    }
}

print "All SRR IDs completed processing through cutadapt!\n";
#die "cutadapt testing complete";

foreach $SRRs (@SRRs)
{
    #Run BWA to align human reads to hg38 genome, pipe sam file to samtools (sam file not saved)
    #samtools to produce human and microbial bam file 
    #convert sam to bam file and output alignments with flag 13(read paired, unmapped, mate unmapped) 
    #output is uncompressed bam, input is sam, include header in output
    #all other alignments should be output by -U flag (should I do samtools twice to specify out?)
    if( system ("$BWA mem -M /mnt/blsm/PublicData/HG38/chrAll.standardChroms.fa $qcOUTpath/$SRRs\.1.fastq $qcOUTpath/$SRRs\.2.fastq -t 20 | $SAMTOOLS view - -f 0xD -u -h -o $qcOUTpath/$SRRs\.bam -U $qcOUTpath/$SRRs\.human.bam"))
    {
	die "bwa or samtools failed\nError is: $!\n";
    }

    if (! -f "$qcOUTpath/$SRRs\.bam" || ! -f "$qcOUTpath/$SRRs\.human.bam")
    {
	die "bwa or samtools failed\nError is: $!\n";
    }
}

print "All SRR bam files have been created using BWA and SAMTOOLS!\n";
#die "bwa and samtools testing complete! SUCCESS!!!";

foreach $SRRs (@SRRs)
{
    #samtools to convert microbial.bam to fastq, pulls out paired end reads into separate files


	if (system ("$SAMTOOLS fastq $qcOUTpath/$SRRs\.bam -1 $OhqcFastqs/$SRRs\.1.qc.fastq -2 $OhqcFastqs/$SRRs\.2.qc.fastq") )
	{
	    die "samtools failed\nError is: $!\n";
	}
	
	if (! -f "$OhqcFastqs/$SRRs\.1.qc.fastq" || ! -f "$OhqcFastqs/$SRRs\.2.qc.fastq")
	{
	    die "You do not have fastqs for $SRRs\nError is: $!\n";
	}
	
	if (system ("bzip2 $OhqcFastqs/$SRRs\.1.qc.fastq $OhqcFastqs/$SRRs\.2.qc.fastq") )
	{
	    die "bz2 compression of .qc.fastq files has failed\nError is: $!\n";
	}
	
	
	if (system ("rm $qcOUTpath/$SRRs\.1.fastq $qcOUTpath/$SRRs\.2.fastq")) {
	    die "removal of intermediate fastq files has failed\nError is: $!\n";
	}
	
}


print "Read preprocessing on all designated SRR IDs complete!\n";
print "You now have bz2 compressed, quality-controlled fastq files!\n"; 
