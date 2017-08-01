#!/usr/bin/env perl

############################################################
#Created by Sarah Schmedes
#Script will pull out reads aligning to markers from shared species across samples in a bodysite
#and create mpileup files, identify variants, and calculate pi (ND) across markers
#(parse through an .mpileup file (run with -Os flags))
############################################################
#Date released: 08/01/2017

use strict;
use warnings;
use Getopt::Long;
use Genomix qw(:constants :meta);

my ($bodysite, $samfilepath, $pileuppath, $samplelist, $db_markers, $help);
my $metafile = $OH_METAFILE;
my (%markerND, %markerNDlength, %finalmarkers);
#setflags
GetOptions('bodysite=s' => \$bodysite, #group of SRSs to run for bodysite
	   'sampath=s' => \$samfilepath,
	   'pileuppath=s' => \$pileuppath,
	   'samplelist=s' => \$samplelist,
	   'db_markers=s' => \$db_markers,
	   'help' => \$help)
    or die "Error with flags. Designatre -bodysite or see -help for more information\n";
if ($help) {
    die "This script can take in 4 flags. -bodysite, -sampath, -pileuppath, and -help [to see this message].\n";
}
#Get desired samples to analyze
open SAMPLES, "<$samplelist" or die "Could not open revised sample list for reading! $!\n";
my %samplenames;
while (<SAMPLES>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^Bodysite/) {
	next;
    }
    my ($body, $ind, $time, @data) = split /\t/, $line;
    if ($body eq $bodysite) {
	my $samplename = "$ind\_$body\_$time";
	$samplenames{$samplename} = 1;
    }
}
close(SAMPLES);
my @sortednames = sort {$a cmp $b} keys %samplenames;
my %bodysamplename;
if ($bodysite ne "foot") {
    my @bodySRS = queryDB($metafile, $bodysite, "bodysiteID", "srsID");
    my $samplename_ref = samplename(@bodySRS); #argument is an array of SRS numbers
    %bodysamplename = %$samplename_ref;
}
my @orderedSRS;
#Get an ordered array of SRSs by samplenames
if ($bodysite eq "foot") {
    @orderedSRS = @sortednames;
} else {
    foreach my $samplename (@sortednames) {
	foreach my $key (keys %bodysamplename) {
	    if ($bodysamplename{$key} eq $samplename) {
		push(@orderedSRS, $key);
	    }
	}
    }
}
#Get taxon path and length for each marker
my ($taxa_ref, $taxaLevel_ref, $marker_ref) = gettaxa();
my %taxamarker = %$taxa_ref; #hash of each marker and associated full taxon path
my %taxaLevel = %$taxaLevel_ref; #hash of each full taxon path and associated
                                 #leaf clade tab separated corresponding taxon level
my %markerlength = %$marker_ref; #hash of each marker and associated length
my @totalmarkers = sort keys %taxamarker;
foreach my $file (@orderedSRS) {
    if (! -f "$samfilepath/$file\.sorted.bam") {
	if (system("bzcat $samfilepath/$file\.sam.bz2 | samtools-1.3.1 view -h -u - | samtools-1.3.1 sort -m 10G -o $samfilepath/$file\.sorted.bam -")) {
	    die "Samtools-1.3.1 error creating .sorted.bam file: $!";
	}
    } else {
	next;
    }
}

#Run samtools-1.3.1 mpileup on designated samples
my @files;
foreach my $file (@orderedSRS) {
    push(@files, "$samfilepath/$file\.sorted.bam");
}

my $filestring = join("\t", @files);
if (! -f "$pileuppath/$bodysite\_$analysis\.mpileup") {
    if (system("samtools-1.3.1 mpileup -Os -d 1000 -f $db_markers/all_markers.fasta $filestring > $pileuppath/$bodysite\.mpileup")) {
	die "Samtools-1.3.1 mpileup error: $!";
    }
}

open MPILE, "<$pileuppath/$bodysite\.mpileup" or die "Could not open input file for reading\n";
open OUT, ">$pileuppath/$bodysite\_ND.txt\n" or die "Could not open output file for writing\n";
my $linenum = 1;
my $samplelabel;
my $nsamples = scalar(@orderedSRS);

while (my $record = <MPILE>) {
    my ($site, $total);
    my (@sites, @bases, @coverage);
    my @samplelabels;
    my (%basestring, %bases, %stotal, %totalcount, %samplecount);
    my $lineoutput;
#    warn "$.\n";

    chomp $record;
    #Store reference fields
    my ($marker, $pos, $refbase, @sampledata) = split/\t/, $record;
    $site = join('_', $marker, $pos);
    push(@sites, $site);

    #Determine number of samples in input file/correct any empty string or * for samples
    my $data = scalar(@sampledata);
    if ($data % 5 != 0) {
	if ($sampledata[$#sampledata] == 0) { #correct empty strings
	    for (1..4) {
		push(@sampledata, '');
	    }
	    $data = scalar(@sampledata);
	} else {
	    die "First:There are not 5 data columns per sample!\nLine number $linenum\n$record\n";
	}
    }
    if ($data % 5 != 0) { #double check the corrected empty strings
	die "Second:There are not 5 data columns per sample!\nLine number $linenum\n$record\n";
    }

    if ($linenum == 1) { #Print header line for output
	my $sampleheader;
	foreach my $samplename (@sortednames) {
	    $sampleheader .= "Sample\t$samplename"."ND\t$samplename"."TotalCov\t";
	}
	print OUT "Marker\tPosition\tTaxonPath\tTaxon\tLevel\t$sampleheader\n";
    }
    my $coverage;
    foreach my $sample (@sortednames) { #hash of all bases for each sample
	my @fields = splice (@sampledata, 0, 5);
	if ($fields[0] < 2) {
	    $coverage++;
	}
	$basestring{$sample} = $fields[1];
    }

    if (defined $coverage && ($coverage == $nsamples)) { #skip over any sites that did not align for all samples
	$linenum++;
	next;
    }
    
    #count bases for each sample and across both samples
    foreach my $sample (keys %basestring) {#hash of base strings for each sample
	my $samplebases = $basestring{$sample};
	$samplebases =~ s/[+-](\d+)(??{"[ACGTN]{$1}"})//gi; #removes indels
	my @samplebases = split //, $samplebases;
	foreach my $base (@samplebases) {
	    if ($base eq '.' || $base eq ',') {
		$base = $refbase;
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'A' || $base eq 'a') {
		$base = 'A';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'T' || $base eq 't') {
		$base = 'T';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'C' || $base eq 'c') {
		$base = 'C';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'G' || $base eq 'g') {
		$base = 'G';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    else { #not interested in other characters
		next;
	    }
	}
    }
    my @basecalls = keys %totalcount; #Bases observed across samples
    if (!@basecalls) {
	$linenum++;
	next;
    }

    my $samplecheck = 0;
    foreach my $sample (@sortednames) {
	if (! exists $stotal{$sample}) {
	    $samplecheck++;
	} else {
	    next;
	}
    }
    if ($samplecheck > 0) {
	$linenum++;
	next;
    }
    
    #Calculate pi (nucleotide diversity) for each sample compared to reference allele
    my %nucdiv;
    my $check;
    foreach my $basecall (@basecalls) {
	if ($basecall eq $refbase) {
	    foreach my $sample (@sortednames) {
		my $p;
		if (exists $samplecount{$basecall}{$sample}) {
		    $p = $samplecount{$basecall}{$sample}/$stotal{$sample};
		} else {
		    $p = 0;
		}
		my $q = 1-$p;
		$nucdiv{$sample} = 2*$p*$q;
		$markerND{$marker}{$sample} +=$nucdiv{$sample};
		$markerNDlength{$marker}{$sample}++;
		$finalmarkers{$marker}++;
	    } 
	} else {
	    $check++;
	    next;
	}
    }
    #Check if no basecalls matched reference call
    if (defined $check && $check == scalar(@basecalls)) {
	foreach my $sample (@sortednames) {
	    $nucdiv{$sample} = 0;
	    $markerND{$marker}{$sample} +=$nucdiv{$sample};
	    $markerNDlength{$marker}{$sample}++;
	    $finalmarkers{$marker}++;
	}
    }
    
    #Output Pi values for each sample at each loci/site
    #Get taxon for each loci
    my $taxonpath = $taxamarker{$marker};
    my $cladelevel = $taxaLevel{$taxonpath};
    my $sampleout;

    foreach my $sample (@sortednames) {
	$sampleout .= "\t$sample\t$nucdiv{$sample}\t$stotal{$sample}";
    }

    $lineoutput = "$marker\t$pos\t$taxonpath\t$cladelevel$sampleout";
    print OUT "$lineoutput\n";
    $linenum++;
}
close(OUT);
close(MPILE);
