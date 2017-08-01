#!/usr/bin/env perl
#############################################
#Created by Sarah Schmedes
#This script executes the entire panphlan workflow
# and classification using 1NN and RMLR in Weka
#Date Released: 08/01/2017
#############################################
use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

my $metafile = $OH_METAFILE;
my ($bodysite, $clade, $cladedb, $species, $samplelist, $pilepath, $fastqpath);
#Set flags
GetOptions('bodysite=s' => \$bodysite,
	   'clade=s' => \$clade,
	   'cladedb=s' => \$cladedb,
	   'species=s' => \$species,
	   'output=s' => \$pilepath,
	   'samplelist=s' => \$samplelist,
	   'fastqpath=s' => \$fastqpath)
    or die "Error with flags. Must designate -bodysite, -clade, -species, -output, and -samplelist\n";

#Get SRSids for bodysite
my @SRS = queryDB($metafile, $bodysite, "bodysiteID", "srsID");
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
my $samplename_ref = samplename(@SRS);
my %allsamplenames = %$samplename_ref;
my @allSRS = sort {$a cmp $b} keys %allsamplenames;
my @orderedSRS;
#Get ordered array by samplenames
foreach my $samplename (@sortednames) {
    foreach my $srs (@allSRS) {
	if ($samplename eq $allsamplenames{$srs}) {
	    push(@orderedSRS, $srs);
	}
    }
}
#Concatenate .bz2 files containing SRRs for each SRS and pipe into panphlan    
#Run copy bodysite panphlan outputs to new subdirectories
foreach my $SRS (@orderedSRS) {
    #Collect all SRR numbers for each SRS sample
    my @sampleSRRs = queryDB($metafile, $SRS, "srsID", "srrID");
    my $SRRstring;
    foreach my $SRR (@sampleSRRs) {
	$SRRstring .= "$fastqpath/$SRR.1.qc.fastq.bz2 $fastqpath/$SRR.2.qc.fastq.bz2 ";
    }
    if (! -f "$pilepath/$bodysite/$allsamplenames{$SRS}\_$species\.csv.bz2") {
	if (system("bzcat $SRRstring | python panphlan_map.py -c $species --fastx fastq --i_bowtie2_indexes $cladedb -p 20 -m 10 -o $pilepath/$bodysite/$allsamplenames{$SRS}\_$species --verbose")) {
	    die "Error while running panphlan_map.py: $!\n";
	}
    }
}
print "panphlan_map.py complete!\n";
#Run panphlan profile and hclust2.py
if (system("python panphlan_profile.py -c $species --i_bowtie2_indexes $cladedb -i $pilepath/$bodysite/ --o_dna $pilepath/$bodysite/$bodysite\_$species\_gene_presence_absence_noref --verbose")) {
    die "Error while running panphlan_profile.py: $!\n";
}
print "panphlan_profile.py complete for $bodysite\!\n";
if (system("python hclust2.py -i $pilepath/$bodysite/$bodysite\_$species\_gene_presence_absence_noref -o $pilepath/$bodysite/$bodysite\_$species\_heatmap.png --legend_file /$pilepath/$bodysite/$bodysite\_$species\_legend.png --image_size 30 --cell_aspect_ratio 0.01 --dpi 300 --slabel_size 24 --f_dist_f jaccard --s_dist_f jaccard")) {
    die "Error while running hclust2.py for $bodysite\n";
}
print "hclust2.py complete for $bodysite\n";

if (! -f "$pilepath/$bodysite/$bodysite\_$species\_pangenePA_fv.arff") {
    if (system(" Rscript /home/sarah/src/panphlanOUT_to_fv_arff.R -bodysite $bodysite -clade $species -output $pilepath")) {
	die "Error while trying to run R file to make arff\n";
    }
}
#Get number of samples
my $lines = `wc -l $pilepath/$bodysite/$bodysite\_$species\_pangenePA_fv.txt`;
my @value = split /\s+/, $lines;
my $nsamples = $value[0] - 1; 
#    die "check .arff file creation\n"; #And double check to make sure the .arff file is saved as above
#weka commands here
if (! -f "$pilepath/$bodysite/$bodysite\_weka_logistic_$species\_PA.txt") {
#Logistic
    if(system("java -Xmx8g -classpath weka.jar weka.classifiers.functions.Logistic -R 1.0E-8 -M -1 -num-decimal-places 4 -x $nsamples -t $pilepath/$bodysite/$bodysite\_$species\_pangenePA_fv.arff > $pilepath/$bodysite/$bodysite\_weka_logistic_$species\_PA.txt")) {
	die "Error executing java and weka logistic. $!\n";
    }
}
if (! -f "$pilepath/$bodysite/$bodysite\_weka_knn_$species\_PA.txt") {
#KNN
    if(system("java -Xmx8g -classpath weka.jar weka.classifiers.lazy.IBk -K 1 -x $nsamples -t $pilepath/$bodysite/$bodysite\_$species\_pangenePA_fv.arff > $pilepath/$bodysite/$bodysite\_weka_knn_$species\_PA.txt")) {
	die "Error executing java and weka knn. $!\n";
    }
}
if (! -f "$pilepath/$bodysite/$bodysite\_weka_attselectLogistic_$species\_PA.txt") {
#AttSelectLog
    if(system("java -Xmx8g -classpath weka.jar weka.classifiers.meta.AttributeSelectedClassifier -x $nsamples -t $pilepath/$bodysite/$bodysite\_$species\_pangenePA_fv.arff -W weka.classifiers.functions.Logistic -- -R 1.0E-8 -M -1 -num-decimal-places 4 > $pilepath/$bodysite/$bodysite\_weka_attselectLogistic_$species\_PA.txt")) {
	die "Error executing java and weka attselect logistic. $!\n";
    }
}
if (! -f "$pilepath/$bodysite/$bodysite\_weka_attselectknn_$species\_PA.txt") {
#AttSelectKNN
    if(system("java -Xmx8g -classpath weka.jar weka.classifiers.meta.AttributeSelectedClassifier -x $nsamples -t $pilepath/$bodysite/$bodysite\_$species\_pangenePA_fv.arff -W weka.classifiers.lazy.IBk -- -K 1 > $pilepath/$bodysite/$bodysite\_weka_attselectknn_$species\_PA.txt")) {
	die "Error executing java and weka attselect knn euclidean. $!\n";
    }
}

my @classifiers =("logistic", "knn", "attselectLogistic", "attselectknn");
my (%attributes, %correctclass, %accuracy, %incorrectclass, %fvmarkers);
foreach my $class (@classifiers) {
    open CLASS, "</$pilepath/$bodysite/$bodysite\_weka_$class\_$species\_PA.txt" or die "Cannot open weka results for reading. $!\n";
    if ($class eq "logistic" || $class eq "knn") {
	while (<CLASS>) {
	    my $line = $_;
	    chomp $line;
	    if ($line =~ /=== Stratified cross-validation ===/) {
		my $blank = <CLASS>;
		chomp $blank;
		my $correctline = <CLASS>;
		chomp $correctline;
		my @correctinfo = split /\s+/, $correctline;
		$correctclass{$class} = $correctinfo[3];
		$accuracy{$class} = $correctinfo[4];
		my $incorrectline = <CLASS>;
		chomp $incorrectline;
		my @incorrectinfo = split /\s+/, $incorrectline;
		$incorrectclass{$class} = $incorrectinfo[3];
	    }
	}
    } else {
	while (<CLASS>) {
	    my $line = $_;
	    chomp $line;
	    if ($line =~ /Selected attributes:/) {
		my @att = split /\s+/, $line;
		my $natt = $att[$#att];
		for my $n (1..$natt) {
		    my $m = <CLASS>;
		    chomp $m;
		    my @marker = split /\s+/, $m;
		    $attributes{$marker[1]}{$class}++;
		    $fvmarkers{$marker[1]}++;
		}
	    }
	    if ($line =~ /=== Stratified cross-validation ===/) {
		my $blank = <CLASS>;
		chomp $blank;
		my $correctline = <CLASS>;
		chomp $correctline;
		my @correctinfo = split /\s+/, $correctline;
		$correctclass{$class} = $correctinfo[3];
		$accuracy{$class} = $correctinfo[4];
		my $incorrectline = <CLASS>;
		chomp $incorrectline;
		my @incorrectinfo = split /\s+/, $incorrectline;
		$incorrectclass{$class} = $incorrectinfo[3];
	    }
	}
    }	
    close(CLASS);
}
my @fvmarkers = keys %fvmarkers;
my $nfvmarkers = scalar(@fvmarkers);
#Prepare line out
open SHORT, ">$pilepath/$bodysite/$bodysite\_$species\_PA_weka_compare_short.txt" or die "Cannot open master output for writing. $!\n";
print SHORT "Bodysite\tNumofMarkersFS\tLogistic\tKNN\tAttSelectLog\tAttSelectKNN\tAttSelectRF500\tMultiClassLog\n";
my $accuracyout;
foreach my $class (@classifiers) {
    if (! exists $accuracy{$class}) {
	die "$class\n";
    }
    $accuracyout .= "\t$accuracy{$class}";
}
my $out = "$bodysite\t$nfvmarkers\t$accuracyout";
print SHORT "$out\n";
close(SHORT);
