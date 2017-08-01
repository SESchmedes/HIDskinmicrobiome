#!/usr/bin/env perl
#########################################
#########################################
#Authored by Sarah E. Schmedes
#Date released: 08/01/2017

#weka_compare.pl will run samtools mpileup on .sam files,
# parse through the mpileup outputs, calculate the ND at 
# each shared marker, create feature vectors, .arff files,
# and run Weka classifiers logistic regression and 1NN with
# and without attribute selection
#########################################
#########################################

use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

my $metafile = $OH_METAFILE;
my ($bodysite, $threshold, $pilepath, $samplelist, $metaphlanoutpath, $wekaoutpath, $wekaexepath);
#Set flags
GetOptions('bodysite=s' => \$bodysite, 
	   'threshold=i' => \$threshold,
	   'pilepath=s' => \$pilepath,
	   'samplist=s' => \$samplelist,
	   'metaphlanoutpath=s' => \$metaphlanoutpath,
	   'wekaoutpath=s' => \$wekaoutpath,
	   'wekaexepath=s' => \$wekaexepath)
    or die "Error with flags.\n";

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
my ($taxa_ref, $taxalevel_ref, $markerlength_ref) = gettaxa();
my %taxamarker = %$taxa_ref;
my %taxalevel = %$taxalevel_ref;
my %markerlength = %$markerlength_ref;
my @totalDBmarkers = sort {$a cmp $b} keys %taxamarker;
#Run pileparseMarkerND.pl to produce new mpileup, ND calculations, and master output
if (! -f "/$pilepath/$bodysite\_all_ND.txt") {
    if (system("perl pileparseMarkerND.pl -bodysite $bodysite")) {
	die "problem executing pileparseMarkerND.pl. $!\n";
    }
}
my $nsamples = scalar(@sortednames);

open IN, "<$pilepath/$bodysite\_all_ND.txt" or die "Could not open ND.txt file for reading\n";

my %markercov;
my %samplemarkercov;
my %sampleND;
while (<IN>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^Marker/) {
	next;
    }
    my ($marker, $pos, $taxonpath, $taxon, $level, @sampledata) = split /\t/, $line;
    my $covcheck = 0;
    for (my $i=2; $i<=$#sampledata; $i+=3) {
	if ($sampledata[$i] < $threshold) {
	    $covcheck++;
	}
    }
    if ($covcheck > 0) {
	next;
    }
    $markercov{$marker}++;
    foreach my $sample (@sortednames) {
	my ($sampleID, $sampleND, $samplecov) = splice (@sampledata, 0, 3);
	if ($sample eq $sampleID) {
	    $samplemarkercov{$marker}{$sample} += $samplecov;
	    $sampleND{$marker}{$sample} += $sampleND;
	} else {
	    die "Sample do not align\n";
	}
    }
}
close(IN);

my @markers = sort {$a cmp $b} keys %markercov;

open FV, ">$pilepath/$bodysite\_markerND_fv_gteq$threshold\.txt" or die "Cannot open fv file for writing\n";
my @fvmarkers;
foreach my $marker (@markers) {
    my $samplecheck;
    foreach my $sample (@sortednames) {
	my $avgcov = $samplemarkercov{$marker}{$sample}/$markercov{$marker};
	if ($avgcov < $threshold) {
	    $samplecheck++;
	} else {
	    next;
	}
    }
    if (defined $samplecheck && $samplecheck > 0) {
	next;
    } else {
	push(@fvmarkers, $marker);
    }
}
my $markerstring = join("\t", @fvmarkers);
print FV "$markerstring\tIndividual\n";

foreach my $sample (@sortednames) {
    my $out;
    foreach my $marker (@fvmarkers) {
	my $ND = $sampleND{$marker}{$sample}/$markercov{$marker};
	$out .= "$ND\t";
    }
    my ($ind, $body, $time) = split /_/, $sample;
    $out .= $ind;
    print FV "$out\n";
}
close(FV);


if (! -f "$pilepath/$bodysite\_markerND_fv_gteq$threshold\.arff") {
#Run R script to create .arff
    if(system("Rscript /home/sarah/src/markerND_fv_to_arff.R --bodysite $bodysite --threshold $threshold --inputpath $pilepath")) {
	die "Error executing rscript to create .arff file\n";
    }
}

#    die "check .arff file creation\n"; #And double check to make sure the .arff file is saved as above
#weka commands here
if (! -f "$pilepath/$bodysite\_weka_logistic_markerND_gteq$threshold\.txt") {
#Logistic
    if(system("java -Xmx8g -classpath $wekaexepath/weka.jar weka.classifiers.functions.Logistic -R 1.0E-8 -M -1 -num-decimal-places 4 -x $nsamples -t $pilepath/$bodysite\_markerND_fv_gteq$threshold\.arff > $pilepath/$bodysite\_weka_logistic_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka logistic. $!\n";
    }
}
if (! -f "$pilepath/$bodysite\_weka_knn_markerND_gteq$threshold\.txt") {
#1NN
    if(system("java -Xmx8g -classpath $wekaexepath/weka.jar weka.classifiers.lazy.IBk -K 1 -x $nsamples -t $pilepath/$bodysite\_markerND_fv_shared_gteq$threshold\.arff > $pilepath/$bodysite\_weka_knn_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka knn. $!\n";
    }
}
if (! -f "$pilepath/$bodysite\_weka_attselectLogistic_markerND_gteq$threshold\.txt") {
#AttSelectLog
    if(system("java -Xmx8g -classpath $wekaexepath/weka.jar weka.classifiers.meta.AttributeSelectedClassifier -x $nsamples -t $pilepath/$bodysite\_markerND_fv_shared_gteq$threshold\.arff -W weka.classifiers.functions.Logistic -- -R 1.0E-8 -M -1 -num-decimal-places 4 > $pilepath/$bodysite\_weka_attselectLogistic_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka attselect logistic. $!\n";
    }
}
if (! -f "$pilepath/$bodysite\_weka_attselectknn_markerND_gteq$threshold\.txt") {
#AttSelect1NN
    if(system("java -Xmx8g -classpath $wekaexepath/weka.jar weka.classifiers.meta.AttributeSelectedClassifier -x $nsamples -t $pilepath/$bodysite\_markerND_fv_shared_gteq$threshold\.arff -W weka.classifiers.lazy.IBk -- -K 1 > $pilepath/$bodysite\_weka_attselectknn_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka attselect knn euclidean. $!\n";
    }
}
my @classifiers =("logistic", "knn", "attselectLogistic", "attselectknn");
my (%attributes, %correctclass, %accuracy, %incorrectclass);
foreach my $class (@classifiers) {
    open CLASS, "</$pilepath/$bodysite\_weka_$class\_markerND_gteq$threshold\.txt" or die "Cannot open weka results for reading. $!\n";
    if ($class eq "logistic" || $class eq "knn" || $class eq "multiclasslog") {
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
#Get all taxon paths for each feature vector marker
#and number of markers per clade
my %paths;
foreach my $marker (@fvmarkers) {
    my $taxonpath = $taxamarker{$marker};
    $paths{$taxonpath}++;
}
my @paths = sort {$a cmp $b} keys %paths;
#Get clade abundances per sample
my %cladeabund;
open CLADE, "<$metaphlanoutpath/$bodysite\_metaphlanlongout.txt" or die "Cannot open metaphlanlong out file for reading. $!\n";
while (<CLADE>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^TaxonPath/) {
	next;
    }
    my ($taxonpath, $taxon, $level, @sampleabund) = split /\t/, $line;
    foreach my $path (@paths) {
	if ($taxonpath eq $path) {
	    for (my $i=0; $i<=$#sampleabund; $i++) {
		$cladeabund{$path}{$sortednames[$i]} = $sampleabund[$i];
	    }
	}
    }
}
close(CLADE);

foreach my $path (@paths) {
    foreach my $sample (@sortednames) {
	if (! exists $cladeabund{$path}{$sample}) {
	    $cladeabund{$path}{$sample} = 0;
	}
    }
}
#Get number of markers in database per wanted clade
my %clademarkerDBcount;
foreach my $path (@paths) {
    foreach my $marker (@totalDBmarkers) {
	if ($path eq $taxamarker{$marker}) {
	    $clademarkerDBcount{$path}++;
	} else {
	    next;
	}
    }
}
#Get number of FS markers per clade
my %featurecheck;
foreach my $marker (@fvmarkers) {
    my $revise = $marker;
    $revise =~ s/[|:-]/./g;
    foreach my $class (@classifiers) {
	if (exists $attributes{$revise}{$class}) {
	    $featurecheck{$marker}++;
	} else {
	    next;
	}
    }
}
my @FSmarkers = sort {$a cmp $b} keys %featurecheck;
foreach my $marker (keys %featurecheck) {
    if ($featurecheck{$marker} != 3) {
	die "There are differences in the marker attributes selected for each AttSelect classifier\n";
    }
}
my %fsmarkcount;
foreach my $fsmark (@FSmarkers) {
    my $path = $taxamarker{$fsmark};
    $fsmarkcount{$path}++;
}
foreach my $path (@paths) {
    if (! exists $fsmarkcount{$path}) {
	$fsmarkcount{$path} = 0;
    }
}

#Prepare line out
open RESULTS, ">$wekaoutpath/$bodysite\_markerND_weka_compare_$threshold\_MASTER.txt" or die "Cannot open output for writing. $!\n";
print RESULTS "Marker\tMarkerLengthDB\tMarkerCoverage\tTaxonPath\tTaxon\tLevel\tTaxonMarkerCountDB\tTaxonAbund\tSampleID\tSampleAvgMarkerCov\tBodysite\tTotalCladeMarkersAligned\tSelectedFeatures\tLogistic\tKNN\tAttSelectLog\tAttSelectKNN\tAttSelectRF500\tMultiClassLogistic\n"; 
open SHORT, ">$wekaoutpath/$bodysite\_markerND_weka_compare_$threshold\_short.txt" or die "Cannot open master output for writing. $!\n";
print SHORT "Bodysite\tTaxonPath\tTaxon\tLevel\tNumofMarkersDB\tNumofMarkersAligned\tNumofMarkersFSClade\tLogistic\tKNN\tAttSelectLog\tAttSelectKNN\tAttSelectRF500\tMultiClassLog\n";
my $accuracyout;
foreach my $class (@classifiers) {
    $accuracyout .= "\t$accuracy{$class}";
}
foreach my $marker (@fvmarkers) {
    my $FS;
    my $length = $markerlength{$marker};
    my $cov = $markercov{$marker};
    my $taxonpath = $taxamarker{$marker};
    my $cladelevel = $taxalevel{$taxonpath};
    my ($clade, $level) = split /\t/, $cladelevel;
    if (exists $featurecheck{$marker}) {
	$FS = 1;
    } else {
	$FS = 0;
    }
    foreach my $sample (@sortednames) {
	my $abund = $cladeabund{$taxonpath}{$sample};
	my $avgsitecov = $samplemarkercov{$marker}{$sample}/$cov;
	my $clademarkers = $paths{$taxonpath};
	my $lineout = "$marker\t$length\t$cov\t$taxonpath\t$clade\t$level\t$clademarkerDBcount{$taxonpath}\t$abund\t$sample\t$avgsitecov\t$bodysite\t$clademarkers\t$FS$accuracyout";
	print RESULTS "$lineout\n";
    }
}
#Condensed output
foreach my $path (@paths) {
    my $cladelevel = $taxalevel{$path};
    my ($clade, $level) = split /\t/, $cladelevel;
    my $markersforcladeDB = $clademarkerDBcount{$path};
    my $markersforclade = $paths{$path};
    my $fsmarkerforclade = $fsmarkcount{$path};
    my $out = "$bodysite\t$path\t$clade\t$level\t$markersforcladeDB\t$markersforclade\t$fsmarkerforclade$accuracyout";
    print SHORT "$out\n";
}
close(RESULTS);
close(SHORT);
