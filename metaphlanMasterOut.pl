#!/usr/bin/env perl

#################################################
#Created by Sarah Schmedes
#This script reads in metaphlan_longout.txt file and _marker_pres_table.txt
#from metaphlan2 output
#################################################
#Date released: 08/01/2017
use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

my $metafile = $OH_METAFILE;
my ($body, $markerpath, $longpath, $featureout, $masterout, $markerfv);
my @SRS;
my $SRSlist;
my (%samplemarkers, %snamemarkers, %markers, %clade, %allclades);

#Set flags
GetOptions('bodysite=s' => \$body,
	   'markerpath=s' => \$markerpath,
	   'longpath=s' => \$longpath,
	   'featureout=s' => \$featureout,
	   'output=s' => \$masterout,
	   'markerfv=s' => \$markerfv)
    or die "Error with flag. Designate -bodysite with bodysite symbol\n";

if ($body) {
    @SRS = queryDB($metafile, $body, "bodysiteID", "srsID");
}

#Get taxon levels for each markerID
my ($taxa_ref, $taxaLevel_ref, $markerlength_ref) = gettaxa();
my %taxamarker = %$taxa_ref; #hash of each marker and associated full taxon path
my %taxaLevel = %$taxaLevel_ref; #hash of each full taxon path and associated 
                                 #leaf clade tab separated corresponding taxon level
my %markerlength = %$markerlength_ref;

my $samplename_ref = samplename(@SRS);
my %samplename = %$samplename_ref;
my @sortmeta = sort{$a cmp $b} values %samplename;
my @SRSlist;
#Get an ordered array by samplename
foreach my $samplename (@sortmeta) {
    foreach my $key (keys %samplename) {
	if ($samplename{$key} eq $samplename) {
	    push(@SRSlist, $key);
	}
    }
}

foreach my $SRS (@SRSlist) {
    open MARKER, "<$markerpath/$SRS\_marker_pres_table.txt\n" or die "Cannot open marker file for reading\n";
    my @samplemarkers;
    my $samplename = $samplename{$SRS};
    while (<MARKER>) {
	my $line = $_;
	chomp $line;
	if ($line =~ /^#SampleID/) { #This is the header line
	    next;
	}
	my @content = split /\s+/, $line;
	my $marker = $content[0]; #hash with all markers across all samples tested
	$markers{$marker} = 1;
	push(@samplemarkers, $marker);
    }
    $samplemarkers{$SRS} = \@samplemarkers; #hash with array of all marker hits for each sample
    $snamemarkers{$samplename} = \@samplemarkers; #hash with array of marker hits for each samplename
    close(MARKER);
}

open LONG, "<$longpath/metaphlanlongout_$body\.txt" or die "Cannot open metaphlanlongout.txt file for reading\n";
my @samples;
while (<LONG>) {
    my ($head1, $head2, $head3);
    my $line = $_;
    chomp $line;
    if ($line =~ /^TaxonPath/) { #This is the header line
	($head1, $head2, $head3, @samples) = split /\s+/, $line;
	next;
    }
    my ($taxonpath, $clade, $level, @abundance) = split /\s+/, $line;
    if (@sortmeta != @samples) {
	die "Samples chosen from flag differ from samples in metaphlanlongout.txt file\n";
    }
    if (scalar(@samples) != scalar(@abundance)) {
	die "There is not an abundance value for each sample\n";
    }
    my $nsamples = scalar(@samples);
    for (my $i=0; $i<=$#samples; $i++) {
	    $clade{$clade}{$samples[$i]} = $abundance[$i];
	    $allclades{$clade}=1;
    }
}
close(LONG);

#Generate feature vector for clades
open CLADEFV, ">$featureout/featurevector_clade_$body\.txt" or die "Cannot open file for writing\n";
my @clades = sort {$a cmp $b} keys %allclades;
my $clade_s = join("\t", @clades);
print CLADEFV "$clade_s\tIndividual\n";
foreach my $sample (@samples) {
    my $out;
    foreach my $clade (@clades) {
	my $abund = $clade{$clade}{$sample};
	$out .= "$abund\t";
    }
    my @samplename = split /_/, $sample;
    my $ind = $samplename[0];
    $out .= $ind;
    print CLADEFV "$out\n";
}
close(CLADEFV);


#Generate masteroutput
open OUT, ">$masterout/metaphlanMasterOut_$body\.txt" or die "Cannot open masterout file for writing\n";

print OUT "Marker\tTaxonPath\tLevel\tClade\tCladeRelAbund\tsrsID\tSubjectID\tBodysiteID\tTimepoint\n"; #header line

my @totalmarkers = sort {$a cmp $b} keys %markers;
my %markerfv;
foreach my $marker (@totalmarkers) {
    foreach my $sample (@samples) {
	my $markers_ref = $snamemarkers{$sample};
	my @samplemarkers = @$markers_ref; #Retrieve all markers for sample
	foreach my $samplemarker (@samplemarkers) {
	    my $srsID;
	    if ($samplemarker eq $marker) {
		my $taxonpath = $taxamarker{$marker};
		my $cladelevel = $taxaLevel{$taxonpath};
		my @cladelevel = split /\t/, $cladelevel;
		my $clade = $cladelevel[0];
		if (exists $clade{$clade}{$sample}) { #check if clade is counted as present
		    $markerfv{$marker}{$sample} = 1;
		    my $level = $cladelevel[1];
		    my $abund = $clade{$clade}{$sample};
		    my ($subjectID, $bodysiteID, $timepoint) = split /_/, $sample;
		    foreach my $SRS (@SRS) {
			if ($sample eq $samplename{$SRS}) {
			    $srsID = $SRS;
			} else {
			    next;
			}
		    }
		    print OUT "$marker\t$taxonpath\t$level\t$clade\t$abund\t$srsID\t$subjectID\t$bodysiteID\t$timepoint\n";
		} else {
		    $markerfv{$marker}{$sample} = 0;
		}
	    } else {
		next;
	    }
	}
    }
}
close(OUT);

#Generate feature vector for all markers indicating presence/absence for each sample
open MARKERFV, ">$markerfv/$body\_allmarkers_fv.txt" or die "Cannot open marker feature vector for writing\n";

my $markerhead = join("\t", @totalmarkers);
print MARKERFV "$markerhead\tIndividual\n";

foreach my $sample (@samples) {
    my $fvout;
    foreach my $marker (@totalmarkers) {
	if (exists $markerfv{$marker}{$sample}) {
	    $fvout .= "$markerfv{$marker}{$sample}\t"; 
	} else {
	    $fvout .= "0\t";
	}
    }
    my ($subjectID, $bodysiteID, $timepoint) = split /_/, $sample;
    $fvout .= "$subjectID";
    print MARKERFV "$fvout\n";
}
close(MARKERFV);
