#!/usr/bin/env perl
##############################################################
#Created by Sarah Schmedes and August Woerner
# this reads some metaphlan _profile.txt files
# outputs in long format and reports taxa that are common to all samples
# files in the ARGV or with bodysite flag
##############################################################
#Date released: 08/01/2017
use strict;
use warnings;
use Getopt::Long;
use Genomix qw(:constants :meta);

my $metafile = $OH_METAFILE;
my ($bodysite, $samplelist, $metaphlanoutpath, $outpath);
my @SRS;
my @SRSlist;
my $filename;
my %bodysamplename;
my @sortedbodynames;
#Set flags
GetOptions('bodysite=s' => \$bodysite,
	   'samplelist=s' => \$samplelist,
	   'metaphlanoutpath=s' => \$metaphlanoutpath,
	   'outpath=s' => \$outpath)
    or die "Error with flag. Designate -bodysite with bodysite symbol\n";

@SRS = queryDB($metafile, $bodysite, "bodysiteID", "srsID");
my $samplename_ref = samplename(@SRS);
%bodysamplename = %$samplename_ref;
@sortedbodynames = sort {$a cmp $b} values %bodysamplename;

open SAMPLES, "<$samplelist" or die "Could not open sample list for reading! $!\n";
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
foreach my $samplename (@sortednames) {
    foreach my $key (keys %bodysamplename) {
	if ($bodysamplename{$key} eq $samplename) {
	    push(@SRSlist, $key);
	}
    }
}

my @tlevels = ("",
	       "Kingdom",
	       "Phylum",
	       "Class", 
	       "Order", 
	       "Family",
	       "Genus", 
	       "Species",
	       "Subspecies");

my %allTogether = ();
my %files;

foreach my $SRS (@SRSlist) {
    open IN, "<$metaphlanoutpath/$SRS\_profile.txt" or die "Cannot open _profile.txt file for reading\n";
    while (<IN>) {
	chomp;
	next if /^#/;
	
	my @s = split /\s+/;
	$allTogether{$s[0]}{$SRS} = $s[1];
	$files{$SRS}=1;
    }
    close(IN);
}
my @taxa = sort keys %allTogether;
my @f = sort keys %files;

open OUT, ">$outpath/$bodysite\_metaphlanlongout.txt" or die "Cannot open file for writing\n";

print OUT "TaxonPath\tTaxon\tLevel\t" , join("\t", @sortednames) , "\n";

foreach my $taxon (@taxa) {
    my $s = $taxon;
    my $taxonomicLevel= $s =~ s/\|/|/g; # count the number of pipes (| separate the taxon labels)

# surely this will break...    
    if ($s =~ m/\|?\w_+([^\|]+)$/) {
	$s = $s . "\t" . $1;
    } else {
	die "Regex broke!\n";
    }
    ++$taxonomicLevel;
    if ($taxonomicLevel > $#tlevels) {
	die "Unknown taxonomic level with: $s $taxonomicLevel\n";
    }
    my $t = $tlevels[$taxonomicLevel];
    $s .= "\t" . $t;
    foreach (@f) {
	if (exists $allTogether{$taxon}{$_}) {
	    $s .= "\t" . $allTogether{$taxon}{$_};
	} else {
	    $s .= "\t0";
	}
    }

    print OUT $s ,  "\n";
}
close(OUT);

open LONG, "<$outpath/$bodysite\_metaphlanlongout.txt" or die "Cannot open file for reading\n";

open COMMON, ">$outpath/$bodysite\_metaphlancommonout.txt" or die "Cannot open file for writing\n";

while (<LONG>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^TaxonPath/) { #This is the header line
	print COMMON "$line\n";
	next;
    }
    my @field = split /\s+/, $line;
    my $field;
    my $t = 0;
    foreach $field (@field) {
	if ($field =~ /^\d/) {
	    if ($field != 0) {
		$t++;
	    }
	}
    }
    if ($t == ((scalar @field)-3)) { #All taxa shared in each sample
	print COMMON "$line\n";
    } else {
	next;
    }
}
close (LONG);
close (COMMON);
