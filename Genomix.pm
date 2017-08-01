#!/usr/bin/env perl

#Authored by August Woerner and Sarah Schmedes
#Date released: 08/01/2017


package Genomix;

use Exporter;         # Gain export capabilities

our (@ISA, @EXPORT,   # Declare some global variables
    @EXPORT_OK, %EXPORT_TAGS,
     # constants you want to define go below here:
     $BWA, $SAMTOOLS, $OH_METAFILE, $TAXA_MARKER, $MARKER_LENGTH);


# Define the absolute path of the binaries
#INPUT YOUR PATHS HERE####
$BWA = '/usr/bin/bwa';
$SAMTOOLS = '/usr/local/bin/samtools-1.3.1';
$OH_METAFILE = "/home/sarah/OhDatasets/OhMeta.txt";
$TAXA_MARKER = "/home/sarah/strainphlan/db_v20/MarkerPerTaxon.txt";
$MARKER_LENGTH = "/home/sarah/strainphlan/db_v20/MarkerLengths.txt";


@ISA =qw(Exporter);   # Take advantage of Exporter's capabilities

# constants you can ask for by name
@EXPORT_OK =          # Exported, if explicitly requested
    qw($BWA $SAMTOOLS $OH_METAFILE $TAXA_MARKER $MARKER_LENGTH queryDB getDBHeaders gettaxa samplename);

# or by group (preferred)
%EXPORT_TAGS =        # Tag groups
     ("constants" =>[qw($BWA $SAMTOOLS $OH_METAFILE $TAXA_MARKER $MARKER_LENGTH)],
     "meta" => [qw(queryDB getDBHeaders gettaxa samplename)]);


# script globals
my @dbHeader= (); # the column labels
my @file = (); # the file (sans the header)
my $dbFilename=''; # the "database"
my $sep="\t"; # the column separator


# Subroutine definitions go here 
sub getDBHeader {
    return @dbHeader;
}

sub loadDB {
    my $file = shift;

# if it's a new file, read it in to memory
    
    open IN, $file or
	die "Failed to open $file for reading!\n";
    
    $_ = <IN>;
    chomp;
    @dbHeader = split /$sep/;
    my %uniq = ();
    foreach (@dbHeader) {
	$_ = uc $_; # for case insensitive comparisons
	$uniq{$_}++;
	if ($uniq{$_} > 1) {
		die "Error with file: $file\nThe column $_ appears at least twice!\nThat cannot be!\n";
	}
    }
    
    @file = ();
    my $errs=0;
    while (<IN>) {
	chomp;
	my @s = split /$sep/;
	if (@s != @dbHeader) { # wrong number of columns...
	    $errs++;
	}
	push(@file, \@s);
    }
    
    close IN;
    if ($errs) {
	die "Your file format is unexpected.\nBased on the header, I was expecting " , scalar(@dbHeader) , " columns to be present, but there were " , $errs , " rows where that was not true!\n";
    }
    $dbFilename = $file;

}

#queryDB allows the user to obtain any metadata fields in an array from the OhMeta.txt doc
#If 2 arguments are used, the first argument is the $metafile and the 
#second argument is the field you want to retreive Ex: @myquery = queryDB($metafile, "srsID");
#If 4 arguments are used, such as @myquery = queryDB($metafile, "Mb", "bodysiteID", "srsID"), will
#result in an array with all srsIDs for the bodysite, Mb.
sub queryDB {
    my ($file, $val, $searchCol, $resultCol) = @_;
# if it's a new file, read it in to memory
    if ($file ne $dbFilename) {
	loadDB($file);
    }

# from here on, the file has been loaded into @file

# used to get all of the uniq column values
    if (@_ == 2) {
	my $outIndex=-1;
	for (0..$#dbHeader) { # figure out which column we want (by name, get the index)
	    if ($dbHeader[$_] eq uc($val)) {
		$outIndex = $_;
	    }
	}
	
	if ($outIndex == -1) {
	    die "Oops! I tried to find column $resultCol, but I couldn't find it in:\n@dbHeader!\n";
	}
	
	my %out = ();
	foreach (@file) {
	    if (defined $_->[ $outIndex]) {
		$out{$_->[$outIndex]}=1; # used to make the array unique
	    }
	}
	return sort keys %out;
    } else {

# we search $searchCol for $val, and return the unique values from the $resultCol
	my $inIndex=-1;
	my $outIndex=-1;
	for (0..$#dbHeader) {
	    if ($dbHeader[$_] eq uc($searchCol)) {
		$inIndex = $_;
	    }
	    if ($dbHeader[$_] eq uc($resultCol)) {
		$outIndex = $_;
	    }
	}
	
	my $errs=0;
	if ($inIndex == -1) {
	    ++$errs;
	    warn "Oops! I tried to find column $searchCol, but I couldn't find it in:\n@dbHeader!\n";
	}
	if ($outIndex == -1) {
	    ++$errs;
	    warn "Oops! I tried to find column $resultCol, but I couldn't find it in:\n@dbHeader!\n";
	}
	die if $errs;

	my %out = ();
	foreach (@file) {
	    if (defined $_->[ $inIndex] && uc($_->[ $inIndex]) eq uc($val) ) {
		$out{$_->[$outIndex]}=1;
	    }
	}
	return sort keys %out;
    }
}

#samplename function result in a hash of samplenames. Keys are srsIDs and samplenames are in the 
#format subjectID_bodysiteID_timepoint. Argument to samplename() is an array of srsIDs
#If only one SRSid is provided as an argument a scalar will be returned for that particular samplename
sub samplename {
    my $metafile = $OH_METAFILE;
    my @SRSs = @_;
    if (scalar(@SRSs) > 1) {
	my %metadata;
	foreach my $SRS (@SRSs) {
	    my @subjectID = queryDB($metafile, $SRS, "srsID", "subjectID");
	    my @bodysiteID = queryDB($metafile, $SRS, "srsID", "bodysiteID");
	    my @timepoint = queryDB($metafile, $SRS, "srsID", "timepoint");
	    $metadata{$SRS} = "$subjectID[0]\_$bodysiteID[0]\_$timepoint[0]";
	}
	return (\%metadata);
    }
    if (scalar(@SRSs) == 1) {
	my $SRS = $SRSs[0];
	my @subjectID = queryDB($metafile, $SRS, "srsID", "subjectID");
	my @bodysiteID = queryDB($metafile, $SRS, "srsID", "bodysiteID");
	my @timepoint = queryDB($metafile, $SRS, "srsID", "timepoint");
	my $samplename = "$subjectID[0]\_$bodysiteID[0]\_$timepoint[0]";
	return ($samplename);
    }
}

#gettaxa() will return three hashes; no arguments required. This pulls data from the MetaPhlAn2 database
#%taxa{marker} = full taxonpath
#%taxaLevel{full taxonpath} = terminalclade \t taxonomic level
#%markerlength{marker} = markerlength in metaphlan2 database
sub gettaxa {
    my (%taxa, %taxacheck, %taxaLevel, %markerlength);
    my $fullpath;
    my $Taxon;
    my @tlevels = ("", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies");
    open TAXA, "<$TAXA_MARKER" or die "Could not open $TAXA_MARKER file for reading";
    while (<TAXA>) {
	my $line = $_;
	chomp $line;
	my ($marker, $taxonpath) = split /\s+/, $line;
	$taxa{$marker} = $taxonpath;
	$taxacheck{$marker}++;
    }
    foreach my $key (keys %taxacheck) {
	if ($taxacheck{$key} >1) {
	    die "$key is not specific to a taxonpath! Check taxa marker db!\n";
	} else {
	    next;
	}
    }
    my @allmarkers = sort keys %taxa;
    foreach my $markerID (@allmarkers) {
	$fullpath = $taxa{$markerID};
	my $taxLevel = $fullpath =~ s/\|/|/g;
	if ($fullpath =~ m/\|?\w_+([^\|]+)$/) {
	    $Taxon = $1;
	}
	else {
	    die "Regex broke!\n";
	}
	++$taxLevel;
	if ($taxLevel > $#tlevels) {
	    die "Unknown taxonomic level with: $fullpath $taxLevel\n";
	}
	my $t = $tlevels[$taxLevel];
	$Taxon .= "\t".$t;
	$taxaLevel{$fullpath} = $Taxon;
    }
    close(TAXA);
    open LENGTH, "<$MARKER_LENGTH" or die "Cannot open $MARKER_LENGTH for reading\n";
    while (<LENGTH>) {
	my $line = $_;
	chomp $line;
	my ($marker, $length, $score) = split /\t/, $line;
	$markerlength{$marker} = $length;
    }
    close(LENGTH);
    return (\%taxa, \%taxaLevel, \%markerlength);
}







# 1 is necessary for perl packages
1;
