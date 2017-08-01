#!/usr/bin/env perl
######################################################
#Created by Sarah Schmedes
#Will run entire strainphlan workflow, including images
#Designate -clade flag
#designated clade when prompted in STDIN
#Date Released: 08/01/2017
#####################################################
use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

#select samples to run through strainphlan
my $metafile = $OH_METAFILE;
my ($bodysite, $samplelist, $sampath, $straindir, $straindbpath);
#set flags
GetOptions('bodysite=s' => \$bodysite,
	   'samplelist=s' => \$samplelist,
	   'sampath=s' => \$sampath,
	   'output=s' => \$straindir
	   'straindb_path=s' => \$straindbpath) 
    or die "Must use bodysite, sampelist, sampath, and output flags\n"; 

#collect all SRSs for a particular body site
my @allbodySRS = queryDB($metafile, $bodysite, "bodysiteID", "srsID");

my $samplename_ref = samplename(@allbodySRS);
my %bodysamplename = %$samplename_ref;
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
my @SRSlist;
foreach my $samplename (@sortednames) {
    foreach my $key (keys %bodysamplename) {
	if ($bodysamplename{$key} eq $samplename) {
	    push(@SRSlist, $key);
	}
    }
}

#run strainphlan on all samples
#all metaphlan2/strainphlan scripts must be in the environmental path
foreach my $sample (@SRSlist) {
    #Generate a marker file for each sample.
    #The marker file contains the consensus of unique marker genes for each species found in the sample
    #This marker file will be used for SNP profiling
    if (! -f "$straindir/$bodysite/consensus_markers/$sample\.markers") {
	if (system("sample2markers.py --ifn_samples $sampath/$sample\.sam.bz2 --input_type sam --output_dir $straindir/$bodysite/consensus_markers --nprocs 20 1> $straindir/$bodysite/consensus_markers/$sample\_log.txt 2> $straindir/$bodysite/consensus_markers/$sample\_error.txt"))
	{
	    die "Strainphlan sample2markers.py ERROR: $!\n";
	}
    }
}

    #Run strainphlan to identify clades that were detected in all samples
    #providing the marker files generated in the prior step
    #to see which clades can be SNP-profiled
if (! -f "$straindir/$bodysite/clades/$bodysite\_clades.txt") {
    if (system("strainphlan.py --ifn_samples $straindir/$bodysite/consensus_markers/*.markers --output_dir $straindir/$bodysite/clades --print_clades_only --nprocs_main 20 1> $straindir/$bodysite/clades/$bodysite\_clades.txt 2> $straindir/$bodysite/clades/$bodysite\_errorlog.txt"))
    {
	die "strainphlan.py ERROR: $!\n";
    } 
}
print "Clade to analyze\n";
my $clade = <STDIN>;
chomp $clade;

print "Please specify reference genome NAME (NAME.fna.bz2)\n";
my $refgen = <STDIN>;
chomp $refgen;

if (! -f "$straindir/db_markers/$clade\.markers.fasta") {
    #Build reference database for the designated clade
    #This step only needs to be done once for each species for all projects
    if (system("extract_markers.py --mpa_pkl $straindbpath/mpa_v20_m200.pkl --ifn_markers $straindir/db_markers/all_markers.fasta --clade $clade --ofn_markers $straindir/db_markers/$clade\.markers.fasta"))
    {
	die "Strainphlan extract_markers.py ERROR: $!\n";
    }
}
else {
    print "$clade\.markers.fasta already exits\n";
}

#Build the multiple sequence alignment and phylogenetic tree
#Will align and clean sample-reconstructed strains (stored in .markers)
#and reference-genome-reconstructed strains (from clade.markers.fasta)
#Builds tree using RAxML
#If a reference genome is not specified or if no clade is specified then
#strainphlan.py will build the tree for all species it can detect
if (system("strainphlan.py --mpa_pkl $straindbpath/mpa_v20_m200.pkl --ifn_samples $straindir/$bodysite/consensus_markers/*.markers --ifn_markers $straindir/db_markers/$clade\.markers.fasta --ifn_ref_genomes $straindir/reference_genomes/$refgen\.fna.bz2 --output_dir $straindir/$bodysite/output --relaxed_parameters2 --nprocs_main 5 --clades $clade 1> $straindir/$bodysite/output/log_full.txt 2> $straindir/$bodysite/output/error_full.txt"))
{
    die "strainphlan.py ERROR: $!\n";
}

#Add metadata to the tree
#must of metadata file in strainphlan group directory
#multiple trees and multiple metadata files can be used (space separated, and wild card can be used)
#metadata file (tab separated, can have multiple columns)
my $metadata = "SubjectID"; #change based on what metadata you want listed on the tree
if (system("add_metadata_tree.py --ifn_trees $straindir/$bodysite/output/RAxML_bestTree.$clade\.tree --ifn_metadatas $straindir/$bodysite/$bodysite\_metadata.txt --metadatas $metadata"))
{
    die "Strainphlan add_metadata_tree.py ERROR: $!\n";
}

#Plot tree using Graphlan (graphlan scripts must be in path)
if (system("plot_tree_graphlan.py --ifn_tree $straindir/$bodysite/output/RAxML_bestTree.$clade\.tree.metadata --colorized_metadata $metadata --leaf_marker_size 60 --legend_marker_size 60"))
{
    die "Graphlan plot_tree_graphlan.py ERROR: $!\n";
}

#Create dendrogram using ggtree script
#breadcrumbs directory must be in path
if (system("strainphlan_ggtree_Mod.R $straindir/$bodysite/output/RAxML_bestTree.$clade\.tree $straindir/$bodysite/$bodysite\_metadata.txt $straindir/$bodysite/output/$clade\.fasta $straindir/$bodysite/output/$bodysite\_$clade\_ggtree_1.png $straindir/$bodysite/output/$bodysite\_$clade\_ggtree_2.png"))
{
    die "strainphlan_ggtree_Mod.R ERROR: $!\n";
}

#Create a distance matrix
if (system("distmat -sequence $straindir/$bodysite/output/$clade\.fasta -nucmethod 2 -outfile $straindir/$bodysite/output/$clade\.distmat"))
{
    die "distmat ERROR: $!\n"; 
}
if (system("strainphlan_ordination_Mod.R $straindir/$bodysite/output/$clade\.distmat $straindir/$bodysite/$bodysite\_metadata.txt $straindir/$bodysite/output/$bodysite\_$clade\_strainord.png"))
{
    die "strainphlan_ordination_Mod.R ERROR: $!\n";
}

