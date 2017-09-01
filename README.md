# HIDskinmicrobiome

Documentation for scripts and analysis from Schmedes et. al 2017, (AEM, under review)
using data from Oh et al., Cell, 2016 for shotgun metagenomic analyses

## General Workflow

- Fastq readpreprocessing
- Taxonomic classification using MetaPhlAn2
- Strain characterization using StrainPhlAn and PanPhlAn
- Statistical classification using pangenome gene presence/absence
- Statistical classification using nucleotide diversity of clade-specific markers

### General files for all workflows  
Genomix.pm is a perl module that is used for all scripts in this repository. This module relies on the OhMeta.txt file
compilied with metadata from the Sequence Read Archive and supplemental material from Oh et al. as well as
database metadata from MetaPhlan2 from pkl file. There are 3 main functions in Genomix.pm.   

queryDB() . 
Uses 2 or 4 arguments and returns an array of of values from the last argument category.  
If 2 arguments are used: $metafile = "OhMeta.txt".  
@srrs = queryDB($metafile, "srrID"); returns an array of all srrIDs from OhMeta.txt . 
If 4 arguments are used:  
@srrs = queryDB($metafile, "Mb", "bodysiteID", "srrID"); returns an array of all srrIDs for all samples from the manubrium bodysite.  

gettaxa() . 
Uses no arguments and returns 3 hashes with MetaPhlAn2 marker info . 
$taxamarker{$markername} = $fulltaxonpathformarker . 
$taxalevel{$fulltaxonpath} = "terminalclade\tcladelevel" . 
$markerlength{$markername} = $markerlengthinDB . 

samplename(@SRS) . 
Uses 1 argument, an array of SRSids. If only element is in the array a scalar will be returned with the associated samplename (HV01_Bodysite_Timepoint) . 
If the array has more than 1 element a hash will be returned. Keys=SRSid, values=samplename . 

