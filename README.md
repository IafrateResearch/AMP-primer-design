# AMP-primer-design

**AMP Primer Design** is an automated pipeline which can be used to design primers for detection of gene fusion (targeted RNA-seq) and mutations (targeted DNA-seq for SNV, indel and CNV) using anchored multiplex PCR.

## Reference publication
[Zheng Z, et al. Anchored multiplex PCR for targeted next-generation sequencing. Nat Med. 2014](http://www.nature.com/nm/journal/v20/n12/full/nm.3729.html)
### Date created: 2012-04-13

## How does AMP-primer-design work?  

The design involves 5 computational steps:

1. Generate template sequences, both sense and antisense for the 'fusion' assay, and plus and minus strands for the 'mutation' assay. The sequences are then masked as 'N' for every dbSNPv137 locus. Masked sequences will be output to seq/ folder, and unmasked sequences to seq.noMask/ folder.


2. Use Primer3 to design candidate primers. Primers are filtered against NGS adaptors including 96 Illumina N5s, 12 N7s, 96 IonTorrent adaptors and human repetitive sequences. Up to 6 design iterations will be attempt. The first 3 rounds of primer design will use masked template sequences, with decreased design stringency. If no candidate primers found, unmasked sequences will be used as template for a further 3 rounds of design. If still no candidate primers, the target will be reported.

3. Map candidate primers against human genome using BLAT. BLAT is a better tool than some other tools (e.g. BLAST, BWA) at mapping short sequence (primer) because of its popularity in gap alignment, which is useful for mapping cDNA reads (primer) against genomic reference. Thus, filtered primers are essentially filtered against all possible alterative splicing variants. Those candidate primers with the last 12 bases at the 3' end map more than 5 (arbitrary) locations in the genome will be discarded.

4. Pairing candidates for GSP1 and GSP2 based on pair penalty scores. Briefly, the penalty takes into account the distance between GSP2 3' end and target site, avoid the same mispriming targets for GSP1 and GPS2, the chromosomal distance and annealing T difference between GSP1 and GSP2, etc. Candidate pairs are ranked by penalty.

5. Check uniqueness of 12 bases at 3' of all primers. If not unique, alternative candidate pairs (ranked 2, 3...) from 'ranked.all.candidate.pairs' will be retrieved until all tail 12 bases are unique.

Lastly, add GSP tags, generate output files, save primer.bed, save final tail 12 bases of all primers for future use (e.g, to check uniqueness when adding new primers to an existing panel).

## Installation

Required software, libraries and data dependency

- python 2.7
- Primer3 (http://primer3.sourceforge.net/primer3_manual.htm)
- BLAT (https://users.soe.ucsc.edu/~kent/src/)
- R libraries: BSgenome.Hsapiens.UCSC.hg19
- The dependency data (e.g. ~/dependency-AMP/)
- To prepare the refGene transcript and exons coordinates:

    > wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

    > gunzip refGene.txt.gz

    > cut -f 2- refGene.txt | sort -k1,1b > refGene.sort 

    > cat /path/to/AMP-primer-design/refGene.header refGene.sort > /path/to/dependency-AMP/hg19.RefGene.NM

- To prepare the mispriming sequences used for filtering primers during primer3 design. It consists of NGS adaptors (supplied by this pipeline 'NGSadaptors.fa') and Human Repetitive sequences (available [here](http://www.girinst.org/server/RepBase/index.php), requires registration - free for non-profit).

    > cat NGSadaptors.fa humrep.fasta > /path/to/dependency/NGSadaptors_and_humanRep.fa

- To prepare dbSNP and clinically relevant SNPs from the 1000 Genomes Project (20120626 Release), run make.dbsnp.database.sh (please note the reference version GRCh37 or GRCh38):

    > bash make.dbsnp.database.sh


The dependency data (e.g. in '~/dependency-AMP') should contain:

| Filename    | Content |
| :---------: |:-------:|
| hg19.2bit   | Contains the complete hg19 Human Genome in 2bit format, available [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) |
| hg19.RefGene.NM | Contains the refGene file prepared above. |
| NGSadaptors_and_humanRep.fa | Contains mispriming sequences used for filtering primers during primer3 design. |
| dbsnp/  | dbSNP and clinically relevant SNPs from the 1000 Genomes Project. Prepared above. |


## Usage examples
Checkout help doc:

    > python /path/to/AMP-primer-design/design.py --help


### Usage examples 1 - a 'fusion' panel

1. Create a panel gene list file with RefSeq transcript IDs "NM_". See 'example/lung.fusion.genelist.txt' for example. Alternatively, you can use your own template sequences (e.g. virus) under the ~/project-AMP/$panel/seq/ folder (under development).

2. Run the command (see also 'example.command.txt')

    > python /path/to/AMP-primer-design/design.py  
	--assaytype="fusion"  
	--panel="lung.fusion"  
	--genelist="~/project-AMP/lung.fusion.genelist.txt"  
	--pjdir="~/project-AMP"   
	--depdir="~/dependency-AMP"  
	--ampdir="~/repo/AMP-primer-design"  
	--primer3path="~/primer3-2.3.5"  
	--blatdir="~/bin/x86_64"  
	--keep_gfSvr=1  
	--GSP1tag="GGATCTCGACGCTCTCCCT"  
	--GSP2tag="CCTCTCTATGGGCAGTCGGTGAT"  
	--ncpu=8  
	--tempsize=90  
	--NGSadaptors_and_humanRep="NGSadaptors_and_humanRep.fa"  


### Usage example 2, a 'mutation' panel

1. Create a panel gene list file with RefSeq transcript IDs "NM_". See 'example/cancer.v1.genelist.txt' for example. Alternatively, you can use your own template sequences (e.g. virus) under the ~/project-AMP/ panel/seq/ folder (under development).

2. Run the command (see also 'example.command.txt')

    > python /path/to/AMP-primer-design/design.py  
	--assaytype="mutation"  
	--panel="cancer.v1"  
	--genelist="~/project-AMP/cancer.v1.genelist.txt"  
	--pjdir="~/project-AMP"  
	--depdir="~/dependency-AMP"  
	--ampdir="~/repo/AMP-primer-design"  
	--primer3path="~/primer3-2.3.5"  
	--blatdir="~/bin/x86_64"  
	--keep_gfSvr=1  
	--GSP1tag="GGATCTCGACGCTCTCCCT"  
	--GSP2tag="CCTCTCTATGGGCAGTCGGTGAT"  
	--ncpu=32  
	--tempsize=100  
	--NGSadaptors_and_humanRep="NGSadaptors_and_humanRep.fa"  
	--subExonSize=300  
	--leadsize=3  


### Usage example 3, using 'bed' file specified by user  

1. Generate a standard BED (e.g. cidstr.bed) for target regions (note the reference version GRCh37 or GRCh38)

2. Run the command (see also 'example.command.txt')

    > python /path/to/AMP-primer-design/design.py  
	--assaytype="bed"  
	--panel="cidstr"  
	--genelist="~/project-AMP/cidstr.bed"  
	--pjdir="~/project-AMP"  
	--depdir="~/dependency-AMP"  
	--ampdir="~/repo/AMP-primer-design"  
	--primer3path="~/primer3-2.3.5"  
	--blatdir="~/bin/x86_64"  
	--keep_gfSvr=1  
	--GSP1tag="GGATCTCGACGCTCTCCCT"  
	--GSP2tag="CCTCTCTATGGGCAGTCGGTGAT"  
	--ncpu=32  
	--tempsize=100  
	--NGSadaptors_and_humanRep="NGSadaptors_and_humanRep.fa"  
	--subExonSize=300  
	--leadsize=3   
 
