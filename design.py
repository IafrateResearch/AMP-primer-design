#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='AMP primer design input arguments')
    parser.add_argument('--assaytype', default='fusion', help="'fusion' or 'mutation'. \
		    Fusion assay will retrieve exonic sequence template for primer \
		    design, and mutation assay intronic template.")
    parser.add_argument('--genelist', required=True, help="name (with path) of the gene list, see\
		    'example/lung.fusion.genelist.txt' for example.")
    parser.add_argument('--depdir', required=True, help="path to dependency data.")
    parser.add_argument('--panel', required=True, help="the name of panel, e.g. 'lung.fusion'.")
    parser.add_argument('--pjdir', required=True, help="the project folder. A project/panel \
		    folder, e.g. '~/project-AMP/lung.fusion' will be created by the pipeline.")
    parser.add_argument('--ampdir', required=True, help="path to AMP primer design.")
    parser.add_argument('--blatdir', required=True, help="path to BLAT.")
    parser.add_argument('--primer3path', required=True, help="path to Primer3.")
    parser.add_argument('--keep_gfSvr', type=int, default=1, help="keep BALT gfSever in memory.")
    parser.add_argument('--ncpu', type=int, default=4, help="number of available CPUs for \
		    multi-threading.")
    parser.add_argument('--tempsize', type=int, default=90, help="size of template sequence to\
		    retrieve from genome and to design primers on. Default 90 bp considers \
		    degraded RNA in FFPE samples.")
    parser.add_argument('--subExonSize', type=int, default=300, help="size for tiling targets.\
		    For exons larger than subExonSize, the exons will be divided into subExons\
		    with max size of subExonSize a, b, c...")
    parser.add_argument('--leadsize', type=int, default=3, help="the chromosomal distance between\
		    template sequence and target location. This is to avoid GSP2 ends in exon\
		    boundary.")
    parser.add_argument('--utr', type=int, default=0, help="whether or not to target UTR.")
    parser.add_argument('--GSP1tag', default='GGATCTCGACGCTCTCCCT', help="the tag to be appended\
		    to 5' end of GSP1 primers. This tag do not participate in sequencing.")
    parser.add_argument('--GSP2tag', default='CCTCTCTATGGGCAGTCGGTGAT', help="the tag to be \
		    appended to 5' end of GSP2 primers. For Illumina, this tag is \
		    Read2 Sequencing Primer. Default here is Ion Torrent (P23) sequence.\
		    This tag allows for the same primers (hundreds to thousands)\
		    to be used for both Ion Torrent and Illumina platforms.\
		    (For Illumina Miseq, if use this GSP2 tag, in wet-lab:\
		    a. Add 3 ul of 100 uM of Illumina.custom.Index1.sequencing.primer\
		    to Miseq Reagent cartridge position 13 (Index Primer Mix)\
		    b. Add 3 ul of 100 uM of Illumina.custom.Read2.sequencing.primer\
		    to Miseq Reagent cartridge position 14 (Read 2 Primer Mix).\
		    See NGSadaptors.fa for the above primer sequences)")
    parser.add_argument('--NGSadaptors_and_humanRep', default='NGSadaptors_and_humanRep.fa')

    args = vars(parser.parse_args())

    return args



if __name__ == '__main__':
    
    args = parseArgs()
    config = open("config.txt", "w+")
    print(args, file = config, sep="\n")
    config.close()

design = "Rscript " + args['ampdir'] + "/main.R"
os.system(design)

## END
