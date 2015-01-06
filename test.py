#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='AMP primer design input arguments')
    parser.add_argument('--assaytype', default='fusion')
#    parser.add_argument('--genelist', required=True)
#    parser.add_argument('--panel', required=True)
    parser.add_argument('--pjdir', required=True)
#    parser.add_argument('--ampdir', required=True)
#    parser.add_argument('--blatdir', required=True)
#    parser.add_argument('--primer3path', required=True)
    parser.add_argument('--keep_gfSvr', type=int, default=1)
    parser.add_argument('--ncpu', type=int, default=4)
#    parser.add_argument('--tempsize', type=int, default=90)
#    parser.add_argument('--subExonSize', type=int, default=300)
#    parser.add_argument('--leadsize', type=int, default=3)
    parser.add_argument('--utr', type=int, default=0)
    parser.add_argument('--GSP1tag', default='GGATCTCGACGCTCTCCCT')
    parser.add_argument('--GSP2tag', default='CCTCTCTATGGGCAGTCGGTGAT')
    parser.add_argument('--NGSadaptors_and_humanRep', default='NGSadaptors_and_humanRep.fa')

    args = vars(parser.parse_args())

    return args

args = parseArgs()
print(args)
config = open("config.txt", "w+")
for arg in args:
	print(arg, file = config, sep="\n")

if __name__ == '__main__':
    
    args = parseArgs()
    if not os.path.exists(args['pjdir']):
        os.makedirs(args['pjdir'])

os.chdir(args['pjdir'])
