#!/usr/bin/env python2.7

import gzip,time
import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="")
parser.add_argument("--output", help="")
args = parser.parse_args()

def reduce(input, output):
    for line in open(input):
        print line
        time.sleep(3)