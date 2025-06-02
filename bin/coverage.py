#!/usr/bin/env python3

import sys
from multiprocessing import Pool
import subprocess
import json
import pickle
import os
import argparse
import math
import re
import time
import glob

# --- Argument parsing and setup ---

version = "1.0"
parser = argparse.ArgumentParser()
parserO = parser._action_groups.pop()
parserR = parser.add_argument_group('required arguments')
parserO.add_argument('--out',   action='store', default='.',
					help='path of the output files and scratch directory (.)')
parserO.add_argument('-j', '--index', type=int, action='store',
					help='index of the LSF job array. One based')
parserO.add_argument('-k', '--max_index', type=int,
					action='store', help='maximum index of the LSF job array')
parserO.add_argument('-t', '--threads', type=int,
					action='store', default=1, help='number of threads (1)')
parserR.add_argument('-R', '--ref', action='store',
					required=True, help="reference sequence")
parserR.add_argument('-A', '--normal', action='store',
					required=True, help="normal BAM / CRAM")
parserR.add_argument('-B', '--duplex', action='store',
					required=True, help="duplex (tumour) BAM / CRAM")
parserO.add_argument('-v', '--version', action='version', version=version)
parser._action_groups.append(parserO)

# Subcommand arguments
subparsers = parser.add_subparsers(dest='subcommand', help='subcommands')
subparsers.required = True
parser_cov = subparsers.add_parser('cov', help='coverage calculation')
parser_covO = parser_cov._action_groups.pop()
parser_covR = parser_cov.add_argument_group('required arguments')
parser_covO.add_argument('--exclude', action='store', default='MT,GL%%,NC_%,hs37d5',
						help='List of contigs to exclude. Comma separated, %% acts as a wild card. (MT,GL%%,NC_%%,hs37d5)')
parser_covO.add_argument('--include', action='store',
						help='Only include these contigs. Comma separated, %% acts as a wild card.')
parser_covO.add_argument('--larger', type=int, action='store',
						default=0, help='Only include contigs larger than this size. (0)')
parser_covO.add_argument('-w', '--win', type=int, action='store',
						default=100, help='bin size for coverage distribution (100)')
parser_covO.add_argument('-Q', type=int, action='store',
						default=0, help="minimum mapQ to include a duplex read (0)")
parser_cov._action_groups.append(parser_covO)

args = parser.parse_args()

# --- Helper classes and functions ---

class GInterval:
	def __init__(self, chrr, beg, end):
		self.chr = chrr
		if (end < beg):
			raise ValueError("Interval %s: %s - %s is invalid!" %
							(chrr, beg, end))
		self.beg = beg
		self.end = end
		self.l = end - beg + 1

	def __lt__(self, other):
		return (self.chr, self.beg, self.end) < (other.chr, other.beg, other.end)

def runCommand(command):
	if (command is None):
		return
	for ijob in command.rstrip(';').split(';'):
		print(f"\nExecuting: {ijob}\n", file=sys.stdout)
		p = subprocess.Popen(
			ijob, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
		std: tuple = p.communicate()
		if (p.returncode != 0):
			error = std[1].decode()
			print(f"\n!Error processing:  {ijob}\n", file=sys.stderr)
			raise ValueError(error)
	return

def runBamcov(bam, mapQ, window, ichr, out):
	if (bam is None):
		return
	out_bed = out + ".cov.bed"
	runCommand("bamcov -q %s -w %s -r \"%s\" -o %s %s" %
			(mapQ, window, ichr, out_bed, bam))
	runCommand("bgzip -l 2 -f %s" % (out_bed))
	return

def getBAMcontigs(bam):
	contigs = []
	if (bam is None):
		return contigs
	job = "samtools view -H %s" % bam
	p = subprocess.Popen(
		job, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	while True:
		output = p.stdout.readline().decode()
		if output == '' and p.poll() is not None:
			break
		if (re.match("@SQ", output)):
			contigs.append(output.strip())
	rc = p.poll()
	if (rc != 0):
		error = p.stderr.read().decode()
		print(f"\n!Error processing:  {job}\n", file=sys.stderr)
		raise ValueError(error)
	return contigs

# --- Coverage section ---

# For array execution, stagger access for files
if (args.index is not None):
	time.sleep(args.index)

if (args.index is None or args.index == 1):
	with open("args.json", "w") as jsonOut:
		json.dump(args.__dict__, jsonOut)

# --- Coverage section logic starts here ---
if (args.subcommand == 'cov'):
	# Build the chromosome dictionary, list and intervals
	if (args.exclude is None or args.exclude == ""):
		excludes = []
	else:
		excludes = [re.compile(istr + "$")
					for istr in args.exclude.replace("%", ".+").split(',')]
	if (args.include is None or args.include == ""):
		includes = [re.compile(".+")]
	else:
		includes = [re.compile(istr + "$")
					for istr in args.include.replace("%", ".+").split(',')]
	chrList = []
	rnames = {}
	with open(args.ref + '.fai', 'r') as iofile:
		for iline in iofile:
			ichr = iline.split('\t')[0]
			if any(iregx.match(ichr) for iregx in includes):
				if any(iregx.match(ichr) for iregx in excludes):
					continue
				ilength = int(iline.split('\t')[1])
				if (ilength <= args.larger):
					continue
				chrList.append(ichr)
				rnames[ichr] = ilength

	# BAMs vs reference sanity checks for contigs
	bamNContigs = {}
	bamNOrder = []
	for iline in getBAMcontigs(args.normal):
		ichr = iline.split('\t')[1].replace('SN:', '')
		ilength = int(iline.split('\t')[2].replace('LN:', ''))
		bamNContigs[ichr] = ilength
		bamNOrder.append(ichr)

	bamTContigs = {}
	bamTOrder = []
	for iline in getBAMcontigs(args.duplex):
		ichr = iline.split('\t')[1].replace('SN:', '')
		ilength = int(iline.split('\t')[2].replace('LN:', ''))
		bamTContigs[ichr] = ilength
		bamTOrder.append(ichr)

	for icontig in rnames:
		if (not (icontig in bamNContigs)):
			sys.exit(f"Reference contig {icontig} was not found in normal BAM")
		if (bamNContigs[icontig] != rnames[icontig]):
			sys.exit(f"Length of contig {icontig} in normal BAM doesn't match reference ({bamNContigs[icontig]} vs {rnames[icontig]})")

	for icontig in rnames:
		if (not (icontig in bamTContigs)):
			sys.exit(f"Reference contig {icontig} was not found in duplex BAM")
		if (bamTContigs[icontig] != rnames[icontig]):
			sys.exit(f"Length of contig {icontig} in duplex BAM doesn't match reference ({bamTContigs[icontig]} vs {rnames[icontig]})")

	# If order of BAMs is not the same, mpilupe bogus results
	for icontig in rnames:
		if (bamNOrder.index(icontig) != bamTOrder.index(icontig)):
			sys.exit("Contigs in BAM files must have the same order ( check order in headers )")

	gintervals = []
	for ichr in chrList:
		gintervals.append(GInterval(ichr, 2, rnames[ichr]-1))
	gintervals.sort()

	reorderchr = []
	for iint in gintervals:
		reorderchr.append(iint.chr)
	chrList = reorderchr
	print("\nAnalysing contigs: %s\n" % chrList, file=sys.stdout)
	print("Starting cov calculation\n", file=sys.stdout)

	inputs = []
	for ii, ichr in enumerate(chrList):
		cov_prefix = f"{ii + 1}"
		inputs.append((args.duplex, str(args.Q), str(args.win), ichr, cov_prefix))
	if (args.index is None):
		with open("gIntervals.dat", 'wb') as iofile:
			pickle.dump(gintervals, iofile)
		with Pool(args.threads) as p:
			p.starmap(runBamcov, inputs)
	else:
		if (args.index == 1):
			with open("gIntervals.dat", 'wb') as iofile:
				pickle.dump(gintervals, iofile)

		commands = [[] for i in range(args.max_index)]
		jj = 0
		for ii in range(len(inputs)):
			if (ii % args.max_index == 0):
				jj = 0
			commands[jj].append(
				"runBamcov(inputs[%s][0], inputs[%s][1], inputs[%s][2], inputs[%s][3], inputs[%s][4])" % (ii, ii, ii, ii, ii))
			jj += 1
		for icmd in commands[args.index - 1]:
			exec(icmd)

	# Concatenate all .cov.bed.gz files into cov.bed.gz and remove the originals
	cov_files = sorted(glob.glob("*.cov.bed.gz"), key=lambda x: int(x.split('.')[0]))
	with open("cov.bed.gz", "wb") as outfile:
		for fname in cov_files:
			with open(fname, "rb") as infile:
				outfile.write(infile.read())
	for fname in cov_files:
		os.remove(fname)

	print("Completed cov calculation\n", file=sys.stdout)
