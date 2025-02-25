#!/usr/bin/env python3

########## LICENCE ##########
# Copyright (c) 2022 Genome Research Ltd
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of NanoSeq.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
###########################

import argparse
import itertools
import gzip

"""
Extract tags from an IDT adaptor tagged library. 

Read two fastq files, one for R1 and one for R2. Iterate through paired-end
reads (in synchrony). Write new R1 and R2 fastq files, where tags are extracted
from the read sequence to rb and mb auxiliary tags.


***General read structure
  NNN degenerate 3-nucleotide tag
  T   overhang for ligation
  X   spacer designed to avoid 'blooming' issues during sequencing
  R   non-adaptor read


***HpyCH4V/AluI library

  Possible sequences NNNTCAR(145), NNNTCTR(145), NNNXTCAR(144),
  NNNXTCTR(144)

  Divide reads into 3 segments
    3 tag nucleotides
    4 skips [TCAR, TCTR, XTCA, XTCT]
    144 non-adaptor read

***Sonicated library

  Possible sequences NNNT(R147) and NNNXT(R146)

  Divide reads into 3 segments:
    3 tag nucleotides
    2 skips [TR or XT]
    146 non-adaptor read
"""


class Read:
    """Class representing one read from a paired-end read."""

    def __init__(self, read_length):
        """Initializes with an empty header, sequence, quality, and tag."""
        self.read_length = read_length
        self.header = None
        self.sequence = None
        self.quality = None
        self.tag = None

    def add_header(self, header):
        """Add header string minus i5/i7 barcode."""
        self.header = header.strip().split(' ')[0]

    def add_sequence(self, sequence):
        """Add sequence string."""
        self.sequence = sequence.strip()

    def add_quality(self, quality):
        """Add quality string."""
        self.quality = quality.strip()

    def add_tag(self, matches):
        """Add first [matches] bases to tag."""
        self.tag = self.sequence[:matches]
        assert len(self.tag) == matches

    def trim_sequence(self, skips):
        """Trim adaptor from sequence string."""
        self.sequence = self.sequence[skips:]
        assert len(self.sequence) == self.read_length - skips

    def trim_quality(self, skips):
        """Trim adaptor from quality string."""
        self.quality = self.quality[skips:]
        assert len(self.quality) == self.read_length - skips

    def amend_header(self, tag):
        """Add rb and mb auxiliary tags to the header."""
        self.header = f'{self.header}\trb:Z:{self.tag}\tmb:Z:{tag}'


class ExtractTags:
    def __init__(self, in1, in2, out1, out2, matches, skips, read_length):
        """Initiate extract tags object."""
        self.matches = matches
        self.tags = self.canonical_tags()
        self.skips = matches + skips
        self.read_length = read_length

        self.open_files(in1, in2, out1, out2)
        self.extract_tags()

    def open_files(self, in1, in2, out1, out2):
        """Open input and output files."""
        mode_in = 'rt' if in1.endswith('.gz') else 'r'
        mode_out = 'wt' if out1.endswith('.gz') else 'w'
        self.in1 = gzip.open(in1, mode_in) if in1.endswith('.gz') else open(in1, mode_in)
        self.in2 = gzip.open(in2, mode_in) if in2.endswith('.gz') else open(in2, mode_in)
        self.out1 = gzip.open(out1, mode_out, compresslevel=2) if out1.endswith('.gz') else open(out1, mode_out)
        self.out2 = gzip.open(out2, mode_out, compresslevel=2) if out2.endswith('.gz') else open(out2, mode_out)

    def canonical_tags(self):
        """Returns a list of tags of size matches that consists of canonical A, C, G, T bases."""
        tags = [''.join(i) for i in itertools.product('ACGT', repeat=self.matches)]
        assert len(tags) == 4 ** self.matches
        return tags

    def extract_tags(self):
        """Extract tags and re-write fastq files."""
        try:
            for header1 in self.in1:
                r1 = Read(self.read_length)
                r2 = Read(self.read_length)
                r1.add_header(header1)
                r2_header = self.in2.readline()
                r2.add_header(r2_header)
                if r1.header != r2.header:
                    raise ValueError(f"Mismatched headers:\nR1: {r1.header}\nR2: {r2_header.strip()}")

                r1.add_sequence(self.in1.readline())
                r2.add_sequence(self.in2.readline())
                _ = self.in1.readline(), self.in2.readline()  # ignore spacer lines
                r1.add_quality(self.in1.readline())
                r2.add_quality(self.in2.readline())
                r1.add_tag(self.matches)
                r2.add_tag(self.matches)

                if r1.tag in self.tags and r2.tag in self.tags:
                    r1.trim_sequence(self.skips)
                    r2.trim_sequence(self.skips)
                    r1.trim_quality(self.skips)
                    r2.trim_quality(self.skips)
                    r1.amend_header(r2.tag)
                    r2.amend_header(r1.tag)
                    self.out1.write(f"{r1.header}\n{r1.sequence}\n+\n{r1.quality}\n")
                    self.out2.write(f"{r2.header}\n{r2.sequence}\n+\n{r2.quality}\n")
            if self.in2.readline():
                raise ValueError('R2 file has additional reads not present in R1')
        finally:
            self.in1.close()
            self.in2.close()
            self.out1.close()
            self.out2.close()


if __name__ == "__main__":
    msg = ("Extract tags from fasta files. For sonicated library use -m 3 -s "
           "2. For HpyCH4V/AluI library use -m 3 -s 4.")
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-a', required=True, help="Input fastq file R1")
    parser.add_argument('-b', required=True, help="Input fastq file R2")
    parser.add_argument('-c', required=True, help="Output fastq file R1")
    parser.add_argument('-d', required=True, help="Output fastq file R2")
    parser.add_argument('-m', type=int, required=True, help="Match length")
    parser.add_argument('-s', type=int, required=True, help="Skip length")
    parser.add_argument('-l', type=int, required=True, help="Read length")
    args = parser.parse_args()
    ExtractTags(args.a, args.b, args.c, args.d, args.m, args.s, args.l)
