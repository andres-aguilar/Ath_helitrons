#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  find_5pTC.py
#
#  Copyright 2017 Andres Aguilar <andresyoshimar@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

import argparse
import warnings
import pandas as pd

from os import path
from Bio import SeqIO

__author__ = "Andres Aguilar"
__date__ = "14/Dic/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

warnings.filterwarnings("ignore")

helend_final = ["Name", "Strand", "Start_helend", "Helend", "Hairpin",
                "Hairpin_len", "Gap_pos", "Mismatches", "Hairpin_start",
                "Hairpin_end", "CTRR_start", "CTRR_end"]


def _get_coordinates(helend_file):
    coords = dict()
    with open(helend_file, "r") as helend:
        while True:
            line = helend.readline()
            if not line:
                break
            line = line.strip()
            line = line.split("\t")
            coords[line[0]] = int(line[-1])
    return coords


def _load_helitrons(helitrons_file):
    return SeqIO.to_dict(SeqIO.parse(helitrons_file, "fasta"))


def _cut_sequence(sequence,  end):
    return sequence[:end]


def cut_cannonical_sequence(helend_file, helitrons_file):
    """ Function: cut_cannonical_sequence

    Params
     ------
       helend_file    : File with helend information
       helitrons_file : Fasta file with helitron sequences
    """
    if not path.isfile(helend_file):
        raise Exception("helend_file must exist!")

    if not path.isfile(helitrons_file):
        raise Exception("helitrons_file must exist!")

    cannonic_helitrons = dict()
    coords = _get_coordinates(helend_file)
    hels = _load_helitrons(helitrons_file)

    for key in coords.keys():
        ctrr_end = coords.get(key)

        if key in hels.keys():
            sequence = str(hels.get(key).seq)
            cannonic_seq = _cut_sequence(sequence, ctrr_end)
            cannonic_helitrons[key] = (cannonic_seq, ctrr_end)

    return cannonic_helitrons


def main(helends_file, helitrons_file, out_file):
    helends = pd.read_table(helends_file, header=None, names=helend_final)

    # Sort values highest CTRR coord first
    helends.sort_values(["Name", "CTRR_end"], inplace=True, ascending=[True, False])
    helends.drop_duplicates("Name", inplace=True)
    helends.to_csv(helends_file, sep="\t", index=False, header=False)

    hels = cut_cannonical_sequence(helends_file, helitrons_file)

    # Write Fasta file with cannonical sequence
    with open(out_file, "w") as out:
        for key in hels.keys():
            element = hels.get(key)
            out.write(">{0}\n{1}\n".format(key, element[0]))


if __name__ == "__main__":
    args = argparse.ArgumentParser()

    args.add_argument("-i", "--helends_file", help="Helend coordinates file", required=True)
    args.add_argument("-f", "--helitrons_file", help="Fasta with TC helitron sequences", required=True)
    args.add_argument("-o", "--output_file", help="output file", required=True)
    p = args.parse_args()

    main(p.helends_file, p.helitrons_file, p.output_file)