#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  filter_helend_GC30.py
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

from Bio.SeqUtils import GC

# Helendout.txt columns
helend_out = ["Name", "Strand", "Start_helend", "Helend", "Hairpin",
              "Hairpin_len", "Gap_pos", "Mismatches"]


def gc_percent(x):
    """ Calculate %GC """
    stem1, loop, stem2 = x.split("*")
    return GC(stem1)


def main(helends_file):
    # Read table with helitrons structural information
    helends = pd.read_table(helends_file, header=None, names=helend_out)

    helends["GC"] = helends["Hairpin"].apply(gc_percent)
    helends = helends.get(helends["GC"] >= 30)

    helends.to_csv(helends_file, sep='\t', index=False, header=False)


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("-helends", "--helends", help="Helends file", required=True)
    p = args.parse_args()

    main(p.helends)