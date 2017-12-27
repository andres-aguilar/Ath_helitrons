#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  group_genes_helitrons.py
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

from __future__ import print_function

import argparse
import warnings
import pandas as pd

__author__ = "Andres Aguilar"
__date__ = "27/Dic/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

warnings.filterwarnings("ignore")


def split_helitron_list(x):
    """ Function: split_helitron_list """
    tmp = x.split("-")
    tmp = list([y.split(".")[0] for y in tmp])
    return tmp


def main(captures_file, genes_dist_file):
    # read captures file
    captures = pd.read_table(captures_file)
    captures.sort_values("NoHelitrons", ascending=False, inplace=True)
    captures.index = list(range(len(captures)))  # reindex

    captures.to_csv(captures_file, index=False, sep='\t')

    print("Captured genes", captures["Gene"].unique().size)

    # Cannonic helitrons with capture
    captures["Helitron_list"] = captures["Helitrons"].apply(split_helitron_list)

    helitrons_list = list()
    for index in captures.index:
        tmp_hels = list(set(captures["Helitron_list"][index]))
        helitrons_list.extend(tmp_hels)

    helitrons_list = list(set(helitrons_list))
    print("Number of cannonic helitrons with capture", len(helitrons_list))

    # Captures accumulation  by gene
    hels = list()
    genes = captures["Gene"].unique().tolist()
    for gene in genes:
        tmp = captures.get(captures["Gene"] == gene)
        hels.append(tmp["NoHelitrons"].sum())

    captured_genes_dist = pd.DataFrame()
    captured_genes_dist["Gene"] = genes
    captured_genes_dist["NoHelitrons"] = hels

    captured_genes_dist.sort_values("NoHelitrons", ascending=False, inplace=True)
    captured_genes_dist.to_csv(genes_dist_file, sep='\t', index=False)

    print("Total number of captures", captured_genes_dist["NoHelitrons"].sum())


if __name__ == "__main__":
    args = argparse.ArgumentParser()

    args.add_argument("-c", "--captures_file", help="Captures file", required=True)
    args.add_argument("-g", "--genes_file", help="Genes group output file", required=True)
    p = args.parse_args()

    main(p.captures_file, p.genes_file)
