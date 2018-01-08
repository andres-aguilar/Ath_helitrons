#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  classify_captures.py
#
#  Copyright 2018 Andres Aguilar <andresyoshimar@gmail.com>
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

__author__ = "Andres Aguilar"
__date__ = "27/Dic/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

warnings.filterwarnings("ignore")
GFF_COLUMNS = ["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]


def read_gff(gff_file, comment='#'):
    """ Function: read_gff """
    if path.isfile(gff_file):
        return pd.read_table(gff_file, header=None, names=GFF_COLUMNS, comment=comment)


def get_chromosome(x):
    if "AT1" in x:
        return "Chr1"
    elif "AT2" in x:
        return "Chr2"
    elif "AT3" in x:
        return "Chr3"
    elif "AT4" in x:
        return "Chr4"
    else:
        return "Chr5"


def main(captures_file, gff_file, helitrons_file, classification_file):
    captures = pd.read_table(captures_file)
    genes = captures["Gene"].unique().tolist()
    captures["Helitron_list"] = captures["Helitrons"].apply(lambda x: x.split("-"))
    captures["Helitron_list"] = captures["Helitron_list"].apply(lambda x: [a.split(".")[0] for a in x])

    hels = list()
    gn_list = list()
    cap_range = list()
    class_list = list()
    for gene in genes:
        tmp = captures.get(captures["Gene"] == gene)
        for x in tmp.itertuples():
            for j in x.Helitron_list:
                gn_list.append(gene)
                hels.append(j)
                cap_range.append(tmp["CapRange"][tmp.first_valid_index()])
                class_list.append(tmp["Class"][tmp.first_valid_index()])

    df = pd.DataFrame()
    df["Gene"] = gn_list
    df["Helitron"] = hels
    df["CapRange"] = cap_range
    df["Subclass"] = class_list

    df["Chr_gene"] = df["Gene"].apply(get_chromosome)
    df["Chr_hel"] = df["Helitron"].apply(get_chromosome)

    df["Same_chr"] = df.apply(lambda x: x.Chr_gene == x.Chr_hel, axis=1)

    # Get genomic coordinates for helitrons and genes
    annotation = read_gff(gff_file)

    annotation = annotation[annotation["Feature"].str.contains(r"^gene$")]

    annotation["Id"] = annotation["Attribute"].apply(lambda x: x.split(";")[0])
    annotation["Id"] = annotation["Id"].apply(lambda x: x.split("=")[-1])

    df["Gene_start"] = 0
    df["Gene_end"] = 0
    for x in df.itertuples():
        g_tmp = annotation.get(annotation["Id"] == x.Gene)
        df["Gene_start"][x.Index] = g_tmp["Start"][g_tmp.first_valid_index()]
        df["Gene_end"][x.Index] = g_tmp["End"][g_tmp.first_valid_index()]

    helitrons = pd.read_table(helitrons_file)
    helitrons["Name"] = helitrons["Long_name"].apply(lambda x: x.split("|")[0])
    helitrons["Start"] = helitrons["Long_name"].apply(lambda x: int(x.split("|")[2]))
    helitrons["End"] = helitrons["Long_name"].apply(lambda x: int(x.split("|")[3]))

    df.index = df["Helitron"].tolist()
    helitrons.index = helitrons["Name"].tolist()

    df["Hel_start"] = helitrons["Start"].astype(int)
    df["Hel_end"] = helitrons["End"].astype(int)

    df.index = list(range(len(df)))

    # ##############################################################################
    # Gene capture classification
    #
    # Class I   - Gene capture
    #   a) Gene inside a helitron.
    #
    #   b) Gene fragment inside a helitron.
    #
    # Class II  - Helitron inside gene
    #
    # Class III - Helitron - gene intersection
    #
    # #############################################################################

    df["Class"] = ""
    for x in df.itertuples():
        starts = x.Gene_start >= x.Hel_start <= x.Hel_end
        ends = x.Gene_end >= x.Hel_start <= x.Hel_end

        if x.Gene_start >= x.Hel_start and x.Gene_end <= x.Hel_end:
            df["Class"][x.Index] = "Ia"
        elif x.Hel_start >= x.Gene_start and x.Hel_end <= x.Gene_end:
            df["Class"][x.Index] = "II"
        elif starts or ends:
            df["Class"][x.Index] = "III"
        else:
            df["Class"][x.Index] = "Ib"

    print(df["Class"].value_counts())

    df.to_csv(classification_file, index=False, sep='\t')


if __name__ == "__main__":
    args = argparse.ArgumentParser()

    args.add_argument("-i", "--captures_file", help="Captures file", required=True)
    args.add_argument("-g", "--gff_file", help="Gff file", required=True)
    args.add_argument("-f", "--helitrons_file", help="file with corrected headers", required=True)
    args.add_argument("-c", "--classification_file", help="output file", required=True)
    p = args.parse_args()

    main(p.captures_file, p.gff_file, p.helitrons_file, p.classification_file)
