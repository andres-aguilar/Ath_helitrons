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
__date__ = "27/Dic/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

warnings.filterwarnings("ignore")
GFF_COLUMNS = ["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]


def read_gff(gff_file, comment='#'):
    """ Function: read_gff """
    if path.isfile(gff_file):
        return pd.read_table(gff_file, header=None, names=GFF_COLUMNS, comment=comment)


def read_fasta(path_to_file):
    """ Function: read_fasta """
    if path.exists(path_to_file) and path.isfile(path_to_file):
        return SeqIO.parse(path_to_file, "fasta")
    else:
        return None


def reads2df(reads, extra_info=False):
    """ Function: reads2df
    Convert a list of SeqRecords into a pandas DataFrame
    """
    names = list()
    seqs = list()
    desc = list()
    for read in reads:
        name = read.id
        names.append(name)
        seqs.append(str(read.seq))
        desc.append(read.description)

    if extra_info:
        df = pd.DataFrame({"Name": names, "Seq": seqs})
        df["Description"] = desc
    else:
        df = pd.DataFrame({"Name": names, "Seq": seqs})

    return df


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


def main(captures, gff_file, helitrons_fasta, classification_file):
    genes = captures["Gene"].unique().tolist()

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

    helitrons = reads2df(read_fasta(helitrons_fasta))
    helitrons["Len"] = helitrons["Seq"].apply(lambda x: len(x))
    helitrons["Start"] = helitrons["Name"].apply(lambda x: int(x.split("|")[2]))
    helitrons["End"] = helitrons["Name"].apply(lambda x: int(x.split("|")[3]))

    df.index = df["Helitron"].tolist()

    df["Hel_start"] = helitrons["Start"]
    df["Hel_end"] = helitrons["End"]

    df.index = list(range(len(df)))

    # ##############################################################################
    # Clasificación de capturas genicas
    #
    # Clase I   - Gene capture
    #   a) Gene inside a helitron.
    #
    #   b) Gene fragment inside a helitron.
    #
    # Clase II  - Helitron inside gene
    #
    # Clase III - Helitron - gene intersection
    #
    # #############################################################################

    df["Class"] = ""
    for x in df.itertuples():
        starts = x.Gene_start >= x.Hel_start and x.Gene_start <= x.Hel_end
        ends = x.Gene_end >= x.Hel_start and x.Gene_end <= x.Hel_end

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

    args.add_argument("-i", "--input_file", help="Helitron sequences file", required=True)
    args.add_argument("-o", "--output_file", help="output file", required=True)
    p = args.parse_args()

    # main(p.input_file, p.output_file)
