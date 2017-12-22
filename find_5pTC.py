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
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

__author__ = "Andres Aguilar"
__date__ = "14/Dic/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

warnings.filterwarnings("ignore")


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


def df2reads(dataframe, id_col="Name", seq_col="Seq", extra_data=""):
    """ Function: df2reads
    Convert pandas DataFrame into a SeqRecord list
    """
    result = list()
    for i in dataframe.index:
        name = dataframe[id_col][i]
        seq = dataframe[seq_col][i]
        desc = ""
        if extra_data != "":
            temp = extra_data.split(',')
            temp = dataframe[temp].ix[i]
            desc = "|".join([str(x) for x in temp.values.tolist()])
        result.append(SeqRecord(Seq(seq, IUPAC.ambiguous_dna), name=name, id=name, description=desc))
    return result


def write_fasta(reads, path_to_fasta):
    """ Function: write_fasta """
    return SeqIO.write(reads, path_to_fasta, "fasta")


def read_fasta(path_to_file):
    """ Function: read_fasta """
    if path.exists(path_to_file) and path.isfile(path_to_file):
        return SeqIO.parse(path_to_file, "fasta")
    else:
        return None


def main(fasta_file, out_file):
    df = reads2df(read_fasta(fasta_file))

    df["TC_Seq"] = df["Seq"].apply(lambda x: x[x.find("TC"):])
    write_fasta(df2reads(df, seq_col="TC_Seq"), out_file)


if __name__ == "__main__":
    args = argparse.ArgumentParser()

    args.add_argument("-i", "--input_file", help="Helitron sequences file", required=True)
    args.add_argument("-o", "--output_file", help="output file", required=True)
    p = args.parse_args()

    main(p.input_file, p.output_file)
