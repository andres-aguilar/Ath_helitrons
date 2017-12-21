#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  get_helitron_sequences.py
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

from os import path
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, reverse_complement

__author__ = "Andres Aguilar"
__date__ = "14/Sep/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

warnings.filterwarnings("ignore")
_gff_columns = ["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]


def read_gff(gff_file, comment='#'):
    """ Function: read_gff """
    if path.isfile(gff_file):
        return pd.read_table(gff_file, header=None, names=_gff_columns,
                            comment=comment)


def get_sequence_from_chromosomes(elements, genome, start_col="Start", end_col="End", strand_col="Strand",
                                  chr_col="Chr", new_col="Seq"):
    """ Function: get_sequence_from_chromosomes

    Params
    -------
     elements : pandas DataFrame with the elements to extract
                ID, Start and End columns must exist
     genome   : fasta file with chromosomes

    Return
    ------
     elements DataFrame with sequence column
    """
    if not path.isfile(genome):
        print("{0} is not a file!".format(genome))
        return None

    chrs = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    elements[new_col] = ""
    for index in elements.index:
        start = elements[start_col][index]-1  # Python has 0 based arrays
        end = elements[end_col][index]
        strand = elements[strand_col][index]

        chr_seq = chrs.get(elements[chr_col][index])
        element_seq = str(chr_seq.seq)[start:end]

        if strand == '-':
            element_seq = reverse_complement(element_seq)
        elements[new_col][index] = element_seq

    return elements


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
        result.append(SeqRecord(Seq(seq, IUPAC.ambiguous_dna),
                        name=name, id=name, description=desc))
    return result


def write_fasta(reads, path_to_fasta):
    """ Function: write_fasta """
    return SeqIO.write(reads, path_to_fasta, "fasta")


def main(gff_file, genome_file, output_file):
    gff = read_gff(gff_file)
    # Extract transposons
    transposons = gff.get(gff["Feature"] == "transposable_element")

    # Helitron families ATREP[\d+] or HELITRON[\d|\w]
    helitrons = transposons[transposons["Attribute"].str.contains(r'ATREP|HELITRON')]  # noqa

    helitrons["Family"] = helitrons["Attribute"].apply(lambda x: x.split(";")[-1])
    helitrons["Family"] = helitrons["Family"].apply(lambda x: x.split("=")[-1])

    helitrons["ID"] = helitrons["Attribute"].apply(lambda x: x.split(";")[0])
    helitrons["ID"] = helitrons["ID"].apply(lambda x: x.split("=")[-1])

    del helitrons["Score"], helitrons["Frame"]
    del helitrons["Feature"], helitrons["Source"], helitrons["Attribute"]

    # Sort columns and reindex
    helitrons = helitrons[['Chr', 'ID', 'Strand', 'Start', 'End', 'Family']]
    helitrons.index = list(range(len(helitrons)))

    # Get helitron sequences
    helitrons["Chr"] = helitrons["Chr"].apply(lambda x: x.replace("Chr", ""))
    helitrons = get_sequence_from_chromosomes(helitrons, genome_file)

    # Write file
    sequences = df2reads(helitrons, id_col="ID", seq_col="Seq")
    write_fasta(sequences, output_file)


if __name__ == '__main__':
    args = argparse.ArgumentParser()

    args.add_argument("-gff", "--gff", help="GFF file", required=True)
    args.add_argument("-genome", "--genome" , required=True,
                        help="All chromososmes file (genome file)")
    args.add_argument("-out", "--output", help="path to output file",
                        required=True)
    p = args.parse_args()

    main(p.gff, p.genome, p.output)
