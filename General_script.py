#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import sys
import warnings
import pandas as pd
import seaborn as sns
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt

from os import path

if os.name == "posix":
    wdir = path.join(path.sep, "home", os.getenv("USER"), "workspace",
                     "helitrons")
elif os.name == "nt":
    wdir = path.join("D:" + path.sep, "repos", "helitrons")

sys.path.append(path.join(wdir, "scripts", "helitrons"))

from basic_functions import read_gff  # noqa
from basic_functions import run_blastx  # noqa
from basic_functions import reads2df, read_fasta  # noqa
from basic_functions import df2reads, write_fasta  # noqa
from basic_functions import get_helend_coordinates  # noqa
from basic_functions import get_sequence_from_chromosomes  # noqa
from basic_functions import read_blast_6, group_captures_perl  # noqa
from basic_functions import run_blastn, search_helend_structure  # noqa

# #############################################################################
# Workspace configuration
warnings.filterwarnings("ignore")

# Pandas configuration
pd.set_option("display.width", 200)

# Matplotlib config
plt.style.use("ggplot")
plt.set_cmap = plt.get_cmap("coolwarm")

# Seaborn config
sns.set(context="paper", font="monospace", rc={"figure.figsize": (16, 12)})
sns.set_palette("husl")
sns.set_style("darkgrid")

# #############################################################################
# Global variables
__author__ = "Andres Aguilar"
__date__ = "14/Sep/2017"
__version__ = "0.0.1"
__mail__ = "andresyoshimar@gmail.com"

_gff_columns = ["Chr", "Source", "Feature", "Start", "End", "Score",
                "Strand", "Frame", "Attribute"]
_blast_columns = ["query_id", "target_id", "identity", "Len", "mismatch",
                  "gap_open", "q_start", "q_end", "t_start", "t_end",
                  "evalue", "bit_score", "Seq"]
# Helendout.txt columns
helend_out = ["Name", "Strand", "Start_helend", "Helend", "Hairpin",
              "Hairpin_len", "Gap_pos", "Mismatches"]

# Helend final file columns
helend_final = ["Name", "Strand", "Start_helend", "Helend", "Hairpin",
                "Hairpin_len", "Gap_pos", "Mismatches", "Hairpin_start",
                "Hairpin_end", "CTRR_start", "CTRR_end"]


# #############################################################################
# Functions
def gc_percent(x):
    """ Calculate %GC """
    stem1, loop, stem2 = x.split("*")
    return GC(stem1)


def split_helitron_list(x):
    """ Function: split_helitron_list """
    tmp = x.split("-")
    tmp = list([y.split(".")[0] for y in tmp])
    return tmp


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


def write_ids(ids, output_file):
    """ """
    with open(output_file, "w") as out:
        for i in ids:
            out.write(str(i) + "\n")


# #############################################################################
# Configuration
plt.style.use("ggplot")
warnings.filterwarnings("ignore")
pd.set_option("display.width", 150)

# #############################################################################
# Files

araport = {
    "annotation": path.join(wdir, "araport",
                            "Araport11_GFF3_genes_transposons.201606.gff"),
    "genome": path.join(wdir, "araport", "TAIR10_chr_all.fas")
}

# #############################################################################
# Analysis

gff = read_gff(araport["annotation"])

# Extect genes sequences
genes = gff.get(gff["Feature"] == "gene")
genes["Chr"] = genes["Chr"].apply(lambda x: x.replace("Chr", ""))
genes = genes.get(genes["Chr"].isin(["1", "2", "3", "4", "5"]))
genes = get_sequence_from_chromosomes(genes, araport["genome"])
genes["ID"] = genes["Attribute"].apply(lambda x: x.split(";")[0])
genes["ID"] = genes["ID"].apply(lambda x: x.split("=")[-1])

genes_seq = df2reads(genes, id_col="ID", seq_col="Seq",
                     extra_data='Strand')

write_fasta(genes_seq, path.join(wdir, "araport", "Genes_araport.fna"))

# Extract transposons sequences
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

# Number of annotated helitrons
helitrons.shape[0]
# Out[1]: 13606

# Families
helitrons["Family"].unique()
# Out[3]:
# array(['ATREP4', 'ATREP3', 'HELITRONY1D', 'ATREP5', 'ATREP11',
#        'HELITRONY1E', 'HELITRONY3', 'ATREP15', 'ATREP10D', 'HELITRON4',
#        'ATREP14', 'ATREP2', 'HELITRONY1A', 'ATREP1', 'ATREP10B',
#        'HELITRONY2', 'ATREP9', 'ATREP2A', 'ATREP19', 'ATREP11A', 'ATREP12',
#        'ATREP7', 'ATREP6', 'ATREP16', 'HELITRONY1C', 'HELITRONY1B',
#        'HELITRONY3A', 'ATREP13', 'HELITRON2', 'ATREP10A', 'ATREP8',
#        'ATREP10C', 'ATREP18', 'HELITRON1', 'ATREP10', 'HELITRON3',
#        'HELITRON5', 'ATREP17'], dtype=object)

helitrons["Family"].unique().size
# Out[4]: 38

# Helitrons per family
helitrons["Family"].value_counts()
# Out[109]:
# ATREP3         1439
# HELITRONY3     1399
# ATREP10D       1295
# ATREP15        1003
# ATREP11         841
# ATREP4          832
# HELITRONY1D     756
# ATREP5          624
# ATREP1          498
# ATREP10B        491
# HELITRONY1E     447
# HELITRONY1B     414
# ATREP18         391
# ATREP10A        279
# HELITRONY1A     271
# ATREP9          201
# HELITRON4       199
# HELITRONY2      198
# HELITRON2       194
# ATREP19         189
# ATREP6          181
# ATREP7          164
# ATREP2          164
# HELITRON1       130
# ATREP10C        123
# HELITRONY1C     120
# ATREP2A         116
# ATREP13         106
# ATREP11A         95
# ATREP8           75
# ATREP10          63
# HELITRON5        56
# HELITRON3        53
# HELITRONY3A      49
# ATREP14          47
# ATREP16          47
# ATREP17          34
# ATREP12          22
# Name: Family, dtype: int64

# Get helitron sequences
helitrons["Chr"] = helitrons["Chr"].apply(lambda x: x.replace("Chr", ""))
helitrons = get_sequence_from_chromosomes(helitrons, araport["genome"])
# Wall time: 3min 40s i7 6g
# Wall time: 8min 12s i3 3g

sequences = df2reads(helitrons, id_col="ID", seq_col="Seq")
write_fasta(sequences, path.join(wdir, "araport", "araport_helitrons.fna"))

# #############################################################################
# Search cannonic helitrons
# Search helitrons with a conserved helend domain
#
# * There are no strand indicator in the fasta header (strand='n')

helends_file = path.join(wdir, "araport", "araport_helitrons_with_helend.tsv")
search_helend_structure(path.join(wdir, "araport", "araport_helitrons.fna"),
                        helends_file, strand='n')
 
# Araport helitrons 13606
# Wall time: 1h 21min 4s  i7 6g (WSL)
# Wall time: 1h 49min 19s i7 6g

# TAIR10 helitrons 12 945
# Wall time: 1h 32min 51s i7 6g
# Wall time: 2h 26min 9s  i3 3g

# Read table with helitrons structural information
helends = pd.read_table(helends_file, header=None, names=helend_out)

helends["GC"] = helends["Hairpin"].apply(gc_percent)
helends = helends.get(helends["GC"] >= 30)

helends.to_csv(helends_file, sep='\t', index=False, header=False)

# Get helend coordinates
get_helend_coordinates(helends_file, path.join(wdir, "araport",
                                               "araport_helitrons.fna"))

helends["Name"].unique().size  # Helitrons with a stable hairpin
# Out[52]: 4553

cannonic_helitrons_ids = helends["Name"].unique()

# Number of cannonic helitrons
cannonic_helitrons_ids.size
# Out:[3]: 4553

# Cannonic helitrons
cannonic_helitrons = helitrons.get(helitrons["ID"].isin(cannonic_helitrons_ids))  # noqa

cannonic_helitrons.to_csv(path.join(wdir, "araport",
                                    "araport_cannonic_helitrons.tsv"),
                          index=False, sep="\t")

sequences = df2reads(cannonic_helitrons, id_col="ID", seq_col="Seq",
                     extra_data="Strand,Start,End,Family")
write_fasta(sequences, path.join(wdir, "araport",
                                 "araport_cannonic_helitrons.fna"))

# Cannonic helitrons proportioin by family
# Number of families
helitrons["Family"].unique().size
# Out[3]: 38

cannonic_helitrons["Family"].unique().size
# Out[4]: 38

# Helitrons per family
helitrons["Family"].value_counts()
# Out[66]:
# ATREP3         1439
# HELITRONY3     1399
# ATREP10D       1295
# ATREP15        1003
# ATREP11         841
# ATREP4          832
# HELITRONY1D     756
# ATREP5          624
# ATREP1          498
# ATREP10B        491
# HELITRONY1E     447
# HELITRONY1B     414
# ATREP18         391
# ATREP10A        279
# HELITRONY1A     271
# ATREP9          201
# HELITRON4       199
# HELITRONY2      198
# HELITRON2       194
# ATREP19         189
# ATREP6          181
# ATREP2          164
# ATREP7          164
# HELITRON1       130
# ATREP10C        123
# HELITRONY1C     120
# ATREP2A         116
# ATREP13         106
# ATREP11A         95
# ATREP8           75
# ATREP10          63
# HELITRON5        56
# HELITRON3        53
# HELITRONY3A      49
# ATREP16          47
# ATREP14          47
# ATREP17          34
# ATREP12          22

cannonic_helitrons["Family"].value_counts()
# Out[67]:
# ATREP3         507
# ATREP10D       501
# HELITRONY3     422
# ATREP15        334
# ATREP11        273
# ATREP4         238
# ATREP18        185
# ATREP5         181
# ATREP1         179
# HELITRONY1D    172
# ATREP10B       141
# HELITRONY1E    134
# HELITRONY1B    121
# ATREP2         100
# ATREP10A        90
# HELITRON4       79
# HELITRON2       69
# HELITRONY1A     68
# ATREP2A         64
# HELITRON1       62
# ATREP9          61
# ATREP6          59
# HELITRONY2      57
# ATREP7          53
# ATREP13         51
# ATREP19         49
# ATREP8          44
# HELITRONY1C     40
# HELITRON5       32
# ATREP10C        30
# HELITRON3       29
# ATREP11A        28
# ATREP10         23
# ATREP14         20
# ATREP16         18
# ATREP17         16
# HELITRONY3A     15
# ATREP12          8

# Proportion of cannonic helitrons by family
families = pd.DataFrame()
families["helitrons"] = helitrons["Family"].value_counts()
families["cannonic_helitrons"] = cannonic_helitrons["Family"].value_counts()

families["proportion"] = families.apply(lambda x: (x["cannonic_helitrons"]*100.0)/x["helitrons"], axis=1)  # noqa

families
# Out[70]:
#              helitrons  cannonic_helitrons  proportion
# ATREP3            1439                 507   35.232801
# HELITRONY3        1399                 422   30.164403
# ATREP10D          1295                 501   38.687259
# ATREP15           1003                 334   33.300100
# ATREP11            841                 273   32.461356
# ATREP4             832                 238   28.605769
# HELITRONY1D        756                 172   22.751323
# ATREP5             624                 181   29.006410
# ATREP1             498                 179   35.943775
# ATREP10B           491                 141   28.716904
# HELITRONY1E        447                 134   29.977629
# HELITRONY1B        414                 121   29.227053
# ATREP18            391                 185   47.314578
# ATREP10A           279                  90   32.258065
# HELITRONY1A        271                  68   25.092251
# ATREP9             201                  61   30.348259
# HELITRON4          199                  79   39.698492
# HELITRONY2         198                  57   28.787879
# HELITRON2          194                  69   35.567010
# ATREP19            189                  49   25.925926
# ATREP6             181                  59   32.596685
# ATREP2             164                 100   60.975610
# ATREP7             164                  53   32.317073
# HELITRON1          130                  62   47.692308
# ATREP10C           123                  30   24.390244
# HELITRONY1C        120                  40   33.333333
# ATREP2A            116                  64   55.172414
# ATREP13            106                  51   48.113208
# ATREP11A            95                  28   29.473684
# ATREP8              75                  44   58.666667
# ATREP10             63                  23   36.507937
# HELITRON5           56                  32   57.142857
# HELITRON3           53                  29   54.716981
# HELITRONY3A         49                  15   30.612245
# ATREP16             47                  18   38.297872
# ATREP14             47                  20   42.553191
# ATREP17             34                  16   47.058824
# ATREP12             22                   8   36.363636

families.sort_values("proportion", ascending=False, inplace=True)

# Generate bar plot
families.head(20).plot.bar()
plt.ylim(0, 500)
plt.title("Cannonic helitrons per family")
plt.ylabel("Number of cannonic helitrons")
plt.xlabel("Families")
plt.yticks(pd.np.arange(0, 500 , 50))
plt.tight_layout()
plt.savefig(path.join(wdir, "araport", "graphics", "cannonic_helitrons_per_family.png"),
            dpi=400, format="png")
plt.close()

# #############################################################################
# Cut sequences from 5' aTC to 3' Helend domain (CTRR)
# If one helitron have more than one helend domain select most 3' proximal
helends_file = path.join(wdir, "araport",
                         "araport_helitrons_with_helend_coordinates.tsv")

df = pd.read_table(helends_file, header=None, names=helend_final)

df["GC"] = df["Hairpin"].apply(gc_percent)
df = df.get(df["GC"] >= 30)

df["Name"].unique().size
# Out[32]: 4553  Match!

del df["GC"]

# Sort values highest CTRR coord first
df.sort_values(["Name", "CTRR_end"], inplace=True, ascending=[True, False])
df.drop_duplicates("Name", inplace=True)

df.to_csv(path.join(wdir, "araport", "helends_unique.tsv"), index=False,
                    header=False, sep="\t")

cannonic_hels = cut_cannonical_sequence(path.join(wdir, "araport",
                                                  "helends_unique.tsv"),
                                        path.join(wdir, "araport",
                                                  "araport_helitrons.fna"))

cannonical_sequences = path.join(wdir, "araport",
                                 "helitrons_with_cannonical_sequences.tsv")
with open(cannonical_sequences, "w") as out:
    out.write("Name\tSequence\tStart\tEnd\n")
    for key in cannonic_hels.keys():
        element = cannonic_hels.get(key)
        element = [str(x) for x in element]
        out.write(key + "\t" + "\t".join(element) + "\n")

# ###################
# Cut first TC motif

df = reads2df(read_fasta(path.join(wdir, "araport",
                                   "araport_cannonic_helitrons.fna")))

df["TC_Seq"] = df["Seq"].apply(lambda x: x[x.find("TC"):])
write_fasta(df2reads(df, seq_col="TC_Seq"),
            path.join(wdir, "araport", "tc_helitrons.fna"))
# Search helend domain
helends_file = path.join(wdir, "araport", "tc_helitrons_with_helend.tsv")
search_helend_structure(path.join(wdir, "araport", "tc_helitrons.fna"),
                        helends_file, strand='n')

helends = pd.read_table(helends_file, header=None, names=helend_out)

helends["GC"] = helends["Hairpin"].apply(gc_percent)
helends = helends.get(helends["GC"] >= 30)
del helends["GC"]
helends.to_csv(helends_file, sep='\t', index=False, header=False)

# Get helend coordinates
get_helend_coordinates(helends_file, path.join(wdir, "araport",
                                               "tc_helitrons.fna"))

helends_file = path.join(wdir, "araport",
                         "tc_helitrons_with_helend_coordinates.tsv")
helends = pd.read_table(helends_file, header=None, names=helend_final)
helends["Name"].unique().size  # Helitrons with a stable hairpin

# Sort values highest CTRR coord first
helends.sort_values(["Name", "CTRR_end"], inplace=True,
                    ascending=[True, False])
helends.drop_duplicates("Name", inplace=True)
helends.to_csv(helends_file, sep="\t", index=False, header=False)

hels = cut_cannonical_sequence(helends_file, path.join(wdir, "araport",
                                                       "tc_helitrons.fna"))

# Write Fasta file with cannonical sequence
cannonical_sequences = path.join(wdir, "araport",
                                "helitrons_cannonical_sequence.fna")
with open(cannonical_sequences, "w") as out:
    for key in hels.keys():
        element = hels.get(key)
        out.write(">{0}\n{1}\n".format(key, element[0]))


# #############################################################################
# Find captures
df = reads2df(read_fasta(path.join(wdir, "araport",
                                   "araport_cannonic_helitrons.fna")))


def add_length(x):
    return x.Name + "|RC/Helitron|" + str(len(x.Seq))


df["Id"] = df.apply(add_length, axis=1)
write_fasta(df2reads(df, id_col="Id"),
            path.join(wdir, "araport",
                      "araport_cannonic_helitrons.fna"))

db = path.join(wdir, "..", "ref_files", "blast_db", "Ath_araport")
out_file = path.join(wdir, "araport", "helitrons_genes.blast6")

run_blastn(db, path.join(wdir, "araport", "helitrons_cannonical_sequence.fna"),
           out_file, wseq=True)

# Extract alignments with identity == 100
df = read_blast_6(path.join(wdir, "araport", "helitrons_genes.blast6"),
                  wseq=True)

df = df.get(df["identity"] == 100)
df.shape
# Out[52]: (3510, 13)

# Minimum length to consider an alignment as a capture
df = df.get(df["Len"] >= 50)
df.shape
# Out[55]: (906, 13)

df["query_id"].value_counts().head()
# Out[61]: 
# AT1TE41625|-|12768947|12780004|ATREP15|RC/Helitron|11058    23
# AT3TE69475|+|17143391|17157692|ATREP3|RC/Helitron|14302     17
# AT3TE70785|+|17469788|17471246|ATREP12|RC/Helitron|1459      7
# AT5TE84250|+|23431544|23432883|ATREP12|RC/Helitron|1340      7
# AT5TE54260|+|15032076|15032771|ATREP15|RC/Helitron|696       6

# Sort alignments decreasing mode
df.sort_values("Len", inplace=True, ascending=False)

df.to_csv(path.join(wdir, "araport", "alignments_helitrons_genes_100.tsv"),
          sep='\t', index=False, header=False)

# Group captures
captures_file = path.join(wdir, "araport", "cannonic_captures.tsv")
blast_file = path.join(wdir, "araport",  "cannonic_bout.tsv")

group_captures_perl(out_file, captures_file, blast_file)

# read captures file
captures = pd.read_table(captures_file)
captures.sort_values("NoHelitrons", ascending=False, inplace=True)
captures.index = list(range(len(captures)))  # reindex

captures.to_csv(captures_file, index=False, sep='\t')

# Captured genes
captures["Gene"].unique().size
# Out[52]: 1719
# By Cannonic captures: 1,285

# Cannonic helitrons with capture
captures["Helitron_list"] = captures["Helitrons"].apply(split_helitron_list)

helitrons_list = list()
for index in captures.index:
    tmp_hels = list(set(captures["Helitron_list"][index]))
    helitrons_list.extend(tmp_hels)

helitrons_list = list(set(helitrons_list))
len(helitrons_list)  # Number of cannonic helitrons with capture
# Out[14]: 3140
# Cannonic helitrons: 2,603

# #############################################################################
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
genes_dist_file = path.join(wdir, "araport",
                            "cannonic_captured_genes_dist.tsv")
captured_genes_dist.to_csv(genes_dist_file, sep='\t', index=False)

# Total number of captures
captured_genes_dist["NoHelitrons"].sum()
# Out[50]: 7204
# Cannonic: 5,346

# Plotting
captured_genes_dist.index = captured_genes_dist["Gene"].tolist()
captured_genes_dist.head(50).plot.bar()
plt.title("Most captured genes (first 50)")
plt.ylabel("Number of captures")
plt.xlabel("Genes")
plt.tight_layout()
mcg_file = path.join(wdir, "araport", "graphics",
                     "cannonic_captured_genes_dist.png")
plt.savefig(mcg_file, format="png")
plt.close()

# #############################################################################
# Annonate cannonical sequence
df = read_blast_6("annotation/cannonical_sequence.bast6")
df = df.get(df["identity"] == 100)
df = df.get(df["q_start"] == 1)

helitrons = reads2df(read_fasta(path.join(wdir, "araport",
                                          "helitrons_cannonical_sequence.fna")))
helitrons["Len"] = helitrons["Seq"].apply(lambda x: len(x))

df = df.get(df["Len"] >= helitrons["Len"].min())

df.index = df["query_id"].tolist()
helitrons.index = helitrons["Name"].tolist()

df["Helitron_len"] = helitrons["Len"]
df.index = list(range(len(df)))

df["Name"] = df["query_id"].apply(lambda x: x.split("|")[0])
df["Chr"] = df["Name"].apply(get_chromosome)
df["Chr"] = df["Chr"].apply(lambda x: x.replace("Chr", ""))

df["Same_chr"] = df.apply(lambda x: True if x.target_id == x.Chr else False, axis=1)
df = df.get(df["Same_chr"] == True)
del df["Same_chr"]

df["Helitron_start"] = df["query_id"].apply(lambda x: x.split("|")[2])
df["Helitron_end"] = df["query_id"].apply(lambda x: x.split("|")[3])

df["Helitron_start"] = df["Helitron_start"].astype(pd.np.int64)
df["Helitron_end"] = df["Helitron_end"].astype(pd.np.int64)


def check_coords(x):
    start = x.Helitron_start <= x.t_start
    end = x.Helitron_end >= x.t_end

    if start and end:
        return True
    else:
        return False

df["Check"] = df.apply(check_coords, axis=1)

cannonic = df.get(df["Check"] == True)
cannonic = cannonic[["Name", "query_id", "t_start", "t_end", "Helitron_len"]]


def modify_name(x):
    name = x.query_id.split("|")
    # 0:Name,1:Strand,2:Start,3:End:4Family,5:Type,6:Len
    name[2] = str(x.t_start)
    name[3] = str(x.t_end)
    name[6] = str(x.Helitron_len)
    return "|".join(name)


cannonic["Long_name"] = cannonic.apply(modify_name, axis=1)
cannonic[["Long_name"]].to_csv(path.join(wdir, "araport", "annotation",
                                         "header_with_correct_annotation.tsv"),
                               sep="\t", index=False)
# Completed by hand
df = pd.read_table(path.join(wdir, "araport", "annotation",
                                         "header_with_correct_annotation.tsv"))

df["Name"] = df["Long_name"].apply(lambda x: x.split("|")[0])

df.index = df["Name"].tolist()
helitrons.index = helitrons["Hel_name"].tolist()

helitrons["Name"] = df["Long_name"]

reads = df2reads(helitrons)
write_fasta(reads, path.join(wdir, "araport", "annotation",
                            "cannonic_helitrons.fna"))
# #############################################################################
# Classify captures
genes = captures["Gene"].unique().tolist()

hels = list()
gn_list = list()
cap_range = list()
class_list = list()
for gene in genes:
    tmp = captures.get(captures["Gene"] == gene)
    tmp_hels = list()
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

df.shape[0]
# Out[78]: 7204
# Cannonic: 5,346

df["Chr_hel"].value_counts()
# Out[145]:
# Chr1    1904
# Chr5    1575
# Chr2    1329
# Chr3    1255
# Chr4    1141

df["Chr_gene"].value_counts()
# Out[146]:
# Chr4    1915
# Chr1    1660
# Chr5    1322
# Chr3    1168
# Chr2    1139

df["Same_chr"] = df.apply(lambda x: x.Chr_gene == x.Chr_hel, axis=1)
df["Same_chr"].value_counts()
# Out[148]:
# False    5094
# True     2110

# Get genomic coordinates for helitrons and genes
annotation = read_gff(araport["annotation"])

annotation = annotation[annotation["Feature"].str.contains(r"^gene$")]  # noqa
annotation["Feature"].unique()
# Out[151]: array(['gene', 'transposable_element'], dtype=object)

annotation["Id"] = annotation["Attribute"].apply(lambda x: x.split(";")[0])
annotation["Id"] = annotation["Id"].apply(lambda x: x.split("=")[-1])

df["Gene_start"] = 0
df["Gene_end"] = 0
for x in df.itertuples():
    g_tmp = annotation.get(annotation["Id"] == x.Gene)
    df["Gene_start"][x.Index] = g_tmp["Start"][g_tmp.first_valid_index()]
    df["Gene_end"][x.Index] = g_tmp["End"][g_tmp.first_valid_index()]

    # h_tmp = annotation.get(annotation["Id"] == x.Helitron)
    # df["Hel_start"][x.Index] = h_tmp["Start"][h_tmp.first_valid_index()]
    # df["Hel_end"][x.Index] = h_tmp["End"][h_tmp.first_valid_index()]

helitrons["Start"] = helitrons["Name"].apply(lambda x: int(x.split("|")[2]))
helitrons["End"] = helitrons["Name"].apply(lambda x: int(x.split("|")[3]))

df.index = df["Helitron"].tolist()

df["Hel_start"] = helitrons["Start"]
df["Hel_end"] = helitrons["End"]

df.index = list(range(len(df)))

df.shape[0]
# Out[207]: 7204
# Cannonic: 5,346

df.to_csv(path.join(wdir, "araport", "captures_extended.tsv"), sep="\t",
          index=False)

# ##############################################################################
# Clasificación de capturas genicas
#
# Clase I - Captura génica
#   a) Cualquier caso en el cual un gen completo se encuentre dentro de un
#      helitron
#
#   b) Cualquier caso en el cual un fragmento de un gen se encuentre fuera de
#      las coordenadas genómicas del helitron, ej. en otro cromosoma
#
# Clase II - Helitrones dentro de genes
#      Cualquier helitron que se encuentre completamente dentro de las
#      coordenadas de un gen
#
# Clase III - Intersección de helitrones con genes
#      Cualquier caso en el cual una porción de un helitron se traslape con la
#      de un gen sin estar completamente contenido
# #############################################################################

df["Class"] = ""
for x in df.itertuples():
        starts = x.Gene_start >= x.Hel_start and x.Gene_start <= x.Hel_end
        ends = x.Gene_end >= x.Hel_start and x.Gene_end <= x.Hel_end

        # Revisamos si cumple con los requisitos de la clase Ia
        if x.Gene_start >= x.Hel_start and x.Gene_end <= x.Hel_end:
            # El gen está contenido en el helitron
            df["Class"][x.Index] = "Ia"
        elif x.Hel_start >= x.Gene_start and x.Hel_end <= x.Gene_end:
            # El helitron está contenido en el gen
            df["Class"][x.Index] = "II"
        elif starts or ends:
            # las coordenadas se sobrelapan en algún punto pero ningún elemento
            # contiene al otro
            df["Class"][x.Index] = "III"
        else:
            # el gen y el helitron están en el mismo cromosoma pero sus
            # coordenadas no sobrelapan en ningún punto
            df["Class"][x.Index] = "Ib"

df.shape[0]
# Out[210]: 7204
# Cannonic: 5,346

df["Class"].value_counts()
# Out[159]:
# Ib     6454
# III     438
# II      231
# Ia       81

# Cannonic
# Ib     4825
# II      378
# III     114
# Ia       29

classification_file = path.join(wdir, "araport",
                                "cannonic_classification_table2.tsv")
df.to_csv(classification_file, index=False, sep='\t')

ia = df.get(df["Class"] == "Ia")
ib = df.get(df["Class"] == "Ib")
ii = df.get(df["Class"] == "II")
iii = df.get(df["Class"] == "III")

# #############################################################################
# Gene Onthology enrichment
# http://www.webgestalt.org/option.php
write_ids(df["Gene"].unique(),
          path.join(wdir, "araport", "GO", "cannonic_all_ids.txt"))
write_ids(ia["Gene"].unique(),
          path.join(wdir, "araport", "GO", "cannonic_class_Ia.txt"))
write_ids(ib["Gene"].unique(),
          path.join(wdir, "araport", "GO", "cannonic_class_Ib.txt"))
write_ids(ii["Gene"].unique(),
          path.join(wdir, "araport", "GO", "cannonic_class_II.txt"))
write_ids(iii["Gene"].unique(),
          path.join(wdir, "araport", "GO", "cannonic_class_III.txt"))


cols = ["geneset", "description", "OverlapGene_UserID"]
df = pd.read_table(path.join(wdir, "araport", "GO", "cannonic",
                             "enrichment_results.tsv"))
df = df[cols]

df.columns = ["GOId", "Category", "Genes"]
df["Count"] = df["Genes"].apply(lambda x: len(x.split(";")))

df
# Out[7]:
#          GOId                                           Category                                              Genes  Count
# 0  GO:0001906                                       cell killing  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     23
# 1  GO:0031640                 killing of cells of other organism  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     23
# 2  GO:0044364              disruption of cells of other organism  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     23
# 3  GO:0035821  modification of morphology or physiology of ot...  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     23
# 4  GO:0050832                         defense response to fungus  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     24

df = df.drop_duplicates("Genes")
df
# Out[9]:
#          GOId                            Category                                              Genes  Count
# 0  GO:0001906                        cell killing  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     23
# 4  GO:0050832          defense response to fungus  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     24
# 6  GO:0051704              multi-organism process  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     24
# 7  GO:0016042             lipid catabolic process  AT2G03980;AT2G19010;AT2G19060;AT3G14820;AT3G43...     12
# 8  GO:0006952                    defense response  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     35
# 9  GO:0098542  defense response to other organism  AT2G24693;AT5G55565;AT1G13605;AT1G59833;AT2G22...     25

df["Genes"] = df["Genes"].apply(lambda x: x.split(";"))
df.index = list(range(len(df)))

df.to_csv(path.join(wdir, "araport", "GO", "cannonic", "enrichment_results_filtered.tsv"),
          index=False, sep='\t')
# According to the Venn diagram
# Cell Killing category ids are the 'core' of all other categories

df["Category"][0]  # index = 0
# Out[186]: 'cell killing'

df["Genes"][0]
# Out[14]:
# ['AT2G24693', <- Encodes a defensin-like (DEFL) family protein. (Protein of unknown function (PD694200))
#  'AT5G55565', <- Encodes a defensin-like (DEFL) family protein.
#  'AT1G13605', <- Encodes a defensin-like (DEFL) family protein.
#  'AT1G59833', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G22807', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G13542', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G03932', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G03933', <- Encodes a defensin-like (DEFL) family protein. (Cysteine-rich protein)
#  'AT2G24625', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G04045', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G24615', <- Encodes a defensin-like (DEFL) family protein.
#  'AT2G22805', <- Encodes a defensin-like (DEFL) family protein.
#  'AT3G17155', <- Encodes a defensin-like (DEFL) family protein. (Cysteine-rich protein)
#  'AT3G27831', <- Encodes a defensin-like (DEFL) family protein. (Gamma-thionin family protein)
#  'AT3G45093', <- Encodes a defensin-like (DEFL) family protein.
#  'AT4G08869', <- Encodes a defensin-like (DEFL) family protein. (Putative membrane lipoprotein)
#  'AT5G19175', <- Encodes a defensin-like (DEFL) family protein. (Cysteine-rich protein)
#  'AT5G32619', <- Encodes a defensin-like (DEFL) family protein.
#  'AT5G37474', <- Encodes a defensin-like (DEFL) family protein. (Putative membrane lipoprotein)
#  'AT5G43401', <- Encodes a defensin-like (DEFL) family protein.
#  'AT5G62627', <- Encodes a defensin-like (DEFL) family protein. (Putative membrane lipoprotein)
#  'AT2G20465', <- Encodes a defensin-like (DEFL) family protein. (Molecular chaperone Hsp40/DnaJ family protein)
#  'AT3G04540'] <- Encodes a defensin-like (DEFL) family protein. (Cysteine-rich protein)
# Source: araport.org

write_ids(df["Genes"][0],
          path.join(wdir, "araport", "GO", "cannonic", "Core_ids.txt"))

# #############################################################################
# Duplicates gene
# Helitron    AT3TE62100 (HELITRON2)
# Gene        AT3G43340  (Pseudouridine synthase family protein)
# Duplicated  AT2G39140  (pseudouridine synthase family protein)
# Gooooooood!

# #############################################################################
# Thionins
# Helitron   AT1TE41625  (ATREP15)
# Gene1      AT1G34795   (Plant thionin family protein)
# Gene2      AT1G34800   (Plant thionin family protein)
# Gene3      AT1G34805   (Plant thionin family protein)

df = read_blast_6(path.join(wdir, "araport", "helitrons_genes.blast6"),
                  wseq=True)

df = df.get(df["query_id"].str.contains("AT1TE41625"))
df.to_csv(path.join(wdir, "araport", "antimicrobial_peptides",
                    "blast_out6_thionins.tsv"),
          index=False, header=None, sep='\t')

thionins = ["AT1G34795|-", "AT1G34800|-", "AT1G34805|-"]
df = df.get(df["target_id"].isin(thionins))
df.to_csv(path.join(wdir, "araport", "antimicrobial_peptides",
                    "blast_out6_only_thionins.tsv"),
          index=False, header=None, sep='\t')

captures_file = path.join(wdir, "araport", "antimicrobial_peptides", "captures.tsv")
blast_file = path.join(wdir, "araport", "antimicrobial_peptides", "bout.tsv")

group_captures_perl(path.join(wdir, "araport", "antimicrobial_peptides",
                    "blast_out6_only_thionins.tsv"), captures_file, blast_file)

# #############################################################################
# Search RepHel proteins

df = reads2df(read_fasta(path.join(wdir, "..", "ref_files",
                                   "Araport11_genes.201606.pep.fasta")))
df = df[df["Name"].str.contains(r"helicase|replicationprotein|replicase")]

rephel_proteins_fasta = path.join(wdir, "araport", "autonomous",
                                  "rephel_proteins.fna")

write_fasta(df2reads(df), rephel_proteins_fasta)

# Make blast database

db = path.join(wdir, "..", "ref_files", "rephel", "rephel")
query = path.join(wdir, "araport", "annotation", "cannonic_helitrons.fna")
rephel_blastx = path.join(wdir, "araport", "autonomous",
                                  "rephel_helitrons.blastx")

# evalue 10x10-5
blastx = run_blastx(db=db, query=query, out_file=rephel_blastx, wseq=True)

df_blastx = read_blast_6(blastx, wseq=True)
df_blastx.shape
# Out[97]: (298, 13)

df_blastx.describe()
# Out[100]:
#          identity         Len    mismatch    gap_open       q_start         q_end     t_start       t_end         evalue   bit_score
# count  298.000000  298.000000  298.000000  298.000000    298.000000    298.000000  298.000000  298.000000   2.980000e+02  298.000000
# mean    49.105604  112.741611   53.436242    1.308725   4946.694631   5281.352349  111.765101  218.067114   2.501136e-06   99.785906
# std     12.447537   71.617698   36.170904    1.449018   5024.623990   5075.272690   66.307151   99.497198   1.143704e-05   58.889796
# min     23.530000   18.000000    0.000000    0.000000      1.000000     60.000000    1.000000   40.000000  4.000000e-158   25.400000
# 25%     41.827500   53.250000   23.000000    0.000000    504.500000    762.250000   58.000000  141.000000   3.250000e-52   49.100000
# 50%     46.755000   89.000000   48.000000    1.000000   3076.000000   3569.500000  116.500000  224.500000   1.500000e-34   81.800000
# 75%     52.380000  165.500000   78.750000    2.000000   9324.750000   9586.250000  144.000000  327.750000   3.750000e-17  147.000000
# max    100.000000  324.000000  163.000000    8.000000  15933.000000  16184.000000  315.000000  344.000000   8.000000e-05  417.000000

df_blastx["query_id"].unique().size
# Out[103]: 75

helicase = df_blastx[df_blastx["target_id"].str.contains("helicase")]
replicase = df_blastx[df_blastx["target_id"].str.contains(r"replication|replicase")]

# Helitrons with helicase motif (with duplicates)
helicase.shape[0]
# Out[107]: 244

# Helitrons with replicase motif (with duplicates)
replicase.shape[0]
# Out[108]: 54

# Sort alignements and remove duplicates
helicase.sort_values(["Len", "mismatch", "evalue"],
                     ascending=[False, True, True], inplace=True)
replicase.sort_values(["Len", "mismatch", "evalue"],
                      ascending=[False, True, True], inplace=True)

# Remove duplicates and keep the best alignment
helicase.drop_duplicates("query_id", keep='first', inplace=True)
replicase.drop_duplicates("query_id", keep='first', inplace=True)

helicase.shape[0]
# Out[113]: 54

replicase.shape[0]
# Out[114]: 29

helicase[["identity", "Len", "evalue"]].describe()
# Out[115]:
#          identity         Len         evalue
# count   54.000000   54.000000   5.400000e+01
# mean    48.145926  161.777778   1.297000e-06
# std     13.410316   70.366069   9.525696e-06
# min     27.860000   20.000000  4.000000e-158
# 25%     40.657500  104.750000   2.037500e-56
# 50%     46.410000  165.000000   1.900000e-44
# 75%     50.857500  200.500000   2.262500e-20
# max    100.000000  324.000000   7.000000e-05

replicase[["identity", "Len", "evalue"]].describe()
# Out[116]:
#          identity         Len        evalue
# count   29.000000   29.000000  2.900000e+01
# mean    53.182759   66.758621  7.496831e-06
# std     15.832915   32.188790  1.714900e-05
# min     33.800000   27.000000  6.000000e-65
# 25%     42.370000   42.000000  4.000000e-26
# 50%     50.560000   59.000000  8.000000e-09
# 75%     53.930000   87.000000  3.000000e-06
# max    100.000000  154.000000  6.000000e-05

# Get helitrons with helicase and replicase
autonomous = helicase.get(helicase["query_id"].isin(replicase["query_id"]))["query_id"].tolist()
autonomous
# Out[119]: 
# ['AT2TE25020|+|6148550|6161787|HELITRON1|RC/Helitron|13238',
#  'AT5TE47085|-|13252607|13244535|HELITRON1|RC/Helitron|8073',
#  'AT1TE43740|+|13371475|13386480|HELITRON1|RC/Helitron|15006',
#  'AT3TE58655|-|14285452|14272958|HELITRON5|RC/Helitron|12495',
#  'AT2TE08320|-|1821216|1804018|HELITRON4|RC/Helitron|17199',
#  'AT5TE52795|+|14659366|14675167|HELITRON4|RC/Helitron|15802',
#  'AT3TE18040|-|4288576|4276853|HELITRON1|RC/Helitron|11724',
#  'AT1TE51985|+|15772990|15789174|HELITRON4|RC/Helitron|16185']

len(autonomous)
# Out[122]: 8


# #############################################################################
# Get captures by class and perform a Blast search to find if we can find it
# in another location (True capture). Only for class Ia, II and III.

# bout file format
# Gene Helitron CapSubClass Attributes
# Gene|Aln_len|Hel_coords|Gene_coords|Hel_strand|Gene_strand|Cap_direction_hel|Cap_direction_gene
bout = pd.read_table(path.join(wdir, "araport", "cannonic_bout.tsv"), header=None)
cols = ["Gene", "Helitron", "CClass", "Attributes"]
bout.columns = cols

cols2 = ["Gene", "Aln_len", "Hel_coords", "Gene_coords", "Hel_strand",
         "Gene_strand", "Cap_direction_hel", "Cap_direction_gene"]

# Start with class III
classes = pd.read_table(path.join(wdir, "araport", "cannonic_classification_table2.tsv"))
df = classes.get(classes["Class"] == "III")

df.shape[0]
# Out[18]: 114 elements

df["Cap_len"] = 0
df["Gene_coords"] = ""
for cap in df.itertuples():
    temp = bout[(bout["Gene"] == cap.Gene) & (bout["Helitron"] == cap.Helitron)]
    temp = temp.to_dict()
    k = temp["Attributes"].keys()
    k = k[0]

    val = temp["Attributes"][k]
    val = val.split("|")
    df["Cap_len"][cap.Index] = val[1]
    df["Gene_coords"][cap.Index] = val[3]

# read blast out file
blast = read_blast_6(path.join(wdir, "araport", "cannonic_helitrons_genes.blast6"), wseq=True)
blast["Helitron"] = blast["query_id"].apply(lambda x: x.split("|")[0])
blast["Gene"] = blast["target_id"].apply(lambda x: x.split("|")[0])

df["Seq_cap"] = ""
for cap in df.itertuples():
    temp = blast[(blast["Gene"] == cap.Gene) & (blast["Helitron"] == cap.Helitron) & (blast["Len"] == cap.Cap_len)]
    df["Seq_cap"][cap.Index] = temp["Seq"][temp.first_valid_index()]


df["Name"] = df.apply(lambda x: "|".join([x.Gene, x.Helitron, x.Class, str(x.Cap_len)]), axis=1)

write_fasta(df2reads(df, seq_col="Seq_cap"), path.join(wdir, "araport", "Classes", "class_III.fna"))
del df["Name"]
df.to_csv(path.join(wdir, "araport", "Classes", "class_III.tsv"), sep="\t", index=False)

# BLAST

df = read_blast_6(path.join(wdir, "araport", "Classes", "classIII.blast6"))
df["query_len"] = df["query_id"].apply(lambda x: int(x.split("|")[-1]))



temp = df.get(df["identity"] == 100)
temp.get(df["identity"] == 100)["query_id"].value_counts().head()
# Out[118]:
# AT1G17147|AT1TE18955|III|403    3 <- Self hit
# AT5G57840|AT5TE84250|III|94     2 <- same sequence 2 chromosomes (both with a helitron)
# AT5G47950|AT5TE69870|III|61     1

temp.get(temp["query_id"] == "AT1G17147|AT1TE18955|III|403")[["query_id", "target_id", "identity", "Len"]]
# Out[128]:
#                          query_id  target_id  identity  Len
# 21   AT1G17147|AT1TE18955|III|403          1     100.0  403 <- Self hit
# 104  AT1G17147|AT1TE18955|III|403          2     100.0   38 <-| 
# 144  AT1G17147|AT1TE18955|III|403          3     100.0   28 <-|

temp.get(temp["query_id"] == "AT5G57840|AT5TE84250|III|94")[["query_id", "target_id", "identity", "Len"]]
# Out[134]:
#                         query_id  target_id  identity  Len
# 227  AT5G57840|AT5TE84250|III|94          5     100.0   94 <- Self hit
# 230  AT5G57840|AT5TE84250|III|94          3     100.0   94 <- fragment of antoher helitron

# Delete self hits and duplications
temp = df.get(df["identity"] != 100)

# Delete small hits
temp = temp[temp["Len"] >= (temp["query_len"]-10)]

# Sort by identity
temp.sort_values("identity", inplace=True, ascending=False)
# Drop duplicates
temp.drop_duplicates("query_id", inplace=True)

temp = temp.get(temp["identity"] >= 90)

# In [147]: temp
# Out[147]: 
#                     query_id  target_id  identity  Len  mismatch  gap_open  q_start  q_end   t_start     t_end         evalue  bit_score 
# AT1G17147|AT1TE18955|III|403         2     98.77  405         3         1        1    403   2138980   2139384   0.000000e+00      719.0   
# AT5G04255|AT5TE33945|III|54          4     98.15   54         1         0        1     54  14979213  14979160   8.000000e-20       95.3   
# AT5G57840|AT5TE84250|III|94          2     97.75   89         1         1        7     94   8893626   8893538   1.000000e-36      152.0   
# AT3G31550|AT3TE52765|III|488         3     96.52  488        17         0        1    488  13222099  13222586   0.000000e+00      808.0   
# AT5G43518|AT5TE63200|III|311         5     95.21  313        10         4        3    311  17482199  17482510  8.000000e-138      490.0   
# AT1G57777|AT1TE70660|III|57          1     94.74   57         3         0        1     57  21400094  21400150   4.000000e-18       89.8   
# AT1G14518|AT1TE16080|III|207         2     93.72  207        11         2        2    207    795431    795636   1.000000e-83      309.0   
# AT5G32619|AT5TE43605|III|215         5     93.49  215        12         2        1    215  12419858  12420070   3.000000e-86      318.0   

# AT1G17147|AT1TE18955|III|403 2:2138980..2139384
# Other transposon AT3TE09015|+|SIMPLEGUY1|908 
# DNA/Harbinger  with putative helitron cannonic structure.
# TC <- 4,5b
# Helend <- 863b ATACTAAATAATCAAACAAACATTATTCAAATGCATCATCTAA

# AT5G04255|AT5TE33945|III|54 4:14979160..14979213
# Other cannonic helitron  AT4TE71140|-|ATREP6|1178

# AT5G57840|AT5TE84250|III|94 2:8893538..8893626
# Other cannonic helitron AT2TE37435|-|ATREP12|1702 

# AT3G31550|AT3TE52765|III|488 3:13222099..13222586
# Other cannonic helitron AT3TE54065|+|HELITRONY2|1030

# AT5G43518|AT5TE63200|III|311 5:17482199..17482510
# Other cannonic helitron AT5TE63190|+|ATREP1|1765  <--- Pseudogen with transcription

# AT1G57777|AT1TE70660|III|57 1:21400094..21400150
# 3' UTR Gene: AT1G57775|-|ECA1 gametogenesis family protein (DUF784)

# AT1G14518|AT1TE16080|III|207 2:795431..795636
# 5' UTR Gene: AT2G02800|-|protein kinase 2B

# AT5G32619|AT5TE43605|III|215 5:12419858..12420070
# Other non cannonic helitron AT5TE44020|+|ATREP3|1824

# #############################################################################
# Phylogeny with helend domains
df = pd.read_table(path.join(wdir, "araport", "tc_helitrons_with_helend.tsv"),
                   header=None)
df.columns = ["Id", "Strand", "Pos", "Helend", "Hairpin", "HP_len",
              "Gap_pos", "Gaps"]

df["Name"] = df.apply(lambda x: "|".join([str(x.Id), str(x.Pos)]), axis=1)

write_fasta(df2reads(df, seq_col="Helend"),
            path.join(wdir, "araport", "Classes", "helends.fna"))
