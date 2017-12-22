#!/usr/bin/env perl

# MODULE 3 FOR TALLO-ASA ALGORITHM
#
# Author: Pablo E. García
#
# Date: oct 19, 2012
#
# Description: Finds the coordinates of the hairpin-loop structure
# and CTRR obtained by the helend.pl program
#
# Inputs: 
#	helendout.txt: archivo de salida de helend.pl
#	original_fasta.fas: archivo de entrada de helend.pl
#
# Outputs:
#	helendout.txt: the same file but with two more columns (start, end)
#
# Example:
#	helend_coordinates.pl helendout.txt original_fasta.fas
#
#

use strict;

my $helend = shift;
my $fasta = shift;
my $id;
my %fasta;  # Stores de FASTA sequences in a hash with keys = IDs

open (FASTA, "$fasta") or die "cannot open the fasta file \n";
while (<FASTA>) {
    chomp;
    if (/>/) {
	$id = $_;
	$id =~ s/>//g;
	
	$fasta{$id} = "";
    }
    else {
	$fasta{$id}.=$_
    }
}
close (FASTA);

# Retrieves the exact coordinates of the haripin loop structure,
# stored in the 5 column of the helendout.txt file
system("cp $helend temporal");

open(HELEND, "temporal") or die "cannot open the helendout file \n";
open(OUT, ">$helend") or die "cannot write the file\n";

my @helend;
my ($start, $end, $wholesequence, $wholesequence2,  $gap_pos, $miss_number, $hairloop, $CTRR_end, $CTRR_start, $strand) = "";

while(my $row = <HELEND>) {
    chomp($row);
    
    # Extracts information from helendout.txt
    @helend = split("\t", $row);
    $id = $helend[0];
    $strand = $helend[1];
    $wholesequence = $helend [3];
    $wholesequence2 = $helend [3];
    $hairloop = $helend[4];
    $hairloop =~ s/\*//g;
    $gap_pos = $helend[6];
    $miss_number = $helend[7];
    
    # Obtains the coordinates of CTRR
    if($strand eq "P") {
	$CTRR_end = $helend[2];
	$CTRR_start = $CTRR_end - 4;
    }
    else {
	$CTRR_end = $helend[2];
	$CTRR_start = $CTRR_end + 4;		
    }

    # Obtains the coordinates of the hairpin depending on the number of mismatches
    if($miss_number == 0) {
	$fasta{$id} =~ /$hairloop/;
	$start = $-[0];
	$end = $+[0];
	print OUT "$row\t$start\t$end\t$CTRR_start\t$CTRR_end\n";
    }

    if($miss_number == 1) {
	$gap_pos = $helend[6];
	
	# Obtains the positions of the hair-loop within the upstream region of CTRR
	$wholesequence2 =~ s/(\w{$gap_pos})\w(.*)/$1$2/;
	$wholesequence2 =~ /$hairloop/;
	$start = $-[0];
	$end = $+[0] + 1;
	
	# Obtains the positions of the hair-loop within the entire fasta sequence
	$fasta{$id} =~ /$wholesequence/;
	$start += $-[0];
	$end += $-[0];
	print OUT "$row\t$start\t$end\t$CTRR_start\t$CTRR_end\n";
    }
	
    if($miss_number == 2) {
	my @gap_pos2 = split(/\|/, $helend[6]);
	my $pos2 = $gap_pos2[1] - $gap_pos2[0] - 1;
	
	# Obtains the positions of the hair-loop within the upstream region of CTRR
	$wholesequence2 =~ s/(\w{$gap_pos2[0]})\w(\w{$pos2})\w(.*)/$1$2$3/;
	$wholesequence2 =~ /$hairloop/;
	$start = $-[0];
	$end = $+[0] + 2;

	# Obtains the positions of the hair-loop within the entire fasta sequence
	$fasta{$id} =~ /$wholesequence/;
	$start += $-[0];
	$end += $-[0];
	print OUT "$row\t$start\t$end\t$CTRR_start\t$CTRR_end\n";
    }
}

close(HELEND);
close(OUT);

system ("rm temporal");
