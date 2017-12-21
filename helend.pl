# MODULE 1 FOR TALLO-ASA ALGORITHM 

# Modified by Pablo Eduardo Garcia, Isaac Rodríguez
# AT oct-18-2012
# AT sep-08-2017

# Description: In a multi-fasta file the prgram searches for a CTRR (R stands for A or G) 
# sequence and finds any posible harpin-loop structure in the 30 bp upstream region.

# Input: a fasta file with one or more sequences.

# helend.pl 
# step 1
# Helend: program searching for helitron ending structure
# Author: Lixing Yang, Department of Genetics, University of Georgia, Athens GA, 30602, USA

#output file name: helendout.txt

# Out put file format:
# A: Sequence ID
# B: Positive or Negative sequence
# C: CTRR position
# D: Whole sequence
# E: Palindrome structure
# F: Palindrome length (pair)
# G: Gap(s) position within the Whole sequence (D)
# H: Mistake number

# Example: perl helend.pl fasta_file.fas n/y output.txt
use strict;
use Getopt::Long;

# initialize file
my $input = shift;
my $is_strand = shift;
my $outputfile = shift;

open (INPUT, "$input") or die "Cannot open the input file \n";
open (OUTPUT,">$outputfile");
print OUTPUT "";
close OUTPUT;

#define argument and initilize
my $hairpinnummin = 5;
my $hairpinnummax = 11;
my $misnum = 2;
my $fraglen = 4;
my $distance = 5;
my $buffer = '';
my $CTRRpos = 0;
my $CTRRposr = 0;
my $position = 0;
my $id;
# Added by Pablo
my $gap_pos = 0;
my %trunked1;
my %trunked2;
my $counter_seq=1;
my %seqArray;
my $missLess = 0;

# #############################################################################
# Table of equivalences or dictionary.
# It will tell how many positions can be assessed in the window size selected
my %equiDIC = ("28" => "4", "29" => "5", "30" => "6", "31" => "7", "32" => "8", "33" => "9", "34" => "10");

# #############################################################################
# Read fasta file and load all DNA sequences into a hash
while (<INPUT>) {
    chomp $_;
    next if($_=~/^#/);
    
    if ($_ =~ /^>/) {
	$id = $_;
	$id =~ s/^>//;
	$|=1;
    } else{
	$seqArray{$id} .= $_;
    }
}
close(INPUT);
print "\n";

# variable that will help define the longest hairpin when the search includes
# mismatches
my $hairPSel;

# #############################################################################
# Main part:
# Each read saved in the hash will be processed looking for a CTRRs and
# hairpins structures
for my $Masterkey(sort {$a<=>$b} keys %seqArray){
    print "\rAnalyzing sequence number: $counter_seq";
    $counter_seq++;
    $buffer = $seqArray{$Masterkey};
    $position = 0;
    $CTRRpos = 0;
    $CTRRposr = 0;	
    $id = $Masterkey;
    $hairPSel = 0;
    # Maximize the window where to look for a hairpin. 
    # 34 = 11 per hairpin arm + 6 bubble + 11 linker - 5 fixed of the linker
    my $upstreamnum = 34;
	
    # Search in positive sequence for perfect hairpin
    my $upstream=0;
    my $wholesequence = 0;
    my ($palindrome, $palinlength, $match, $revmatch, $loop) = (0,0,0,0,0);
    
    my $GCcont = 0;
    my $GCcontPast = 0;
	
    while ($buffer =~ /CT[AG][AG]/gi) {
	$CTRRpos = $position + $-[0] + $fraglen;
	
	# Check if the helitron is larger than 28nt but smaller than 34, 
	# if so, it will determine the max size of the window to look for 
	# the hairpin
	# Re-define the value of $upstreamnum
	my $redf = $-[0] - $distance;
	if($redf < $upstreamnum && $redf >= 28){
	    $upstreamnum = $redf;
	}
		
	$upstream = substr($buffer, ($-[0] - $upstreamnum - $distance), $upstreamnum) if (($-[0] - $upstreamnum - $distance +1) > 0);
	$wholesequence = substr($buffer, ($-[0] - $upstreamnum - $distance), ($upstreamnum + $fraglen + $distance)) if (($-[0] - $upstreamnum - $distance +1) > 0);
	
	# Skip to the next sequence if the helitron to be analyzed is 
	# really small
	next if(!$upstream);
	($palindrome, $palinlength, $match, $revmatch, $loop) = palindrome($hairpinnummin, $hairpinnummax, $upstreamnum, $upstream, $missLess, $Masterkey,$wholesequence);
	
	$gap_pos = 0;
	my $upstreamtrunk = 0;
	my $check = 0;
	my $checkseq = 0;
	$upstreamtrunk = $upstream;
				
	my $PalEval = $palindrome;
	$PalEval =~ s/\*//g;
	if($upstream =~ /$PalEval([CcAaGgTt]{0,15})/){
	    $check = length($1);
	    $checkseq = $1;
	} else{
	    $check = 100;       # high random number to allow continue the search
	    $checkseq = "NAN";
	}
	# Check if the distance of the hairpin from the linker is longer than 
	# 11bp, if so, the hairping found is discarded
	# REMEMBER: 
	# minimum distance of the linker = 5bp, maximum = 11, that's why check 
	# has to be smaller or eaqual to 6
	$check = 0;
	
	###########################################################
	# Search for hairpin allowing 1 mismatch 
	
	if (!$palindrome and $misnum>0){
	    $GCcont = 0;
	    $GCcontPast = 0;
	    $missLess = 1;
	    $hairPSel = 0;
	    undef(%trunked1);
	    trunk1($upstream,$upstreamnum);
	    
	    for my $key(sort {$a<=>$b} keys %trunked1) {
		my $value = $key;
		my $posKey = $trunked1{$key};
		my ($palindromet, $palinlengtht, $match, $revmatch, $loop) = (0,0,0,0,0);
		($palindromet, $palinlengtht, $match, $revmatch, $loop) = palindrome($hairpinnummin, $hairpinnummax, $upstreamnum, $value, $missLess);
		
		# More info, check in the section of the 2 mismatches
		my $PalEval = $palindromet;
		$PalEval =~ s/\*//g;
		if($value =~ /$PalEval([CcAaGgTt]{0,15})/){
		    $check = length($1);
		    $checkseq = $1;
		} else{
		    $check = 100;
		    $checkseq = "NAN";
		}		
		
		$palindromet = 0 if $check > 6;		
		
		if ($palindromet and $check <= 6) {
		    my $gc_val = $match =~ tr/GCgc/GCgc/;
		    $GCcont= 100 * $gc_val/length($match);
		    
		    if($GCcont > $GCcontPast){
			$hairPSel = length($palindromet);
			($palindrome, $palinlength) = ($palindromet, $palinlengtht);
			$upstreamtrunk = $value;
			$gap_pos = $posKey;
			$gap_pos =~ s/(\d+)\|(\d+)/($1-1)."|".($2)/e;
			$GCcontPast = $GCcont;
		    } elsif($GCcont == $GCcontPast && length($palindromet) > $hairPSel){
			$hairPSel = length($palindromet);
			($palindrome, $palinlength) = ($palindromet, $palinlengtht);
			$upstreamtrunk = $value;
			$gap_pos = $posKey;
			$gap_pos =~ s/(\d+)\|(\d+)/($1-1)."|".($2)/e;
			$GCcontPast = $GCcont;
		    }
		}
		$check = 0;
	    }
	}
	
	# #############################################################################
	# Serch for hairpin allowing 2 mismatches 
	if (!$palindrome and $misnum>1) {
	    $GCcont = 0;
	    $GCcontPast = 0;
	    $missLess = 2;
	    $hairPSel = 0;
	    undef(%trunked2);
	    trunk2($upstream, $upstreamnum);
	    
	    for my $key ( sort {$a<=>$b} keys %trunked2) {
		my $value = $key;
		my $posKey = $trunked2{$key};
		my ($palindromet, $palinlengtht, $match, $revmatch, $loop) = (0,0,0,0,0);
		($palindromet, $palinlengtht, $match, $revmatch, $loop) = palindrome($hairpinnummin, $hairpinnummax, $upstreamnum, $value, $missLess);
		
		# Fixing a problem... Instead of looking for the second half of
		# the hairpin, which can have more than one match, the regular 
		# expression will search por the position of the complete 
		# hairpin, in order to determine the true length of the linker 
		# segment
		my $PalEval = $palindromet;
		$PalEval =~ s/\*//g;
		if($value =~ /$PalEval([CcAaGgTt]{0,15})/){
		    $check = length($1);
		    $checkseq = $1;
		} else{
		    $check = 100;
		    $checkseq = "NAN";
		}
		
		$palindromet = 0 if $check > 6;
		if ($palindromet and $check <= 6) {
		    my $gc_val = $match =~ tr/GCgc/GCgc/;
		    $GCcont= 100 * $gc_val/length($match);
		    if($GCcont > $GCcontPast){
			$hairPSel = length($palindromet);
			($palindrome, $palinlength) = ($palindromet, $palinlengtht);
			$upstreamtrunk = $value;
			$gap_pos = $posKey;
			$gap_pos =~ s/(\d+)\|(\d+)/($1-1)."|".($2)/e;
			$GCcontPast = $GCcont;
		    } elsif($GCcont == $GCcontPast && length($palindromet) > $hairPSel){
			$hairPSel = length($palindromet);
			($palindrome, $palinlength) = ($palindromet, $palinlengtht);
			$upstreamtrunk = $value;
			$gap_pos = $posKey;
			$gap_pos =~ s/(\d+)\|(\d+)/($1-1)."|".($2)/e;
			$GCcontPast = $GCcont;
		    }
		}
		$check = 0;
	    }
	}
	if ($palindrome) {
	    my $mis = length($upstream) - length($upstreamtrunk);
	    
	    open OUTPUT,">>$outputfile";
	    print OUTPUT "$id\tP\t$CTRRpos\t$wholesequence\t$palindrome\t$palinlength\t$gap_pos\t$mis\n";
	    close OUTPUT;
	} 
    }
}
print "\n\n";
END;

# #############################################################################
# Subroutines

# Truncate upstream sequences by one base to find perfect palindrome
sub trunk1{
    my ($Tupstream,$Tupstreamnum) = @_;
    my $i;
    
    for ($i=2; $i<$Tupstreamnum; $i++) {
	my $var1 = substr($Tupstream,0,$i-1);
	my $var2 = substr($Tupstream,$i);
	my $Tupstreamgap = $var1.$var2;
	$trunked1{$Tupstreamgap} = $i if(!exists($trunked1{$Tupstreamgap}));
    }
}

# #############################################################################
# Truncate upstream sequences by two bases to find perfect palindrome
sub trunk2 {
    my ($Tupstream,$Tupstreamnum) = @_;
    my ($i,$j) = (0,0);
    
    for ($i=2; $i<$Tupstreamnum; $i++) {
	my $var1 = substr($Tupstream,0,$i-1);
	my $var2 = substr($Tupstream,$i);
	my $Tupstreamgap = $var1.$var2;
	
	for ($j=$i+1; $j<($Tupstreamnum-1); $j++) {
	    my $var3 = substr($Tupstreamgap,0,$j-1);
	    my $var4 = substr($Tupstreamgap,$j);
	    my $Tupstreamgapj = $var3.$var4;
	    $trunked2{$Tupstreamgapj} = $i."|".$j if(!exists($trunked2{$Tupstreamgapj}));
	}
    }
}

sub palindrome {
    my ($hairpinnummin, $hairpinnummax, $upstreamnum, $test, $missRess) = @_;
    
    my $prevGC = 0;
    my $GCnow = 0;
    my $match = 0;
    my $revmatch = 0;
    my $palindromesub = 0;
    my $palinlength = 0;
    my $matcht = 0;
    my $revmatcht = 0;
    my $palindromesubt = 0;
    my $palinlengtht = 0;
    my $loopt = 0;
    my $loop = 0;
    
    my ($i, $j);
    # max number of usable positions to look for the hairpin 
    # (when the stem is 11bp long in a window of 34 nt)
    my $comb = $equiDIC{$upstreamnum};
    
    for ($i=$hairpinnummax; $i >= $hairpinnummin; $i--) {
	for ($j = 0; $j<=($comb+2*($hairpinnummax - $i)); $j++) {
	    $matcht = substr($test, $j, $i);
	    $revmatcht = reverse($matcht);
	    $revmatcht =~ tr/CAGTcagt/GTCAgtca/;
	    
	    if ($test =~ /($matcht)([CcAaGgTt]{2,6})($revmatcht)/) {
		$loopt = $2;
		$palindromesubt = $matcht."*".$loopt."*".$revmatcht;
		$palinlengtht = length($palindromesubt)-2;
		$test =~ /($matcht)([CcAaGgTt]{2,6})($revmatcht)([CcAaGgTt]{0,15})/;
		my $lengthFilter = length($4);
		
		my $gc_count = $matcht =~ tr/GCgc/GCgc/;
		$GCnow = 100 * $gc_count/length($matcht);
		
		if($GCnow > $prevGC && $lengthFilter <= 6){
		    $palindromesub = $palindromesubt;
		    $palinlength = $palinlengtht;
		    $match = $matcht;
		    $revmatch = $revmatcht;
		    $loop = $loopt;
		    $prevGC = $GCnow;
		} elsif($GCnow == $prevGC && $lengthFilter <= 6){
		    # If two or more hairpins have the same GC content, 
		    # the longest hairpin will be selected
		    if($palinlengtht > $palinlength){
			$palindromesub = $palindromesubt;
			$palinlength = $palinlengtht;
			$match = $matcht;
			$revmatch = $revmatcht;
			$loop = $loopt;
			$prevGC = $GCnow;
		    }	
		}
	    }
	}
    }
    return ($palindromesub, $palinlength, $match, $revmatch, $loop);   
}
