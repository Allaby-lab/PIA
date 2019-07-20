################ FileMerge.pm ################
#
# Subroutines to use when merging "Original" and "Adjusted" files
#
# Contents:
# Num	Name						Description	
# 1. 	Extract_data				Extract data from "Original" and "Adjusted" files
# 2. 	Name_entries				Name Entries in "Original" and "Adjusted" files
# Written by Roselyn Ware
#
# Last updated 23.09.15
#
# To use this module add the following two lines of code to the program (removing the "#")
#
# use lib './modules/';
# use FileMerge;
#
###############################################

#Declare package in order to preserve namespace
package FileMerge;

use strict;
#use warnings;
use lib './modules/';

### 1. Extract_data ###
# Extract data from "Original" and "Adjusted" files

sub Extract_data {

	my @taxa2 = ();
	unlink "Original1.txt"; #would loop if already present
	unlink "Adjusted1.txt"; #would loop if already present


	open(ORIGINAL,">Original1.txt") or die "Cannot open Original file\n";
	my(@handles);
	for(<Output/CollapsedPercentAdjustedGenomeSize*\.txt>){
		open($handles[@handles],$_) or die "Cannot open input file\n";
	}
	my $atleastone=1;
	while ($atleastone){
		$atleastone=0;
  		for my $op(@handles){
   			if ($_=readline($op)){
      		my @col=split(/\t/);
      		#print ORIGINAL "$col[1]\t";
      		  my $output= $col[1];
      			chomp $output;
      			print ORIGINAL "$output\t";
      		push @taxa2, "$col[0]\t";
      		$atleastone=1;
    		}
 		 }
  		print ORIGINAL"\n";
	} 
	close(ORIGINAL); 
	undef @handles; #closes all files
 
	open(ADJUSTED,">Adjusted1.txt")or die "Cannot open Adjusted file\n";
	my(@handles2); 
	for(<Output/CollapsedPercentAdjustedGenomeSize*\.txt>){
  		open($handles2[@handles2],$_) or die "Cannot open input file\n";
	}
	my $atleastone2=1;
	while ($atleastone2){
 		$atleastone2=0;
  		for my $op2(@handles2){
   			if ($_=readline($op2)){
     			my @col=split(/\t/);
      			my $output= $col[2];
      			chomp $output;
      			print ADJUSTED "$output\t";
      			$atleastone2=1;
   			 }
  		}
		print ADJUSTED"\n";
	} 
	close(ADJUSTED); 
 	undef @handles2; #closes all files 
 	return @taxa2;
 }



### 2. Name_entries ###
# Name Entries in "Original" and "Adjusted" files

sub Name_entries{
	my ($taxaref)=@_;
	my @taxa2= @$taxaref;
	my @unique = do { my %seen; grep { !$seen{$_}++ } @taxa2 };
	unlink "OriginalBasic.txt"; #would loop if already present
	unlink "AdjustedBasic.txt"; #would loop if already present

	open(ORIGINALONE,"Original1.txt") or die "Cannot open Original1 file\n";
	open(ORIGINAL,">Output/OriginalBasic.txt") or die "Cannot open Original file\n";

	my $numberTaxa=(scalar @unique);
	my $i=0;
		while (my $line=<ORIGINALONE>){
			#my $t = ($i-1);
			chomp $line;
			my $row= $unique[$i]."\t".$line."\n";
			if ($row =~ /^[a-zA-Z]/){
				print ORIGINAL $row;
			}
			$i++;
		}

	close(ORIGINALONE);
	close(ORIGINAL); 

	open(ADJUSTEDONE,"Adjusted1.txt") or die "Cannot open adjusted file\n";
	open(ADJUSTED,">Output/AdjustedBasic.txt") or die "Cannot open adjusted file\n";

	$i=0;
		while (my $line=<ADJUSTEDONE>){	 
			#my @array=();
			chomp $line;
			my $row= $unique[$i]."\t".$line."\n";
			if ($row =~ /^[a-zA-Z]/){
				print ADJUSTED $row;
			}
			$i++;
		}

	close(ADJUSTEDONE);
	close(ADJUSTED); 
	unlink "Original1.txt"; #would loop if already present
	unlink "Adjusted1.txt"; #would loop if already present
}


 
# All perl modules must finish with "1;" in order to evaluate to true
1;








