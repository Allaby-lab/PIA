################ FileManipulations.pm ################
# Subroutines for manipulating file contents
#
# Contents:
# Num	Name						Description	
# 1. 	Sort_columns				Sort file by contents of first 5 columns	
# Written by Roselyn Ware
#
# Last updated 23.09.15
#
# To use this module add the following two lines of code to the program (removing the "#")
#
# use lib './modules/';
# use FileManipulations;
#
###############################################

#Declare package in order to preserve namespace
package FileManipulations;

use strict;
use warnings;
use lib './modules/';

### 1. Sort_columns ###
# Sort file by contents of first 8 columns

sub Sort_columns {
	my ($filename)=@_;
	open (UNSORTED, $filename) or die "Cannot open file for sorting\n";
	my @unsorted=<UNSORTED>;
	close UNSORTED;
	my $header= shift(@unsorted);
	$header=~ s/[\t]+/\t/g;#remove multiple tabs
	$header=~ s/\s+$/\n/g;# remove space at end of header (whilst retaining newline character)
	my @unsorted2 = map {chomp; [split /[\t]+/, $_] }@unsorted ; #read each row into an array
	my @sorted=();
	@sorted = sort {$a->[0].$a->[1].$a->[2].$a->[3].$a->[4].$a->[5].$a->[6].$a->[7] cmp $b->[0].$b->[1].$b->[2].$b->[3].$b->[4].$b->[5].$a->[6].$a->[7]} @unsorted2; # sort the rows (numerically) by 3rd column

	open (SORTED, ">$filename") or die "Cannot open file for sorting\n";
	print SORTED $header;
	for (@sorted) {
  		print SORTED join("\t", @$_)."\n"; 
	}
	close SORTED;
	tie my @array, 'Tie::File', $filename or die $!;
	chomp $array[-1];
	untie @array;
}

# All perl modules must finish with "1;" in order to evaluate to true
1;








