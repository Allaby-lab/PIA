################ FileChecks.pm ################
#
# Subroutines to use when working with any file types
#
# Contents:
# Num	Name						Return?		Description	
# 1. 	process_help				N			Decide if help file should be printed. -h flag
# 2.	process_subset_size			Y			Ensure that a subset_size has been given as an argument. -k flag
# 3.	process_disposal_cuttoff	Y			Ensure that a disposal cuttoff has been given as an argument. -t flag
# 4. 	process_filename			Y			Ensure that a filename has been given as an argument. -f flag
# 5. 	check_numeric				N			Check that input is numeric
# 6.	check_exists				N			Ensure that file exists
# 7.	check_name					N			Ensure that file ends in user entered file type- INCOMPLETE
# 8.	open_file					N			Open File and return filehandle
# 9. 	create_output_file			Y			Create an output file name and check if this is accptable
# 10. 	print_output 				N			Print file
# 12. 	output_directory			N			Creates directory named "Output" in current working directory
# Written by Roselyn Ware
#
# Last updated 23.9.2015
#
# To use this module add the following two lines of code to the program (removing the "#")
#
# use lib './modules/';
# use FileChecks;
#
###############################################

#Declare package in order to preserve namespace
package FileChecks;

use strict;
use warnings;

### 1. process_help ###
# Decide if help file should be printed. -h flag

sub process_help {
#include the following lines in main body
#my %options=();
#getopts("h", \%options);
#	if ($options{h}){
#	FileChecks::process_help();															#FileChecks.pm	
#	}	
		my ($helpfile)=@_;
		open(DATA, "$helpfile") or die "Couldn't open $helpfile, $!";
			while(<DATA>){
			print "$_";
		}
		exit;
}



### 2. process_subset_size ###
# Ensure that a subset_size has been given as an argument. -k flag

sub process_subset_size {
#include the following lines in main body
#my %options=();
#getopts("k:", \%options);
#process_help($options{k});

	#subset_size
	my ($subset_size)=@_;
		# Check that the subset size has been included. if not, default = 1000
		# Get and clean subsetsize

	#if no input specified, default = 10
	if (!$subset_size){
		$subset_size=10;
	}
	check_numeric ($subset_size, "subset size");
	print "Subset size:\t $subset_size\n";
	return $subset_size;
}	



### 3. process_disposal_cuttoff ###
# Ensure that a disposal cuttoff has been given as an argument. -t flag

sub process_disposal_cuttoff {
#include the following lines in main body
#my %options=();
#getopts("t:", \%options);
#process_help($options{t});

	#cutoff
	# Get and clean Threshold
	my ($threshold)=@_;

	#if no input specified, default = 10
	if (!$threshold){
		$threshold=10;
	}
	check_numeric ($threshold, "threshold");
	print "Threshold:\t $threshold\n";
	return $threshold;
}



### 4. process_filename ###
# Ensure that a filename has been given as an argument. -f flag

sub process_filename {
#include the following lines in main body
#my %options=();
#getopts("f:", \%options);
#process_help($options{f});

	#File name
	my ($inputFileName)=@_;
	if ($inputFileName) {
		# Check that the filename has been included
		my($USAGE) = "$0 Requres a filename to be included as an argument\n\n";
		unless($inputFileName){
			print $USAGE;
			exit;
		}

		# Get and clean input filename
		chomp $inputFileName;
		#print "$inputFileName\n";
		return $inputFileName;
	}
	unless($inputFileName){
		print "Enter Input Filename. \nFor help with arguments; use the -h flag.\n";
		exit;
	}
}	



### 5. check_numeric ###
# check if input is numeric

sub check_numeric {
	#check if present
	my ($in,$name)=@_;

	#check if numeric
	if ($in =~/^\d+$/){
		#do nothing
	}else{
		print "The $name must be entered as an integer\n";
		exit;
	}
}



### 6. check_exists ###
# Ensure that file exists

sub check_exists {
	my ($inputFile) =@_;
	
	unless (-e $inputFile) {
	    print "Cannot open input file\n";
	    exit;
	}
}



### 7. check_name ###
# Ensure that file ends in user entered file type
###INCOMPLETE###

sub check_name {
	my ($fileName) = @_;
	
	if ($fileName =~ /.fasta/) {
		print "Input file is named correctly\n";
	} 
	elsif ($fileName =~ /.fas/) {
		print "Input file is named correctly\n";
	} else {
		print "Inproper file name; use .fasta or .fas\n";
		print $fileName;
		exit;
	}
}



### 8. open_file ###
# Open File and return filehandle

sub open_file {
my ($closedFile)= @_;

open (INFILE, $closedFile) || die "Could not open input file: $!\n";
my @rawfasta = <INFILE>;
close INFILE;
return @rawfasta;
}



### 9. create_output_file ###
# Create an output file name

sub create_output_file {
#include the following lines in main body
#my %options=();
#getopts("o:c", \%options);
#create_output_file(in.file, subset.size,($options{o}),($options{c}));
# Print output name in format inputfile_subset_k.fasta
# Input is input file name and subset size. addititional options are the output file name (flagged) and the continuous flag
	my ($inFileNameFormat, $k, $outputFileName,$continuous)=@_;

	if ($outputFileName) {
		print "Output file name is: $outputFileName\n";
		return $outputFileName;
	} else {
		$inFileNameFormat=~ s{\.[^.]+$}{};
		my $newOutputFileName= $inFileNameFormat.".".$k.".fasta";
		
		if ($continuous){
			return $newOutputFileName;
		} else {
			#Check that file name is acceptable, else allow them to submit their own filename
			print "Output File Name is: $newOutputFileName\nIs this acceptable? [Y/N]";
		
			my $test=<STDIN>;
			chomp $test;
			if ($test eq "Y"){
				return $newOutputFileName;	
			}if ($test eq "N"){
				print "Please enter required output file:\n";
				$newOutputFileName= <STDIN>;
				print "New output file name is: $newOutputFileName";
				return $newOutputFileName;
			}else{
 	  			print "\n\n=====================================\n";
 		   		print "You Have Entered an Incorrect Option\n";
 		   		print "Valid Options are [Y/N]\n";
  		  		print "=====================================\n\n";
  		  		#Call subroutine again. 
    			&create_output_file;
			}
		}
	}
}

### 10. print_output ###
# Print to output file

sub print_output {
	my (@fastaSubset)= @{$_[0]};
	my $outputFileName= $_[1];
	open my $fh, ">", $outputFileName or die("Could not open file. $!");
	print $fh @fastaSubset;
	close $fh;
}


### 11. output_directory ###
# make output directory named "Output"
sub output_directory {
    my $directory = "Output";
   unless(mkdir $directory) {
   
     print "Overwriting directory:$directory/\n";
    }
}


# All perl modules must finish with "1;" in order to evaluate to true
1;