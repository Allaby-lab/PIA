################ TreeOfLife.pm ################
#
# Subroutines to use when working with names.dmp and nodes.dmp
#
# Contents:
# Num	Name						Description
# 1.	check_taxa_names				
# 2. 	find_code_from_name			Look up taxa name in nodes file, and find code
# 3.	find_classification			Look up taxa code in names file, return classification
# 4. 	find_tax_level				Use old code to find taxonomic level and new code
# 5.	get_lineage_array			Backfill megan file and return lineage array
# 6. 	format_megan_ss				Backfill megan data and return fully formatted tab separated spreadsheet ready for genome size adjustment	
# 7. 	format_output				Format final output files ready for analysis in R
# Written by Roselyn Ware
#
# Last updated 23.9.2015
#
# To use this module add the following two lines of code to the program (removing the "#")
#
# use lib './modules/';
# use TreeOfLife;
#
###############################################

#Declare package in order to preserve namespace
package TreeOfLife;

use strict;
use warnings;
use lib './modules/';
use FileChecks;
use TreeOfLife;
use Data::Dumper qw(Dumper);


### 1. check_taxa_names ###
# If Taxa name does not exist in names file, remove line from analysis and write line to new file

sub check_taxa_names{
	my ($megandata,$namesfile)=@_;
	my @megandata=@$megandata;
	my @megandata2=();
	my @names=();

	print "****Checking taxa names exist in names.dmp file****\n";
	foreach my $line (@megandata){
		open (NAMES, $namesfile)or die "Cannot open names file\n";
		my @data = ();
		@data = split (/\t/, $line);
		my $taxaname= $data[0];
		if(grep{/^[0-9]+\t\|\t$taxaname/} <NAMES> ){
			my $row= join("\t", @data);
			push @megandata2,"$row";
			print $taxaname." exists\n";
		
		}else {
			print "$taxaname not found in names.dmp file.\nIt will be excluded from the analysis\nIf these data are important, check spelling, etc.\n";
			next;
		}
		close NAMES;
	}
	
	return @megandata2;
}	


	
### 2. find_code_from_name ###
# Look up taxa name in nodes file, and find code

sub find_code_from_name {
	my ($name_ref, $namesfile)= @_;
	my @name = $name_ref;
	my $name= join('',@name);	
	open (NAMES,  $namesfile) or die "Cannot open names file\n";

	while (my $line = <NAMES>){
		my @nameline = ();
		@nameline = split (/\t/, $line);
		my $namecode = $nameline[0];
		my $nameclass = $nameline[2];
		if ($nameclass eq $name){
			return ($namecode, $nameclass);
		}
	}
	
	close (NAMES);
}



### 3. find_classification ###
# Look up taxa code in names file, return classification

sub find_classification{
	my ($code_ref, $namesfile)=@_;
	my @code =$code_ref;
	my $code= join('',@code);	
	open (NAMES, $namesfile) or die "Cannot open names file\n";
	while (my $line = <NAMES>){
		my @nameline = ();
		@nameline = split (/\t/, $line);
		my $namecode = $nameline[0];
		my $nameclass = $nameline[2];
		if ($namecode eq $code){
			return ($nameclass);
		}
	}
	close (NAMES);
}



### 4. find_tax_level ###
# Use old code to find taxonomic level and new code

sub find_tax_level {
	my ($code_ref, $nodesfile)=@_;
	my @code =$code_ref;
	my $code= join('',@code);
	open (NODES,  $nodesfile) or die "Cannot open nodes file\n";
	while (my $line = <NODES>){
		my @nodeline = ();
		@nodeline = split (/\t/, $line);
		my $oldcode = $nodeline[0];
		my $newcode = $nodeline[2];
		my $taxlevel= $nodeline[4];
		if ($oldcode eq $code){
			return ($newcode, $taxlevel, $oldcode);
		}	
	}	
	close(NODES);
}


### 5. get_lineage_array ###
# Backfill megan file and return lineage array

sub get_lineage_array {
	my ($megandata_ref,$nodesfile,$namesfile)= @_;
	my @megandata=@$megandata_ref;
	
#### separate header line ####
	my $header= shift @megandata;
	$header =~ s/Taxa\t//g;
	
#### Backfill megan file and return lineage array
	
	my @taxalevels= ( "superkingdom","kingdom","phylum", "class" , "order","family", "genus","species");
	my $taxalevels= "superkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\t";
	my $toptaxa= "superkingdom";
	my $species="species";
	my $genus="genus";
	my $family="family";
	my $order="order";
	my $class="class";
	my $phylum="phylum";
	my $kingdom="kingdom";
	my $superkingdom="superkingdom";
	

	
	my @returnarray=();
	push @returnarray, $taxalevels;
	
	
	foreach my $line  (@megandata){
		my @outarray=("NISS\t","NISS\t","NISS\t","NISS\t","NISS\t","NISS\t","NISS\t","NISS\t");
		my @row= ();
		@row = split (/\t/, $line);
		my $namecode = $row[0];
		#open names, find matching entry, return code and classification
		my ($startcode, $classification) =TreeOfLife::find_code_from_name($namecode, $namesfile);		#find_code_from_name
		#open nodes, find matching code, return new code and tax level
		my ($nextcode, $taxlevel)= TreeOfLife::find_tax_level($startcode, $nodesfile);				#find_tax_level
		#if $taxlevel is in @taxalevels, write classification to array
		print "$namecode\t$taxlevel\n";
	
		if ($taxlevel eq $species){
			chomp $classification;
			splice @outarray,7,1, "$classification\t";
		} elsif ($taxlevel eq $genus){
			chomp $classification;
			splice @outarray,6,1, "$classification\t";
		} elsif ($taxlevel eq $family){
			chomp $classification;
			splice @outarray,5,1, "$classification\t";
		} elsif ($taxlevel eq $order){
			chomp $classification;
			splice @outarray,4,1, "$classification\t";
		} elsif ($taxlevel eq $class){
			chomp $classification;
			splice @outarray,3,1, "$classification\t";
		} elsif ($taxlevel eq $phylum){
			chomp $classification;
			splice @outarray,2,1, "$classification\t";
		} elsif ($taxlevel eq $kingdom){
			chomp $classification;
			splice @outarray,1,1, "$classification\t";
		} elsif ($taxlevel eq $superkingdom){
			chomp $classification;
			splice @outarray,0,1, "$classification\t";
		} else {
		}
		
		if ($taxlevel eq $toptaxa){
		} else{
			do {
				#here is where the repeat code needs to go
				($nextcode, $taxlevel, $startcode)= TreeOfLife::find_tax_level($nextcode, $nodesfile);	#find_tax_level
				#my $newmatch = grep {/$taxlevel/} @taxalevels;
				#if ($newmatch){
						if ($taxlevel eq $species){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,7,1, "$classification\t";
						} elsif ($taxlevel eq $genus){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,6,1, "$classification\t";
						} elsif ($taxlevel eq $family){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,5,1, "$classification\t";
						} elsif ($taxlevel eq $order){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,4,1, "$classification\t";
						} elsif ($taxlevel eq $class){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,3,1, "$classification\t";
						} elsif ($taxlevel eq $phylum){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,2,1, "$classification\t";
						} elsif ($taxlevel eq $kingdom){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,1,1, "$classification\t";
						} elsif ($taxlevel eq $superkingdom){
							$classification= TreeOfLife::find_classification($startcode, $namesfile);			#find_classification
							chomp $classification;
							splice @outarray,0,1, "$classification\t";
						} else {
					}
			
				} until ( $taxlevel eq $toptaxa);
		}
		
		chomp  @outarray;
		

		my $outarray=join('',@outarray);
		chomp $outarray;
		push @returnarray,"$outarray";
		
	}
	return @returnarray;
}



### 6. format_megan_ss ###
# Backfill megan data and return fully formatted tab separated spreadsheet ready for genome size adjustment	

sub format_megan_ss {
	my ($megandata_ref,$count)= @_;
	my @megandata=@$megandata_ref;
	print "$count\n";
	#### Backfill megan file and return fully formatted tab separated spreadsheet ready for genome size adjustment	
	open (OUT, ">Output/FormattedMeganExtraction".$count.".txt") or die "Cannot write outname ".$count." file\n";
	foreach my $line (@megandata){
		print OUT $line;
		print OUT "\n";		
	}
	close OUT;
}



### 7. format_output ###
# Format final output files ready for analysis in R

sub format_output{
	my ($string, $nodesfile,$namesfile,$lineagesarrayfinal, $output)=@_;
	my @lineagesarrayfinal=@$lineagesarrayfinal;
	open(ORIGINAL,$output.$string."Basic.txt") or die "Cannot open ".$string." file\n";
	my @original= <ORIGINAL>;
	my @line=();
	foreach my $line(@lineagesarrayfinal){
		chomp $line;
	}
	foreach my $line(@original){

		chomp $line;
		$line=~ s/^\s+|\s+$//g;
	}
	@original= map {"$lineagesarrayfinal[$_]$original[$_]"} 0 .. $#lineagesarrayfinal;
	close ORIGINAL;
	TreeOfLife::format_megan_ss(\@original,$string);					#TreeOfLife.pm	
}

### 8. format_output_ss ###
# Backfill megan data and return fully formatted tab separated spreadsheet ready for genome size adjustment	

sub format_PIA_ss {
	my ($megandata_ref,$count)= @_;
	my @megandata=@$megandata_ref;
	my @count2=();
	@count2=split (/\_/, $count);
	my $count2= $count2[0];
	#### Backfill megan file and return fully formatted tab separated spreadsheet ready for genome size adjustment	
	open (OUT, ">".$count2."/".$count."_Extended.txt") or die "Cannot write outname ".$count2." file\n";
	foreach my $line (@megandata){
		print OUT $line;
		print OUT "\n";		
	}
	close OUT;
}
### 7. format_output ###
# Format final output files ready for analysis in R

sub format_PIA_output{
	my ($string, $nodesfile,$namesfile,$lineagesarrayfinal, $output)=@_;
	my @lineagesarrayfinal=@$lineagesarrayfinal;
	open(ORIGINAL,$output.$string."_Basic.txt") or die "Cannot open ".$string." file\n";
	my @original= <ORIGINAL>;
	my @line=();
	foreach my $line(@lineagesarrayfinal){
		chomp $line;
	}
	foreach my $line(@original){

		chomp $line;
		$line=~ s/^\s+|\s+$//g;
	}
	@original= map {"$lineagesarrayfinal[$_]$original[$_]"} 0 .. $#lineagesarrayfinal;
	close ORIGINAL;
	TreeOfLife::format_PIA_ss(\@original,$string);					#TreeOfLife.pm	
}

# All perl modules must finish with "1;" in order to evaluate to true
1;








