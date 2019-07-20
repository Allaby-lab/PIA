#!/usr/bin/perl

######################################################
##################### PIA_inner.pl ###################
######################################################
#########  Phylogenetic Intersection Analysis ########
############## Robin Allaby, UoW 2013 ################
######################################################
############## Version 4.0, 2019-07-20 ###############
######################################################

# Edited by Roselyn Ware, UoW 2015
# A method of metagenomic phylogenetic assignation, designed to be robust to partial representation of organisms in the database		
		
# Further edited by Becky Cribdon, UoW 2019.
# - Names searches now take the scientific name from the names file rather than just the first one to match that taxonomic ID.
# - Added a log.
# - Made variable names more consistent and meaningful where appropriate.
# - Removed $blastfile2 and functions that were not fully implemented.
# - Added index files for the names and nodes files.

# Issues:
# - Names search doesn't account for taxonomic trickiness.
# - @scores is used to check whether the current BLAST taxon is new or has already been seen. The taxon name is also used for this. What does @scores add? Why do we not care about BLAST taxa with the same score (BLAST score and E value) as another?
# - Add error messages ($!) to all, uh, error messages.
# - Anything else labelled "to do".
# Please report any problems to r.cribdon@warwick.ac.uk.

# To run, first create a header file from the final FASTA:
# > cat test_sample.fasta | grep ">" > test_sample.header
# Then run PIA_inner.pl giving the header file as -f and a BLAST file as -b.
# > perl PIA_inner.pl -f $filevar.$taxa.collapsed.header -b ../blast/$filevar.$taxa.collapsed.txt
# The header file lists reads by their query number. The BLAST file contains every read (by query number) and a list of reasonable matches (BLAST hits).


	use strict;
	use warnings;
	use lib './Modules';
	use FileChecks; # A bespoke module in /Modules.
	use TreeOfLife; # A bespoke module in /Modules.
	use Getopt::Std;
	use File::Find;
	use Data::Dumper qw(Dumper);
	use FileMerge; # A bespoke module in /Modules.
	use FileManipulations; # A bespoke module in /Modules.
    use DB_File;
    use Fcntl;


######################################################										#####################
########### Check arguments and Input Data ###########										###### Modules ######
######################################################										#####################

##### Get arguments from command line #####
	my %options=();
	getopts('hf:b:c:egn:N:G:', \%options); 														#Getopt::Std
	
	# If other text found on command line, do:
	print "Other things found on the command line:\n" if $ARGV[0];
	foreach (@ARGV)	{
	print "$_\n";
	}

##### Check other commandline arguments and store #####
	#Display help file and exit if -h flag called #		
	my $helpfile="Helpfile_PIA.txt";
	if ($options{h}){
	FileChecks::process_help($helpfile);														#FileChecks.pm	
	}				
							
	#print "\n****Inputs****\n";
	
##### Check header filename and open file #####
	# The header file is simply the read names (headers) extracted from a FASTA file. PIA.pl makes the header file.
	# Ensure that a filename has been given as an argument. -f flag
	# Get and clean input filename
	my $tophitfile =FileChecks::process_filename($options{f});									#FileChecks.pm
	chomp $tophitfile ;
	
##### Check BLAST filename and open file #####
	# Ensure that a filename has been given as an argument. -b flag
	# Get and clean input filename
 	my $blastfile =FileChecks::process_filename($options{b});									#FileChecks.pm
	chomp $blastfile ;
	
##### See if cap input #####
	# The cap impacts the taxon diversity score. Default is 100.	
	# Ensure that a number has been given as an argument. -c flag
	my $capopt =($options{c});																	#FileChecks.pm
	my $cap = 100;
	if ($capopt) { $cap = $capopt;}
	
##### See if extended summary required #####
	#Create extended summary if -e flag called #		
	my $summary;
	if ($options{e}){ $summary = 1;}															#FileChecks.pm		
	
##### Locate nodes.dmp and names.dmp files #####
	# If nodes.dmp and names.dmp aren't in their default location (Reference_files/), use -n and -N to specify where they are.	
	# Check that the path has been included after -l flag
	# Get and clean input path
	my $nodesfile="Reference_files/nodes.dmp";
	if ($options{n}){																			#FileChecks.pm
		$nodesfile=($options{n});
	}
	my $namesfile="Reference_files/names.dmp";
	if ($options{N}){																			#FileChecks.pm
		$namesfile=($options{N});
	}


######################################
########### Start log file ###########
######################################

    my $log_filename = $tophitfile . '_PIA_inner_log.txt';
    open( my $log_filehandle, '>', $log_filename) or die "Cannot open $log_filename for writing.\n$!\n"; # Note: this will overwrite old logs. I didn't see any point in appending old ones when you need to remove the other PIA outputs before running it again anyway.
    use IO::Handle; # Enable autoflush.
    $log_filehandle->autoflush(1); # Set autoflush to 1 for the log filehandle. This means that Perl won't buffer its output to the log, so the log should be updated in real time.
    print $log_filehandle "#####################################\n";
    
    # Print output from the last section to the log.
    print $log_filehandle "Other things found on the command line:\n" if $ARGV[0];
    foreach (@ARGV)	{
    print $log_filehandle "$_\n";
    }
    print $log_filehandle "\n****Inputs****\nHeader file: $tophitfile\nBLAST file: $blastfile\nCap: $cap (default is 100)\n";


##############################################		
########### Set up DBM index files ###########
##############################################
    print "Setting up DBM index files...\n";
    print $log_filehandle "Setting up DBM index files...\n";
    
    # Names index where keys are names
    #---------------------------------
    # This will be used to look up the taxonomic IDs of BLAST hits.
    # TO DO: it does not account for taxonomic trickiness.
    #print "Top hit file for names names: $tophitfile\n";
    #print $log_filehandle "Top hit file for names names: $tophitfile\n";
    my $namesfileDBMnames = $namesfile . $tophitfile . '_names.dbm'; # Including $tophitfile in the DBM filename ensures that every thread makes its own index file. That's a lot of files, but they don't make their own blastfile2.txt any more, and it should prevent read/write conflict. Otherwise I'd have to get into file locks and oh my goodness not now.
        
    if (-e $namesfileDBMnames) {
        print "$namesfileDBMnames already exists. Overwriting.\n";
        print $log_filehandle "$namesfileDBMnames already exists. Overwriting.\n";
        unlink $namesfileDBMnames;
    }
        
    my %namesfileDBMnames = (); # Keys are names (scientific and not, because BLAST contains all sorts) and values are the corresponding taxonomic ID. It is built as a hash and then written to the DBM file.
    tie (%namesfileDBMnames, "DB_File", $namesfileDBMnames, O_RDWR|O_CREAT) or die "Can't open $namesfileDBMnames.\n$!\n"; # Open the DBM file (in read-write mode; O_RDWR). If it doesn't exist, create one (O_CREAT). Note that these options use the Fcntl module.
        
    # Now populate the hash:
    open (my $names_filehandle, $namesfile) or die "Could not open $namesfile.\n$!\n";
    while (1) { # Run this loop until "last" is called.
            my $line = <$names_filehandle>; # Read the next line from the names file.
            if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
    
            my @line = split(/\|/, $line); # Split the line by | characters.
            my $name = $line[1]; # Pick out the name.
            $name =~ s/^\s+|\s+$//g; # Trim whitespace off the start and end.
            $name = substr ($name, 0, 52); # Only store the first 52 characters of the name (counting from 0) because the vast majority of BLAST hits are only 52 characters long, followed by an ellipsis. Many taxa have very long names that only differ at the very end, but considering the PIA only outputs species, in my opinion the effect of assigning to the wrong one will be negligible.
            my $ID = $line[0]; # Also get the ID.
            $ID =~ s/^\s+|\s+$//g;
            
            $namesfileDBMnames{$name} = $ID; # Just assume names are unique (they're not). In reality, some IDs will overwrite each other.
                
            ## TO DO: may be helpful when I finally make it account for taxonomic trickiness!
            ## Have some sort of table where the <angiosperm> label means it will be accepted for data from Viridiplantae but not Metazoa or something.
            ## Field [2] contains a clarifier if there is taxonomic trickiness, but it also clarifies other things too I think.
            ## Store both $name and $ID in the hash. If $name doesn't exist, create it as a key and give $ID as the value. If it already exists, check whether $position already exists. If not, add $position to its list.
            #unless (defined $namesfileDBMnames{$name}) {
            #    $namesfileDBMnames{$name} = "$ID";
            #} else {#print "Position: $position\tID: $taxonomic_ID\tValues at start: $taxonomic_ID_positions{$taxonomic_ID}\n";
            #        if (index ($namesfileDBMnames{$name}, "$ID") == -1) { # If the value for $taxonomic_ID_positions{$taxonomic_ID} does not already contain this $position,
            #        #print "Value doesn't contain $position\n";            
            #        $namesfileDBMnames{$name} = $namesfileDBMnames{$name} . "\t$ID";
            #        #print "Values at end: $taxonomic_ID_positions{$taxonomic_ID}\n";         
            #        }
            #}
    }
    #print Dumper (\%namesfileDBMnames);
    close $names_filehandle;
    untie %namesfileDBMnames;
    
    
    # Names index where keys are IDs
    #-------------------------------
    # This will be used to look up the scientific names of taxa using their IDs.
    # To do: filter out taxa that we know can't be used: strains and sub-categories of viruses that have "no rank" in the nodes file? A BLAST hit will never match them and, being so specific, they should not ever be an intersect. But might take more work to remove them than to check them each time the hash is consulted.
    #print "Top hit file for names IDs: $tophitfile\n";
    #print $log_filehandle "Top hit file for names IDs: $tophitfile\n";
    my $namesfileDBMids = $namesfile . $tophitfile . '_IDs.dbm';
    
    if (-e $namesfileDBMids) {
        print "$namesfileDBMids already exists. Overwriting.\n";
        print $log_filehandle "$namesfileDBMids already exists. Overwriting.\n";
        unlink $namesfileDBMids;
    }
        
    my %namesfileDBMids = (); # Keys are taxonomic IDs and values are scientific names.
    tie (%namesfileDBMids, "DB_File", $namesfileDBMids, O_RDWR|O_CREAT) or die "Can't open $namesfileDBMids.\n$!\n";
    
    # Now populate the hash:
    open ($names_filehandle, $namesfile) or die "Could not open $namesfile.\n$!\n";
    while (1) { # Run this loop until "last" is called.
            my $line = <$names_filehandle>; # Read the next line from the names file.
            if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
    
            my @line = split(/\|/, $line); # Split the line by | characters.
            
            if ($line[3] eq "\tscientific name\t") {
                my $name = $line[1];
                $name =~ s/^\s+|\s+$//g; # Trim whitespace off the start and end.
                #print "Name: $name\n";
                my $ID = $line[0];
                $ID =~ s/^\s+|\s+$//g;
                #print "ID: $ID\n";
                $namesfileDBMids{$ID} = $name; # Store in the hash.
            } # If that line is not for a scientific name, ignore it and move on.
    }
    #print Dumper (\%namesfileDBMids);
    close $names_filehandle;
    untie %namesfileDBMids;
    
    
    # Nodes index where keys are IDs
    #-------------------------------
    #print "Top hit file for nodes: $tophitfile\n";
    #print $log_filehandle "Top hit file for nodes: $tophitfile\n";
    my $nodesfileDBM = $nodesfile . $tophitfile . '.dbm';
    
    if (-e $nodesfileDBM) {
        print "$nodesfileDBM already exists. Overwriting.\n";
        print $log_filehandle "$nodesfileDBM already exists. Overwriting.\n";
        unlink $nodesfileDBM;
    }
        
    my %nodesfileDBM = (); # Keys are taxonomic IDs and values are parent nodes and taxonomic ranks.
    tie (%nodesfileDBM, "DB_File", $nodesfileDBM, O_RDWR|O_CREAT) or die "Can't open $nodesfileDBM.\n$!\n";
    
    # Now populate the hash:
    open (my $nodes_filehandle, $nodesfile) or die "Could not open $nodesfile.\n$!\n";
    while (1) { # Run this loop until "last" is called.
            my $line = <$nodes_filehandle>; # Read the next line from the nodes file.
            if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
    
            my @line = split(/\|/, $line); # Split the line by | characters.
            my $ID = $line[0];
            $ID =~ s/^\s+|\s+$//g; # Trim whitespace off the start and end.
            my $parent_node = $line[1];
            $parent_node =~ s/^\s+|\s+$//g;
            my $rank = $line[2];
            $rank =~ s/^\s+|\s+$//g;
            $nodesfileDBM{$ID} = $parent_node . "\t" . $rank; # Store in the hash. Tab-separate the parent node and rank.
    }
    #print Dumper (\%nodesfileDBM);
    close $nodes_filehandle;
    untie %nodesfileDBM;
    
    
    print "Finished indexing.\n\n";  
    print $log_filehandle "Finished indexing.\n\n";


######################################################
###################### Run PIA #######################
######################################################
    
	my $min_coverage_perc = 95; # The default coverage cut-off is 95%. Reads where the top BLAST hit doesn't have 95% coverage will be discarded.
	
	print "\n****Starting first round of PIA****\n";
	print $log_filehandle "\n\n****Starting first round of PIA****\n";
	
	my ($corename) = PIA($tophitfile, $blastfile, $cap, $namesfile, $min_coverage_perc); # PIA() returns a base name for this sample file. The base name is [header file]_out.
	print $log_filehandle "\nPIA() subroutine finished.\n\n";
	print "\nPIA() subroutine finished.\n\n";

    # IF NOT ACTUALLY RUNNING THE PIA; FOR TESTING
	#my $corename = '76.header_out';


######################################################	
################## Summarise Data ####################
######################################################

##### Extract simple summary from the intersects file
	my $tophitfile2="$corename"."/"."$corename.intersects.txt";
	my $simplesummary = simple_summary($tophitfile2);

##### Extract extended summary (optional) - simple summary with lineage data
	if ($summary){
		extended_summary($simplesummary,$nodesfile,$namesfile);
	}


######################################################
##################### Tidy Up ########################
######################################################

##### Remove un-necessary files	
	unlink("$corename"."/"."temp_blast_entry.txt");
	unlink("$corename"."/"."TEMP");
	unlink("$corename"."/"."hittempfile");
    #unlink $namesfileDBMnames; # Not necessary to delete the index files if running with PIA.pl because it deletes the copy of Reference_files/ that they're in anyway.
    #unlink $namesfileDBMids;
    #unlink $nodesfileDBM;
	
# Finish the log and move it into the output directory.
    print $log_filehandle "This run of PIA_inner.pl is finished.\n\n\n\n";
    close $log_filehandle;
    my $log_final_destination = $corename . '/' . $corename . '_PIA_inner_log.txt';
    system("mv $log_filename $log_final_destination");
    
    system('echo "This run of PIA_inner.pl is finished.\n"');
	
	
######################################################	
######################################################
#################### SUBROUTINES #####################
######################################################
######################################################	

sub extended_summary{
##### Take simple summary and append lineage information 
	# Optional= only if -e flag supplied
	my ($simplesummary, $nodesfile, $namesfile) = @_;
	open (SUM, $simplesummary) or die "Cannot write simplesummary file\n";
	my @original=<SUM>;
	close SUM;
	my @simplesum=();
	@simplesum= split (/\//, $simplesummary);
	my $name= $simplesum[0];
	my @lineagesarrayfinal= TreeOfLife::get_lineage_array(\@original,$nodesfile,$namesfile);	#TreeOfLife.pm
	my $output =$name."/";
	#$name=$name."_Extended_Summary";
	TreeOfLife::format_PIA_output($name,$nodesfile,$namesfile, \@lineagesarrayfinal,$output);				#TreeOfLife.pm

}


sub simple_summary {
#### Create simple summary of output- This is the default. An extended summary can be created if -e flag is supplied
	my ($PIAinputname) = @_; # $PIAinputname is $tophitfile2: basically the header file.
	
	my @intersects_file= FileChecks::open_file($PIAinputname);										#FileChecks.pm
									
    # Get a list of classification intersects. Exclude reads with a taxonomic diversity score below 0.1.					
	my @intersects=();
	
	foreach my $line (@intersects_file){
		my @row= split(/, classification intersect: |, id confidence class: /,$line); # Split on the classification intersect field first (note that it won't match to "most distant classification intersect"), followed by the ID confidence field. This is not an 'or'. It's one after the other, chopping off text from the left and right sides to leave just the classification intersect value in the middle.
		my @check= split(/, taxa diversity score: |, classification intersect: /,$line); # Similarly, leave just the taxa diversity score. 
		if ($check[1]>=0.1 ){ # $check[1] is the taxa diversity score.
			push @intersects, $row[1]; # $row[1] is the classification intersect. Keep it if the taxa diversity score was at least 0.1.
		} else {
			push @intersects, "none found"; # If the taxa diversity score was below 0.1, the classification intersect from now on is "none found".
		}
	}
									
    # Now count the classification intersects.										
	my @intersects_with_counts = (); # This will contain the name of each intersect followed by the number of times it occurs in the intersects file. Unfortunately, it will contain this information once for every occurrence!
	foreach my $intersect(@intersects){
		my @occurrences = grep (/$intersect/, @intersects); # Gather every occurrence of this $intersect, including the current element.
		my $intersect_with_count = join("\t", $intersect, scalar (@occurrences)); # scalar (@occurrences) is the number of occurrences of this $intersect.
		push @intersects_with_counts, $intersect_with_count."\n";
	}
    
	my @name=split (/_/,$PIAinputname); # Pick out a sample name from $PIAinput to use in the output file.
	my $name=$name[0];
    
	open (OUT, ">".$PIAinputname."_Summary_Basic.txt") or die "Cannot write outname file\n";
	print OUT "#Series:\t$name\n"; # Output $name as a header.
    
    # Then output the unique elements in @intersects_with_counts. We only need one statement that Asteridae has 7 reads.
	my @unique = do { my %seen; grep { !$seen{$_}++ } @intersects_with_counts};
	foreach my $uniq (@unique){
		unless ($uniq =~ /none found|all/){
		print OUT $uniq;
		}
	}
    
	close OUT;	
	my $filename=$PIAinputname."_Summary_Basic.txt";
	return $filename;
}



sub PIA {
##### Run phylogenetic intersection analysis

	# Summary
	#--------
	# List the header of each read
	# Use BLAST entries to check coverage for each read
	# For reads that passed the coverage check:
	# 		Get a list of sequences producing significant alignments (the BLAST hits) and their E values
	#		Find the Identities score (coverage) and E value for the top BLAST hit
	#		Find the taxon for each BLAST hit using the names index file
	#		We have a list of unique BLAST taxa. How different are they?
	#		If so, how different? What's their lowest common rank (intersection)?
	#		How robust is the intersection?
	#		Print all of this information to an intersects.txt file
	# Return $corename
    
	my ($tophitfile, $blastfile, $cap, $namesfile, $min_coverage_perc) = @_;
	
	my $corename = $tophitfile . "_out"; # # Generate a core name: a base name for naming things linked to this sample. For example, "test_sample.header_out".
	`mkdir $corename`; # Create a separate directory for outputs (makes folder organisation tidier).
    
    
	# List the header of each read
	#-----------------------------
	open (TOPS, $tophitfile) or die "Cannot open $tophitfile\n"; # Extract headers from the header file. A header is the ID number given to a read, like a name, so a header represents a read. "Header" and "read" are used somewhat interchangably in this code, but it never deals directly with a read. It always works via the header.

	my @tops = (); # @tops will store the headers.
	
	while (my $header = <TOPS>) {
		substr ($header, 0, 1) = ""; # Remove the first character, which is the ">" symbol. These separate FASTA records. We don't need them.
		chomp $header; # Remove any trailing newlines.
		push (@tops, $header);
		
	}
	close TOPS; # Close the header file.
	
	my $number_of_headers = scalar(@tops); # Count the elements in @tops.
	print $log_filehandle "$number_of_headers reads to process.\n";
	

	# Use BLAST entries to check coverage for each read
	#--------------------------------------------------
	my @tops2 = (); # Will contain reads that pass the coverage check. Remember that each header represents a read.

	foreach my $header (@tops) {
			my $header2="Query= ". $header; # $header2 is the read header prefaced with "Query= ", so it looks like the start of a BLAST entry.
			my $read_length_section = `grep -A 2 "$header2" "$blastfile"`; # Search in $blastfile for $header2. Return the two lines after it: a newline and the query length line.
			
			# Coverage is [match length] / [read length]. First, find the read length.
			my @read_length_section = split(/Length=/, $read_length_section); # Split $read_length_section on Length=. Everything after that is the value for read length.
			my $read_length = $read_length_section[1];
			
			# Then, if there is in fact a sequence, find the match length. This is the length of the top BLAST hit and is also called 'query cover':  how much of the query sequence is covered by the top hit (not necessarily covered perfectly). Identities (number of identical bases) are given as a fraction of coverage. If the identities are 86/92, then coverage is 92 bp, but there are 6 mismatches.
			# So, like BLAST, the PIA defines "coverage" as the total length covered regardless of whether there are mismatches or not. This is fine.
			if ($read_length > 0 ){ # To do: can this ever be not the case? Would a BLAST file contain an entry with no hits?
					my $blast_entry_partial=`sed -n -e '/$header/,/Identities =/ p' "$blastfile"`;
					# sed is used as a selective print. -n means "no print" and p here means "print if matches". It will print all lines beginning with the first pattern ("$header") and ending with the second pattern ("Identities =") regardless of what goes on in between.
					# So, this code searches $blastfile for the BLAST entry of $header, but only prints up to the end of the line containing "Identities =".
					my @blast_entry_partial= split(/Identities = [0-9]*\//, $blast_entry_partial); # Split the partial BLAST entry at "Identities =".
					my $match_length_section=$blast_entry_partial[1]; # Pick out what comes after "Identities =", so the end of the final line. It will be something like "63 (95%), Gaps = 0/63 (0%)".
        
					if ( length( $match_length_section || '' )) { # Uh. || basically means "or". So, if there is anything in $match_length_section? I have no idea. To do: understand.
							my @match_length_section = split(/ \([0-9]*\%\), Gaps/, $match_length_section); # Split $match_length_section into the values for identity and gaps.
							my $match_length=$match_length_section[0]; # Pick out just the numerical value of the match length. It's a single number.
							
                            # Calculate percentage coverage and filter out those less than or equal to the minimum (user defined ARGV[2] or default of 100). We calculate percentage ourselves because the percentage in the BLAST output is rounded.
							my $coverage = $match_length / $read_length;
							my $min_coverage = $min_coverage_perc / 100;
							if ($coverage >= $min_coverage) {
									push @tops2, $header;
							} # If this BLAST entry, this query number, has sufficient coverage, add it to @tops2. This is the same format as @tops: a list of query numbers/headers without ">" symbols. To do: make @tops a hash instead and delete headers that don't pass, instead of making @tops2? Depends whether @tops2 can be a hash.
					}
			} else {
				print "ERROR: FOUND A BLAST ENTRY WITH READ LENGTH 0!\n";
				print $log_filehandle "ERROR: FOUND A BLAST ENTRY WITH READ LENGTH 0!\n";
				next;
				} # If there isn't actually a sequence, move on.
						
	}

	$/ = "\n"; # Return input record separator to default ("\n"). # TO DO: redundant?

	$number_of_headers = scalar(@tops2); # Count the elements in @tops2.
	print "$number_of_headers reads to process after filtering by coverage.\n\nProcessing BLAST hits...\n";
	print $log_filehandle "$number_of_headers reads to process after filtering by coverage.\n\nProcessing BLAST hits...\n";

	# This section checked the coverage of each header in @tops using $blastfile. Headers with sufficient coverage were saved in @tops2.

	
	# For reads that passed the coverage check...
	#--------------------------------------------	
	my $current_read = 1; # This counts the number of reads to process. It's just something for a human to look at, so starts at 1 instead of 0.
	foreach my $header (@tops2) {
		print "\t$current_read of $number_of_headers: $header\n";
		print $log_filehandle "\t$current_read of $number_of_headers: $header\n";
		$current_read ++;
		
		
		# Get a list of sequences producing significant alignments (the BLAST hits) and their E values
		#---------------------------------------------------------------------------------------------
        # To do: combine with check coverage for each read above?
		# First, get the full BLAST entry for this read.
		open (my $blast_filehandle, $blastfile) or die "Cannot open $blastfile\n";
		my $line_contains_header = 0; # A flag. Starts at 0 but actually defaults to 1 once running.
		my $match = "Query= ".$header; # Put the current header back into BLAST format.
		my @blast_entry; # This array will store the BLAST entry for the current read.
		
		while (my $line = <$blast_filehandle>) { # Look at every line of $blastfile. This code is complicated: the second if statement is actually the first one to get activated.
					
					# To do: add a preliminary index that asks whether the first character is Q, or something similar. If it is, search the whole line for $header. Check whether this is faster than just searching every line once for $header. Is it better to have more searches but where most are smaller, or fewer but where all are bigger?
					
					if ($line_contains_header == 1) { # Statement 3! If $line_contains_header is 1, the *previous* $line contained the header. If this line is part of the same BLAST entry, we want it. If it's the start of a new entry, we've got what we need and can exit.
						if ($line =~/Query= /) {
							$line_contains_header = 0;  close OUT; }
					}
					if (index($line, $match) != -1) { # Statement 1! Double negative. Search for the header in this $line of the BLAST file. If there is a match, the output isn't -1, so continue.
						$line_contains_header = 1; # Switch the flag on.
					} # But there is no else: don't switch the flag off, even if $line doesn't contain the header. Only statement 3 switches off.
					if ($line_contains_header == 1) { # Statement 2! Save the current line to @blast_entry.
						push (@blast_entry, $line);
					}
		}
		close $blast_filehandle; # Do close at some point, or it will be kept open. Could make future searching easier if in order? But what happens at the end of this block? To do.
		
		my @blast_hits = (); # Now we will shorten @blast_entry to include just lines under "Sequences producing significant alignments". These are the BLAST hits. To do: could combine this with the above block, but it could get complicated.
		my $line_number = 6; # Start at line 6 (counting from 0, remember). This line always contains the 0th BLAST hit.
		my $next = 0;
		do {
					my $line = $blast_entry[$line_number];
					if ($line eq "\n") {
						$next = 1; # The line after the end of "Sequences producing significant alignments" is just a newline. If we've reach that, exit the loop.
					} else {
						my @line = (); @line = split (/\s\s+/, $line); # Split $line on two or more spaces. This separates the columns.
						my $line_length = scalar @line; # $line_length is the number of elements (columns) in @line.
						my $e_value = pop @line; # The E values are in the final column.
						chomp $e_value;
						push (@blast_hits, join(" ", $line[1], $e_value)); # Save two parts of $line in @blast_hits: the description and the E value.
						$line_number++; # Ready to look at the next in the list of "Sequences producing significant alignments".
					}
		} until ($next == "1");
		
		
		# Find the Identities score (coverage) and E value for the top BLAST hit
		#-----------------------------------------------------------------------
		my $found = 0; my $expect = (); my $identities = (); # $expect is the E (expect) value for the top BLAST hit. $identities is the number of bases in the top BLAST hit that match the query. It's in the format "86/92 (93%)".
		HIT:  foreach my $line (@blast_entry) { # We're looking at every line in the BLAST entry again.
				if ($found == 1) { next HIT;} # If we've already dealt with this BLAST hit, move on to the next.
				if ($line =~/Score = /) { # Search for the line containing "Score = ".
						my @line = (); @line = split (/\=/, $line); # Split the line on "=". The E value is the last thing on the line, so everything after the "=" is relevant.
						$expect = pop (@line); # Save that last chunk as $expect. It will be something like " 3e-29".
						chomp $expect; # Take the newline off the end.
				}
				if ($line =~/Identities = /) { # Search for the line containing "Identities = ".	
						my @line = (); @line = split (/\,/, $line); # The line looks like this: "Identities = 86/92 (93%), Gaps = 0/92 (0%)". Split it on that comma.
						my $identities_section = $line[0];
						my @identities_section = (); @identities_section = split (/\= /, $identities_section); # Now the gaps section is out of the way, split on "=" to isolate the value for identities.
						$identities = $identities_section[1];
						$found = 1; # Mark this BLAST hit as done.
				}
		}
	
		
		# Find the taxon for each BLAST hit using the names index file
		#-------------------------------------------------------------
		# But only $cap unique taxa for each read! $cap is 100 by default.
		my $number_of_blast_hits = scalar @blast_hits; # This is NOT the number of unique BLAST taxa.
		my $number_of_blast_taxa = 0; # *This* is the number of unique BLAST taxa.
        my @blast_info = ();
		my $lastname = ();
		my %hit_names = (); # Contains unique BLAST hit names. Only used to check for unique hits. Unique hits are then processed using the full $blast_hit.
		my @scores = (); # Contains unique BLAST E values. To do: lots of questions about @scores.
        
		BLASTIT:    foreach my $blast_hit (@blast_hits) {
                
				if ($number_of_blast_taxa == $cap) {
					next BLASTIT; # We will only look at the first $cap unique BLAST taxa.
				}
				
				my @blast_hit = split(/ /, $blast_hit); # Split the BLAST hit line on single spaces.
                my $current_hit_ID = 0;
                my $current_e_value = pop @blast_hit; # The E value is always the final element.
                my $original_hit_name = join (' ', @blast_hit); # Join the whole description field together. It should contain the organism name plus extra fluff.
                my $current_hit_name = $original_hit_name;
                
                # Make $current_hit_name more closely resemble names in the names file:
                if (index ($current_hit_name, ',') != -1) { # If $hit contains a comma,
                    my @current_hit_name = split(',', $current_hit_name); # Remove anything after a comma. Sometimes BLAST hits go like "Sphingobium cloacae DNA, complete genome, strain: JCM...", but only two names in the names file (and they're synonyms) contain ", strain", so what comes after a comma is overwhelmingly useless.
                    $current_hit_name = $current_hit_name[0];
                }
                if (index ($current_hit_name, 'virus') != -1) { # If $hit contains 'virus',
                    my @current_hit_name = split('virus', $current_hit_name); # Remove anything after "virus". The species name of the virus is like "Ovine lentivirus" or "Influenza A virus". Any detail beyond that is classed in the nodes file not as subspecies etc., but as "no rank", making it useless to the PIA. So don't worry about matching more specifically than species.
                    $current_hit_name = $current_hit_name[0] . 'virus'; # Keep the pre-'virus' part but add 'virus' back on (split removed it).
                }
                my @current_hit_name = split(' ', $current_hit_name); # Split into words.
                if ($current_hit_name[0] eq "PREDICTED:") { # Remove the word "PREDICTED:" from the start of the hit. We don't care if it's predicted or not.
                    shift @current_hit_name;
                }
                if ($current_hit_name[0] eq 'Uncultured') { # If the description starts with "Uncultured", change to "uncultured". The names file always has it lower-case.
                    $current_hit_name[0] = 'uncultured';
                }
                $current_hit_name = join (' ', @current_hit_name);
                
                ## TESTING: Faster to chop out these with extra hash lookups per name?
                #my %non_org_words;
                #undef $non_org_words{'strain'}; # There are many names containing "strain" in the names file but, like viruses, the ranks all seem to be "no rank" and therefore useless to the PIA. So trim down to species.
                #undef $non_org_words{'chloroplast'};
                #undef $non_org_words{'clone'};
                #undef $non_org_words{'culture-collection'};
                #undef $non_org_words{'genome'};
                #undef $non_org_words{'genomic'};
                #undef $non_org_words{'hypothetical'};
                #undef $non_org_words{'large'};
                #undef $non_org_words{'partial'};
                #undef $non_org_words{'probable'};
                #undef $non_org_words{'protein'};
                #undef $non_org_words{'ribosomal'};
                #
                #my @current_hit_name2;
                #foreach my $word (@current_hit_name) {
                #    if (exists $non_org_words{$word}) {
                #        last;
                #    } else {
                #        push (@current_hit_name2, $word);
                #    }
                #}
                ##print join (" ", @current_hit_name2);
                #@current_hit_name = @current_hit_name2;
                
				if ( exists $hit_names{$current_hit_name} ) { # Check whether we've already seen this tidied name.
					next BLASTIT; # If $current_hit_name is already listed in %hit_names, we've already seen it, so move on.
					}
				if ( grep( /^$current_e_value$/, @scores ) ) {
					#print "Next- Drop-off not reached\n"; # Not my note; what does drop-off mean?
					#print $log_filehandle "Next- Drop-off not reached\n";
					next BLASTIT; #  If $current_e_value is already listed in @scores, move on. Uh, not sure why. To do. Are we assuming that if the taxon is identical, the score will be identical and vice versa? Is this right? @scores is only used for this purpose. What does drop-off mean?
				}
				push @scores, $current_e_value; # Note the current E value in @scores. It only contains unique E values, but I don't know why they have to be unique. To do.
				$hit_names{$current_hit_name} = undef; # Note the current hit name in %hit_names. It only contains unique names.
                
                @current_hit_name = split(' ', $current_hit_name); # Split into words again.

                    
                # Get the taxon ID from the names file.	
                my %namesfileDBMnames = (); # Set up a hash to hold the names DBM file.
                tie (%namesfileDBMnames, "DB_File", $namesfileDBMnames) or die "Can't open $namesfileDBMnames:$!\n";
                
                my $found = 0;
                while (@current_hit_name) { # Do the following until last is called or you run out of words in $current_hit_name.
                    my $current_hit_name = join (' ', @current_hit_name);
                    #print "Match attempt: $current_hit_name\n";
                         
                    if (exists $namesfileDBMnames{$current_hit_name}) {
                        # To do: this does not account for names across classification systems. It does not account for taxonomic trickiness.
                        my @namesfile_info = split ("\t", $namesfileDBMnames{$current_hit_name});
                        $current_hit_ID = $namesfile_info[0]; # Save out the ID. We'll look up the rest of the information shortly.
                        $found = 1;
                        last; # Once a match is found, break out of the loop.
                    } else {
                        pop @current_hit_name; # If no match, take a word off the end of $current_hit_name and try again.
                    }
                }
                if ($found == 0) {
                    print "\t\tFailed to match BLAST hit '$original_hit_name' to names file!\n";
                    print $log_filehandle "\t\tFailed to match BLAST hit '$original_hit_name' to names file!\n";
                    my $current_hit_info = 0 . "\t" . 'undef' . "\t" . 'unassigned'; # Store the default null values.
                    push (@blast_info, $current_hit_info);
                    next BLASTIT;
                }
				$number_of_blast_taxa++; # Add one unique BLAST taxon to the count.				
                untie %namesfileDBMnames;
                
                # Look up the more information about this hit using its ID.
                $current_hit_name = retrieve_name ($current_hit_ID, $namesfileDBMids);
                my $current_hit_rank = retrieve_rank ($current_hit_ID, $nodesfileDBM);
                my $current_hit_info = $current_hit_ID . "\t" . $current_hit_name . "\t" . $current_hit_rank;
                
                push (@blast_info, $current_hit_info); # Each element contains [ID\tname\trank] for one BLAST hit.
		}
        
        #print "Number of unique BLAST taxa for this read: $number_of_blast_taxa\n";
		#print Dumper \@blast_info;

		# This section searched the names file DBM for information about each BLAST taxon.
		# However, it only considered up to $cap unique BLAST taxa. Duplicates and any additional unique taxa are ignored.
        # Names that could not be matched to the names file are carried forward with null values, so still count for diversity. Unfortunately, if they are the top, second-top or bottom hit the intersection cannot be calculated.
		
		
        # We have a list of unique BLAST taxa. How different are they?
        #-------------------------------------------------------------
        my $tax_diversity = scalar @blast_info; # $tax_diversity is the number of unique IDs among these BLAST taxa.
        my @contrastinghit_info = ();  # $contrastinghit is second unique name after the top hit: number 1 if you're counting from 0.
        my $contrastinghit_name = (); 
		if ($tax_diversity > 1) { # If the BLAST taxa aren't all the same:
                @contrastinghit_info = split ("\t", $blast_info[1]);
                $contrastinghit_name = $contrastinghit_info[1]; # Fetch the contrasting hit name from its info array.
		} else { 
				$contrastinghit_name = "no hit"; # If the BLAST taxa are all the same, there isn't a contrasting hit.
		}

		# Calculate the taxonomic diversity score.
		my $tax_diversity_score = ($tax_diversity/$cap) - (1/$cap); # Remember, $tax_diversity is the number of unique BLAST taxa and $cap defaults to 100.
		
        # So, are the top and second BLAST taxa in the same genus? Family? The search down from each taxon stops at class in most cases. If they aren't at least in the same class, it's not a worthy intersection? To do: sounds a little arbitrary.
		my $intersect = (); my $intersect_rank = (); my $intersect_name = "none found"; my $tophit_route = ();
        my @tophit_info = split ("\t", $blast_info[0]);
		if ($tax_diversity == 1) { # If the BLAST taxa are all the same:
					$intersect_rank = "none"; $intersect = 0; # We can't calculate an intersect.
		} else {
                my $tophit_ID = $tophit_info[0];
				if ($tophit_ID == 0){ # If the BLAST taxa aren't all the same but the top taxon has an ID 0, which means a proper ID was never found and it stayed on the default of 0:
					$intersect_rank = "none"; $intersect_name = "none found"; # We also can't calculate an intersect.
				} else { # If the BLAST taxa aren't all the same and the top one is valid:
						$tophit_route = retrieve_taxonomic_structure ($tophit_ID, $nodesfileDBM); # retrieve_taxonomic_structure() returns the route from this taxon down to either the root (if the taxon rank was higher than class) or down to its parent class (if the taxon rank was class or lower). To do: not sure why there's this difference.
						my $contrastinghit_ID = $contrastinghit_info[0]; # $contrastinghit_ID is the ID of the second-best taxon.
						my $contrastinghit_route = retrieve_taxonomic_structure ($contrastinghit_ID, $nodesfileDBM);
                        
						my $intersect = find_taxonomic_intersect ($tophit_route, $contrastinghit_route); # find_taxonomic_intersect() returns the lowest shared rank between the two routes. If one or more routes were undefined, it returns 0.
						if ($intersect == 0) {
								$intersect_rank = "none"; $intersect_name = "none found";
						} else { # If there was an intersect, find its rank and name.
								$intersect_rank = retrieve_rank ($intersect, $nodesfileDBM);
								$intersect_name = retrieve_name ($intersect, $namesfileDBMids);
						}	
				}
		}
        
		# How robust is the intersection?
		#--------------------------------
		# To help evaluate the robustness of the first intersect, ascertain the spread of blast hits phylogenetically by taking the bottom hit below cap and getting that intersect: the top intersect.
		my $topintersect = (); my $topintersect_rank = "none found"; my $bottomhit_name = "none"; my $topintersect_name = "none found"; my $bottomhit_ID = 0;
		if ($tax_diversity == 1) { $topintersect = "none"; $topintersect_rank = "none found";} # Remember, $tax_diversity is the number of unique BLAST taxa. If there's only one, there are no intersections at all.
		if ($tax_diversity == 2) { $topintersect = $intersect_rank; $topintersect_rank = $intersect_name; $bottomhit_name = $contrastinghit_name; # If there were only two BLAST taxa, the top intersect is the same as the intersect.
		}
		if ($tax_diversity > 2) { # If there were more than two BLAST taxa, it gets more complicated.
				my $bottomhit = pop @blast_info; # $bottom is the final unique BLAST taxon ID. The least good BLAST match.
                my @bottomhit_info = split ("\t", $bottomhit);
                my $bottomhit_ID = $bottomhit_info[0];
				my $bottomhit_name = $bottomhit_info[1];
				
				if ($bottomhit_ID==0){
					$topintersect = 0; # If the bottom ID is 0, the taxon is invalid, so we can't calculate a top intersect.
				} else {
					my $bottomhit_route = retrieve_taxonomic_structure ($bottomhit_ID, $nodesfileDBM); # retrieve_taxonomic_structure() returns the route from this taxon to either its parent class or the root.
					$topintersect = find_taxonomic_intersect ($tophit_route, $bottomhit_route); # Returns the lowest common taxon for the top and bottom BLAST taxa.
					
				}
				if ($topintersect == 0) {
					$topintersect_rank = "none"; $topintersect_name = "none found";
				} else {
					$topintersect_rank = retrieve_rank ($topintersect, $nodesfileDBM); $topintersect_name = retrieve_name ($topintersect, $namesfileDBMids);
					# $topintersectionclass is the rank of the lowest common taxon between the top and bottom BLAST taxa. $topranking is the name.
				}
		}
        
		# Print all of this information to an intersects.txt file
		#--------------------------------------------------------
        my $tophit_name = $tophit_info[1];
		open (INTERSECTS, ">>".$corename."/"."$corename".".intersects.txt") or die "Cannot write intersects file ".$corename."\n"; # Open intersect file for appending.

		chomp $header;
		print INTERSECTS "Query: $header, first hit: $tophit_name, expect: $expect, identities: $identities, next hit: $contrastinghit_name, last hit: $bottomhit_name, most distant classification intersect: $topintersect_rank, phylogenetic range: $topintersect_name, number of hits: $number_of_blast_hits, taxa diversity score: $tax_diversity_score, classification intersect: $intersect_name, id confidence class: $intersect_rank\n";
		# Close output file
		close INTERSECTS;

	}
    
	return $corename; # Return the $corename path for use in other subroutines.
}


sub retrieve_name {
##### Use taxid to get classification
	my ($query_ID, $namesfileDBMids) = @_; # $query_ID is the taxonomic ID in question. We want to find its name.
	my $name = undef;	# We want the name for this ID. Default to undef.
    
    unless ($query_ID == 0) { # ID 0 is not in the names file.
        my %namesfileDBMids = (); # Set up a fresh hash to hold the names DBM file.
        tie (%namesfileDBMids, "DB_File", $namesfileDBMids) or die "Can't open $namesfileDBMids: $!\n";
        
        if (exists $namesfileDBMids{$query_ID}) {
            $name = $namesfileDBMids{$query_ID};
        } else {
            print "ERROR: ID $query_ID is not 0 but was not found in names file. Assuming no name.\n";
        }
    }
    untie %namesfileDBMids;
	return $name; 
}


sub find_taxonomic_intersect {
##### Compare both taxonomic routes, find intersect
	my ($first_route, $second_route) = @_;
	my $intersect = 0; # If either the first or second route aren't defined, return the intersect as 0.
	if (defined $first_route && defined $second_route){
		my @first_route = (); my @second_route = ();
		@first_route = split (/\t/, $first_route);
		@second_route = split (/\t/, $second_route);
		my $connected = 0; # The $connected flag activates when the first shared rank is stored as $intercect and prevents it from being overwritten. Because the routes go from lower to higher taxa, higher shared ranks are ignored. To do: there's probably a more efficient way to achieve this.
		foreach my $rank (@first_route) { # Start with each rank in the first route.
			foreach my $otherrank (@second_route) { # Take the ranks in the second route.
				if ($rank == $otherrank) { # If the ranks match,
					unless ($connected == 1) { # (If the routes are already flagged as connected, move on to the next rank in the second route.)
						$connected = 1; # Flag them as connected.
						$intersect = $rank; # $interesect becomes this shared rank.
					}
				}
			}
		}
	}
	return $intersect;	
}


sub retrieve_rank {
##### Retrieve rank of taxonomic ID from nodes DBM file
	my ($query_ID, $nodesfileDBM) = @_; # $query_ID is the taxonomic ID of the name in question. $nodesfileDBM is the path to the nodes index file.
    my $rank = 'unassigned'; # The taxonomic rank, like "family" or "genus". Defaults to 'unassigned'.
    
    unless ($query_ID == 0) { # ID 0 is not in the nodes file.
        my %nnodesfileDBM = (); # Set up a fresh hash to hold the nodes DBM file.
        tie (%nodesfileDBM, "DB_File", $nodesfileDBM) or die "Can't open $nodesfileDBM: $!\n";
        
        if (exists $nodesfileDBM{$query_ID}) {
            my @node_info = split("\t", $nodesfileDBM{$query_ID});
            $rank = $node_info[1]; # The 0th field is parent node. We don't need that here.
        } else {
            print "ERROR: ID $query_ID is not 0 but was not found in nodes file. Rank unassigned.\n";
        }
    }
    untie %namesfileDBMids;
	return $rank; 
}


sub retrieve_taxonomic_structure {
##### Get hierarchy from nodes.dmp file
	my ($query_ID, $nodesfileDBM) = @_; # $query_ID is a non-0 taxonomic ID. $nodesfileDBM is the path to the nodes index file.
    my $route = undef; # If there's a problem, default the route to undef.
    
    unless ($query_ID == 0) { # ID 0 is not in the nodes file.
    
        my $exit = 0; # When to exit the do loop below.
        my $next_level_ID; # The ID of the parent node.
        my $rank; # The rank of the current ID.
        my @route = (); # @route is a list of tab-separated taxonomic IDs moving down from the current node. @route goes down to the root if the first taxon is ranked higher than class. If the first taxon is a class or lower, @route only goes as far as class. To do: why?
        my %nnodesfileDBM = (); # Set up a fresh hash to hold the nodes DBM file.
        tie (%nodesfileDBM, "DB_File", $nodesfileDBM) or die "Can't open $nodesfileDBM: $!\n";
        
        do {
            push (@route, $query_ID); # Add the current ID to @route.
    
            if (exists $nodesfileDBM{$query_ID}) {
                    my @node_info = split("\t", $nodesfileDBM{$query_ID});
                    $next_level_ID = $node_info[0];
                    $rank = $node_info[1];
            } 
            untie %namesfileDBMids;
        
            if ($query_ID == $next_level_ID) { # If the current node is its parent, we're at the root. We have a route.
                $exit = 1;
            } elsif ($rank =~/class/) { # If the current node is ranked class, that's far enough. We have a route.
                $exit = 1;
            }
            $query_ID = $next_level_ID; # If we're not yet at the root or a class, move on to the current parent node.
            
        } until ($exit == 1);
        
        $route = join ("\t", @route);
    }
    return $route;
}
