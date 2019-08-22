#!/usr/bin/perl

######################################################
##################### PIA_inner.pl ###################
######################################################
#########  Phylogenetic Intersection Analysis ########
############## Robin Allaby, UoW 2013 ################
######################################################
############## Version 4.7, 2019-08-22 ###############
######################################################

# Edited by Roselyn Ware, UoW 2015
# A method of metagenomic phylogenetic assignation, designed to be robust to partial representation of organisms in the database		
		
# Further edited by Becky Cribdon, UoW 2019.
# - Names search only returns scientific names.
# - Added a log.
# - Made variable names more consistent and meaningful.
# - Removed $blastfile2 and functions that were not fully implemented.
# - Added index files for the names and nodes files.
# - Names search accounts for duplicate names. It only takes potential names in the expected phylogenetic range, then chooses the one with the highest rank (to be conservative). If multiple names have the highest rank, it chooses from them randomly.

# Edits for indexing re-vamp August 2019:
# - Removed extended summary.
# - Removed $namesfile as an argument for PIA().
# - Removed indexing. Now in PIA_index_maker.pl.
# - Added ', O_RDONLY, 0666, $DB_BTREE' to all tie calls. The DBM files are now in DB_BTREE format instead of DB_HASH.
# - Removed names and nodes input options. Just go straight to the DBM files in References/. The README is clear that they need to be in there. Not finding them via the original names and nodes files means that those can be deleted.
# - Added the section "# Make copies of the DBM index files for this run of PIA_inner.pl to use".
# - Un-commented the unlink to remove the copies at the end.

# Issues:
# - Think about class thing. Why stop at class?
# - At the moment, it accepts things like "Zea mays chloroplast" and "Zea mays full genome" as different BLAST taxa. Should probably take just "Zea mays" instead.
# - About excluding BLAST hits which have E-values that have already been seen: instead of excluding 'repeats', how about averaging them by finding their intersection and treating that as a read? That's what Smith 2015 implies is happening. HOWEVER, taxonomic diversity should take into account the full list. It's a measure of the database, not the read.
# - Re-name "taxa diversity" and "taxa diversity score" in the intersects file to "taxon" or "taxonomic". But I'll have to update the FASTA making script too.
# - Anything else labelled "to do".
# Please report any problems to r.cribdon@warwick.ac.uk.

# To run, first create a header file from the final FASTA:
# > cat test_sample.fasta | grep ">" > test_sample.header
# Then run PIA_inner.pl giving the header file as -f and a BLAST file as -b. -p is the parent taxon you expect your reads to fall into. It's used to help choose the right taxon when a name matches more than one.
# > perl -f [FASTA or corresponding header file] -b [BLAST file] -p [taxonomic ID of expected phylogenetic range] 
# The header file lists reads by their query number. The BLAST file contains every read (by query number) and a list of reasonable matches (BLAST hits). If you use the PIA.pl wrapper, it makes header files for you.


	use strict;
	use warnings;
	use lib './Modules';
	use FileChecks; # A bespoke module in /Modules.
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
	getopts('f:b:c:C:hp:', \%options); 														#Getopt::Std
    
	# If other text found on command line, do:
	print "Other things found on the command line:\n" if $ARGV[0];
	foreach (@ARGV)	{
        print "$_\n";
	}			
	
##### Check header filename and open file #####
	# The header file is simply the read names (headers) extracted from a FASTA file. PIA.pl makes the header file.
	# Ensure that a filename has been given as an argument. -f flag
	my $tophitfile = FileChecks::process_filename($options{f});									#FileChecks.pm
	chomp $tophitfile ;
	
##### Check BLAST filename and open file #####
	# Ensure that a filename has been given as an argument. -b flag
 	my $blastfile = FileChecks::process_filename($options{b});									#FileChecks.pm
	chomp $blastfile ;

##### See if cap input #####
	# The cap impacts the taxon diversity score. Default is 100.	
	# Ensure that a number has been given as an argument. -c flag
	my $capopt = ($options{c});
	my $cap = 100;
	if ($capopt) { $cap = $capopt;} # If there is an option, overwrite the default.

##### See if min % coverage input #####
	# The minimum percentage coverage a top BLAST hit must have for a read to be taken forward. Default is 95.
	# Ensure that a number has been given as an argument. -C flag
	my $min_coverage_perc_opt = ($options{C});
	my $min_coverage_perc = 95;
	if ($min_coverage_perc_opt) { $min_coverage_perc = $min_coverage_perc_opt;} # If there is an option, overwrite the default.

##### Display help file and exit if -h flag called #####
	my $helpfile="Helpfile_PIA.txt";
	if ($options{h}){
		FileChecks::process_help($helpfile);														#FileChecks.pm	
	}

##### See if expected phylogenetic range input #####	
	# Ensure that a number has been given as an argument. -p flag
	my $expected_phylogenetic_range_opt = ($options{p});
	my $expected_phylogenetic_range = 1;
	if ($expected_phylogenetic_range_opt) {$expected_phylogenetic_range = $expected_phylogenetic_range_opt;}


######################################################
######### Make copies of the DBM index files #########
######################################################
    # TO DO: replace with file locking eventually
    my $nodesfileDBM = 'Reference_files/nodes.dmp.dbm' . "_$tophitfile";
    system ("cp Reference_files/nodes.dmp.dbm $nodesfileDBM");
    my %nodesfileDBM = ();
    
    my $namesfileDBMnames = 'Reference_files/names.dmp_names.dbm' . "_$tophitfile";
    system ("cp Reference_files/names.dmp_names.dbm $namesfileDBMnames");
    my %namesfileDBMnames = (); # And initialise a hash to hold the DBM later.
    
    my $namesfileDBMids = 'Reference_files/names.dmp_IDs.dbm' . "_$tophitfile";
    system ("cp Reference_files/names.dmp_IDs.dbm $namesfileDBMids");
    my %namesfileDBMids = ();
    
    my $abbreviationsDBM = 'Reference_files/Species_abbreviations.txt.dbm' . "_$tophitfile  ";
    system ("cp Reference_files/Species_abbreviations.txt.dbm $abbreviationsDBM");
    my %abbreviationsDBM = ();



######################################
########### Start log file ###########
######################################

    my $log_filename = $tophitfile . '_PIA_inner_log.txt';
    open( my $log_filehandle, '>', $log_filename) or die "Cannot open $log_filename for writing.\n$!\n"; # Note: this will overwrite old logs. I didn't see any point in appending old ones when you need to remove the other PIA outputs before running it again anyway.
    use IO::Handle; # Enable autoflush.
    $log_filehandle -> autoflush(1); # Set autoflush to 1 for the log filehandle. This means that Perl won't buffer its output to the log, so the log should be updated in real time.
    print $log_filehandle "#####################################\n";
    
    # Print run parameters to log.
	my $min_taxdiv_score = 0.1; # The minimum taxonomic diversity score is 0.1. Reads must have at least eleven unique taxa. This is not an option because it depends on $cap: they are related by maths that I don't know. But it's helpful to see what they both are. 
    
    print $log_filehandle "\n****Inputs****\nHeader file:\t$tophitfile\nBLAST file:\t$blastfile\nExpected phylogenetic range (ID):\t$expected_phylogenetic_range\nMinimum coverage for top BLAST hit:\t$min_coverage_perc %\nCap of BLAST taxa to examine (default is 100):\t$cap\nMinimum taxonomic diversity score:\t$min_taxdiv_score\n\n";
    

######################################################
###################### Run PIA #######################
######################################################
    
	print "\n****Starting first round of PIA****\n";
	print $log_filehandle "\n\n****Starting first round of PIA****\n";
	
	my ($corename) = PIA($tophitfile, $blastfile, $cap, $min_coverage_perc); # PIA() returns a base name for this sample file. The base name is [header file]_out.
	print $log_filehandle "\nPIA() subroutine finished.\n\n";
	print "\nPIA() subroutine finished.\n\n";
    
	#my $corename = '50.header_out'; # IF NOT ACTUALLY RUNNING THE PIA; FOR TESTING


######################################################	
################## Summarise Data ####################
######################################################

##### Extract simple summary from the intersects file
	my $tophitfile2 = "$corename"."/"."$corename.intersects.txt";
	my $simplesummary = simple_summary($tophitfile2, $min_taxdiv_score);

###### Extract extended summary (optional) - simple summary with lineage data
#	if ($summary){
#		extended_summary($simplesummary, $nodesfile, $namesfile);
#	}


######################################################
##################### Tidy Up ########################
######################################################

##### Remove un-necessary files	
	unlink("$corename"."/"."temp_blast_entry.txt");
	unlink("$corename"."/"."TEMP");
	unlink("$corename"."/"."hittempfile");
    unlink $namesfileDBMnames;
    unlink $namesfileDBMids;
    unlink $nodesfileDBM;
    unlink $abbreviationsDBM;
	
# Finish the log and move it into the output directory.
    print $log_filehandle "****This run of PIA_inner.pl is finished.****\n\n\n\n";
    close $log_filehandle;
    my $log_final_destination = $corename . '/' . $corename . '_PIA_inner_log.txt';
    system("mv $log_filename $log_final_destination");
    
    system('echo "This run of PIA_inner.pl is finished.\n"');
	
	
######################################################	
######################################################
#################### SUBROUTINES #####################
######################################################
######################################################	

#sub extended_summary{
###### Take simple summary and append lineage information 
#	# Optional= only if -e flag supplied
#	my ($simplesummary, $nodesfile, $namesfile) = @_;
#	open (SUM, $simplesummary) or die "Cannot write simplesummary file\n$!\n";
#	my @original=<SUM>;
#	close SUM;
#	my @simplesum=();
#	@simplesum= split (/\//, $simplesummary);
#	my $name= $simplesum[0];
#	my @lineagesarrayfinal= TreeOfLife::get_lineage_array(\@original,$nodesfile,$namesfile);	#TreeOfLife.pm
#	my $output =$name."/";
#	#$name=$name."_Extended_Summary";
#	TreeOfLife::format_PIA_output($name,$nodesfile,$namesfile, \@lineagesarrayfinal,$output);				#TreeOfLife.pm
#
#}


sub simple_summary {
#### Create simple summary of output- This is the default. An extended summary can be created if -e flag is supplied
	my ($PIAinputname, $min_taxdiv_score) = @_; # $PIAinputname is $tophitfile2: basically the header file.
    
	my @intersects_file= FileChecks::open_file($PIAinputname);										#FileChecks.pm					
	my %intersects = (); # Keys are intersect names and IDs in the format "name (ID)". Values are the number of times that name (ID) occurs.
    
    # Get a list of classification intersects where the taxa diversity score was at least $min_taxdiv_score.
	foreach my $line (@intersects_file){
		my @row= split(/, classification intersect: |, id confidence class: /, $line); # Split on the classification intersect field first (note that it won't match to "most distant classification intersect"), followed by the ID confidence field. This is not an 'or'. It's one after the other, chopping off text from the left and right sides to leave just the classification intersect value in the middle.
		my @check= split(/, taxa diversity score: |, classification intersect: /, $line); # Similarly, leave just the taxa diversity score.
        
        if ($check[1] >= $min_taxdiv_score ){ # $check[1] is the taxa diversity score. 'none found' means there was no intersect.
            
            # Change the format of the name field ($row[1]) for outputting.
            my $name_field = $row[1];
            my @name_field = split (/ /, $name_field);
            my $ID = pop @name_field; # The ID is the last word.
            chomp $ID;
            $ID =~ tr/()//d; # Remove the parentheses from it (this is transliterate with delete).
            $name_field = join (" ", @name_field); # Join the remaining words back together. These are the taxon name.
            my $name_and_ID = $name_field . "\t" . $ID; # Join the name and the ID with a tab.
            
            if (exists $intersects{$name_and_ID}) {
                $intersects{$name_and_ID} = $intersects{$name_and_ID} + 1;
            } else {
                $intersects{$name_and_ID} = 1;
            }
        }
	}
    
	my @name=split ("\/",$PIAinputname); # Pick out a sample name from $PIAinputname to use in the output file.
	my $name=$name[0];
    
    my $summary_basic_filename = $name."_Summary_Basic.txt";
	open (OUT, ">", $PIAinputname . "_Summary_Basic.txt") or die "Cannot write output file\n$!\n";
    #print "PIA input name: $PIAinputname\n";
	print OUT "#Series:\t$name\n"; # Output $name as a header.

    foreach my $intersect (keys %intersects) {
        unless ($intersect eq "none found\t0") {
            print OUT $intersect . "\t" . $intersects{$intersect} . "\n";
        }
    }

	close OUT;
	return $summary_basic_filename;
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
    
	my ($tophitfile, $blastfile, $cap, $min_coverage_perc) = @_;
	
	my $corename = $tophitfile . "_out"; # # Generate a core name: a base name for naming things linked to this sample. For example, "test_sample.header_out".
	`mkdir $corename`; # Create a separate directory for outputs (makes folder organisation tidier).
    
    
	# List the header of each read
	#-----------------------------
	open (TOPS, $tophitfile) or die "Cannot open $tophitfile\n!$\n"; # Extract headers from the header file. A header is the ID number given to a read, like a name, so a header represents a read. "Header" and "read" are used somewhat interchangably in this code, but it never deals directly with a read. It always works via the header.

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
			my $read_length_section = `grep -A 5 "$header2" "$blastfile"`; # Search in $blastfile for $header2. Return the two lines after it: a newline and the query length line.
            if (index($read_length_section, '***** No hits found *****') != -1) { next; } # If BLAST did not find any hits, skip this read. The warning is printed just under length.
			
            # Coverage is [match length] / [read length]. First, find the read length.
			my @read_length_section0 = split(/Length=/, $read_length_section); # Split $read_length_section on Length= followed by whitespace. Everything in between is the value for read length.
            my @read_length_section = split(/\s/, $read_length_section0[1]);
            my $read_length = $read_length_section[0];
            
			# Then, if there is in fact a sequence, find the match length. This is the length of the top BLAST hit and is also called 'query cover':  how much of the query sequence is covered by the top hit (not necessarily covered perfectly). Identities (number of identical bases) are given as a fraction of coverage. If the identities are 86/92, then coverage is 92 bp, but there are 6 mismatches.
			# So, like BLAST, the PIA defines "coverage" as the total length covered regardless of whether there are mismatches or not. This is fine.
			if ($read_length > 0 ){ # Just in case there's a read length of 0... which I'm not sure is possible.
					my $blast_entry_partial=`sed -n -e '/$header/,/Identities =/ p' "$blastfile"`;
					# sed is used as a selective print. -n means "no print" and p here means "print if matches". It will print all lines beginning with the first pattern ("$header") and ending with the second pattern ("Identities =") regardless of what goes on in between.
					# So, this code searches $blastfile for the BLAST entry of $header, but only prints up to the end of the line containing "Identities =".
					my @blast_entry_partial= split(/Identities = [0-9]*\//, $blast_entry_partial); # Split the partial BLAST entry at "Identities =".
					my $match_length_section = $blast_entry_partial[1]; # Pick out what comes after "Identities =", so the end of the final line. It will be something like "63 (95%), Gaps = 0/63 (0%)".
                    
					if ( length( $match_length_section || '' )) { # Uh. || basically means "or". So, if there is anything in $match_length_section? I have no idea. To do: understand.
							my @match_length_section = split(/ \([0-9]*\%\), Gaps/, $match_length_section); # Split $match_length_section into the values for identity and gaps.
							my $match_length = $match_length_section[0]; # Pick out just the numerical value of the match length. It's a single number.
                            
                            # Calculate percentage coverage and filter out those less than or equal to the minimum (user defined ARGV[2] or default of 100). We calculate percentage ourselves because the percentage in the BLAST output is rounded.
							my $coverage = $match_length / $read_length;
							my $min_coverage = $min_coverage_perc / 100;
							if ($coverage >= $min_coverage) {
									push @tops2, $header;
							} # If this BLAST entry, this query number, has sufficient coverage, add it to @tops2. This is the same format as @tops: a list of query numbers/headers without ">" symbols. To do: make @tops a hash instead and delete headers that don't pass, instead of making @tops2? Depends whether @tops2 can be a hash.
					}
			} else {next;} # If there isn't actually a sequence, move on.
						
	}

	$/ = "\n"; # Return input record separator to default ("\n"). # TO DO: redundant?

	$number_of_headers = scalar(@tops2); # Count the elements in @tops2.
	print "$number_of_headers reads to process after checking for BLAST hits and filtering by coverage.\n\nProcessing BLAST hits...\n";
	print $log_filehandle "$number_of_headers reads to process after checking for BLAST hits and filtering by coverage.\n\nProcessing BLAST hits...\n";

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
		open (my $blast_filehandle, $blastfile) or die "Cannot open $blastfile\n$!\n";
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
	
		
		# Find the taxonomic ID for each BLAST hit using the names index file
		#--------------------------------------------------------------------
		# But only $cap unique taxa for each read! $cap is 100 by default.
		my $number_of_blast_hits = scalar @blast_hits; # This is NOT the number of unique BLAST taxa.
		my $number_of_blast_taxa = 0; # *This* is the number of unique BLAST taxa.
        my @blast_info = ();
		my $lastname = ();
		my %hit_names = (); # Contains unique BLAST hit names. Only used to check for unique hits. Later hit processing uses the full $blast_hit.
        my @scores = (); # Contains unique BLAST E values. To do: lots of questions about @scores.
        
        
		BLASTHIT:    foreach my $blast_hit (@blast_hits) {

				if ($number_of_blast_taxa == $cap) {
					next BLASTHIT; # We will only look at the first $cap unique BLAST taxa.
				}
                
				my @blast_hit = split(/ /, $blast_hit); # Split the BLAST hit line on single spaces.
                my $current_hit_ID = 0;
                my $current_e_value = pop @blast_hit; # The E value is always the final element.
                my $original_hit = join (' ', @blast_hit); # Join the whole description field together. It should contain the organism name plus extra fluff.
                
                
                # Clean up the name to make it more closely resemble names in the names file
                #---------------------------------------------------------------------------
                my $clean_hit = $original_hit;
                if (index ($clean_hit, ',') != -1) { # If $hit contains a comma,
                    my @clean_hit = split(',', $clean_hit); # Remove anything after a comma. Sometimes BLAST hits go like "Sphingobium cloacae DNA, complete genome, strain: JCM...", but only two names in the names file (and they're synonyms) contain ", strain", so what comes after a comma is overwhelmingly useless.
                    $clean_hit = $clean_hit[0];
                }
                
                #if (index ($clean_hit, 'virus') != -1) { # If $clean_hit contains 'virus',      # TO DO: is this really wise? Forcing every virus name to species or higher?
                #    my @clean_hit = split('virus', $clean_hit); # Remove anything after "virus". The species name of the virus is like "Ovine lentivirus" or "Influenza A virus". Any detail beyond that is classed in the nodes file not as subspecies etc., but as "no rank", making it useless to the PIA. So don't worry about matching more specifically than species.
                #    $clean_hit = $clean_hit[0] . 'virus'; # Keep the pre-'virus' part but add 'virus' back on (split removed it).
                #}
                
                my @clean_hit1 = split(' ', $clean_hit); # Split into words.
                my @clean_hit2 = (); # An intermediate for further cleaning.
                foreach my $word (@clean_hit1) {
                    if ($word eq 'UNVERIFIED:' or $word eq 'PREDICTED:' or $word eq 'TPA:' or $word eq 'TPA_asm:') { # Remove these words. They are about the sequence, not the organism. TPA stands for third party annotation, by the way.
                        next;
                    }
                    push (@clean_hit2, $word);
                }
                $clean_hit = join (' ', @clean_hit2);
                
                
                # Check whether this hit or E-value has already been seen
                #---------------------------------------------------------
				if ( exists $hit_names{$clean_hit} ) { # If $clean_hit_name is already listed in %hit_names, we've already seen it, so move on.
					next BLASTHIT; 
				}
                if ( grep( /^$current_e_value$/, @scores ) ) {
					#print "Next- Drop-off not reached\n"; # To do: not my note; what does drop-off mean?
					#print $log_filehandle "Next- Drop-off not reached\n";
					next BLASTHIT; #  If $current_e_value is already listed in @scores, move on. Uh, not sure why. To do. Are we assuming that if the taxon is identical, the score will be identical and vice versa? Is this right?
				}
                
                $number_of_blast_taxa++; # Not already seen? Add one unique BLAST taxon to the count.
                push @scores, $current_e_value; # Note the current E value in @scores. It only contains unique E values. To do.
				$hit_names{$clean_hit} = undef; # Note the hit name in %hit_names. It only contains unique hit names.
                
                
                # Search for a matching taxon
                #----------------------------
				$found = 0; # Reset the found flag.
                
                # First, two hard-coded difficult names that won't be matched otherwise.
                if (index ($clean_hit, 'Genomic sequence from Human') != -1 ) { # Badly-named human sequences.
                    my $current_hit_info = 9606 . "\t" . 'Homo sapiens' . "\t" . 'species';
                    push (@blast_info, $current_hit_info);
                    next BLASTHIT;
                }
                
                if (index ($clean_hit, 'MACACA MULATTA BAC clone') != -1 || index ($clean_hit, 'Rhesus Macaque BAC') != -1 || index ($clean_hit, 'Rhesus Macaque Centromeric ') != -1 ) { # There are 400 NCBI sequences for the Rhesus macaque where the organism name is in all-caps and an unknown number where the species epithet is capitalised.
                    my $current_hit_info = 9544 . "\t" . 'Macaca mulatta' . "\t" . 'species';
                    push (@blast_info, $current_hit_info);
                    next BLASTHIT;
                }
                if ($original_hit eq 'F.odoratum large subunit ribosomal RNA') { # This exact hit is to Flavobacterium odoratum (a synonym). Could be confused with Funastrum odoratum. But since it's a synonym, it's unlikely to be submitted again, so I've just hard-coded it.
                    my $current_hit_info = 256 . "\t" . 'Myroides odoratus' . "\t" . 'species';
                    push (@blast_info, $current_hit_info);
                    next BLASTHIT;
                }
                if ($original_hit eq 'complete chromosome Acholeplasma palmae') { # This exact hit is to Acholeplasma palmae.
                    my $current_hit_info = 38986 . "\t" . 'Acholeplasma palmae' . "\t" . 'species';
                    push (@blast_info, $current_hit_info);
                    next BLASTHIT;
                }
                
                
                my $current_hit = $clean_hit; # Keep a copy of $clean_hit and save a working copy as $current_hit.
                my @current_hit = split(' ', $current_hit); # Split the name into words again.
                
                # A small number of hits start with gramatically-incorrect species abbreviations. These are those I've found. Check for them before starting the huge hash lookups.
                if ($found == 0) { # A small number of hits start with gramatically-incorrect species abbreviations. These are those I've found. Check for them before starting the huge hash lookups.
                    tie (%abbreviationsDBM, "DB_File", $abbreviationsDBM, O_RDONLY, 0666, $DB_BTREE) or die "Can't open $abbreviationsDBM: $!\n";
                    
                    if (exists $abbreviationsDBM{$current_hit[0]}) {
                                $current_hit_ID = $abbreviationsDBM{$current_hit[0]};
                                $found = 1;
                                print "\t\tFound '$original_hit' in the abbreviations DBM and assigned to $current_hit_ID\n";
                    }
                    untie %abbreviationsDBM;
                }
                
                
                my $not_in_expected_range = 0; # A flag for when the name does have at least one match in the names file, but none in the expected phylogenetic range.
                if ($found == 0) {
                    ($current_hit_ID, $not_in_expected_range) = search_names_file($namesfileDBMnames, $nodesfileDBM, $log_filehandle, $original_hit, $expected_phylogenetic_range, \@current_hit);
    
                    if ($not_in_expected_range == 1) { # If there were matches but not in the right range, make the hit information null and move on to the next BLAST hit.
                        print "\t\tBLAST hit '$original_hit' matched names file, but no matches in the expected phylogenetic range. Setting values to null.\n";
                        print $log_filehandle "\t\tBLAST hit '$original_hit' matched names file, but no matches in the expected phylogenetic range. Setting values to null.\n";
                        my $current_hit_info = 0 . "\t" . 'undef' . "\t" . 'unassigned';
                        push (@blast_info, $current_hit_info);
                        next BLASTHIT;
                    }
                    
                    if ($current_hit_ID != 0) { # If there was a good match, mark as found.
                        $found = 1;
                    }
                }
                
                
                if ($found == 0) { # We've tried all hard-coded options and we've tried iterating word by word. Now try uncapitalising the first word and iterating again.
                    my $first_word = $current_hit[0];
                    my @first_word = split ('', $first_word);
                    my $first_letter = $first_word[0];
                    # Is the first letter upper-case?
                    if ($first_word[0] =~ /^\p{Uppercase}/) { # If the first letter of the first word is indeed a capital,
                        $first_word[0] = lc $first_word[0]; # Make it lower case.
                        $first_word = join ('', @first_word); # Join the word back together.
                        $current_hit[0] = $first_word; # Overwrite the first word of @current_hit.
                        
                        # And run the iterative names file search again:
                        ($current_hit_ID, $not_in_expected_range) = search_names_file($namesfileDBMnames, $nodesfileDBM, $log_filehandle, $original_hit, $expected_phylogenetic_range, \@current_hit);
                        
                        if ($not_in_expected_range == 1) { # If there were matches but not in the right range, make the hit information null and move on to the next BLAST hit.
                            print "\t\tBLAST hit '$original_hit' matched names file after making the first word lower-case, but no matches in the expected phylogenetic range. Setting values to null.\n";
                            print $log_filehandle "\t\tBLAST hit '$original_hit' matched names file after making the first word lower-case, but no matches in the expected phylogenetic range. Setting values to null.\n";
                            my $current_hit_info = 0 . "\t" . 'undef' . "\t" . 'unassigned';
                            push (@blast_info, $current_hit_info);
                            next BLASTHIT;
                        }
                        
                        if ($current_hit_ID != 0) { # If there was a good match, mark as found.
                            $found = 1;
                        }
                    
                    } # If the first word was lower case already, never mind.
                }
                    
                
                if ($found == 0) { # If the name still hasn't been matched to an ID, make the hit information null and move on to the next BLAST hit.
                    print "\t\tFailed to match BLAST hit '$original_hit' to names file!\n";
                    print $log_filehandle "\t\tFailed to match BLAST hit '$original_hit' to names file!\n";
                    my $current_hit_info = 0 . "\t" . 'undef' . "\t" . 'unassigned';
                    push (@blast_info, $current_hit_info);
                    next BLASTHIT;
                }
					
                untie %namesfileDBMnames;
                
                # Look up the more information about this hit using its ID.
                my $current_hit_name = retrieve_name ($current_hit_ID, $namesfileDBMids);
                my $current_hit_rank = retrieve_rank ($current_hit_ID, $nodesfileDBM);
                #print "\t\tMatched '$original_hit_name' to '$current_hit_name ($current_hit_ID, $current_hit_rank)'\n"; # TESTING
                #print $log_filehandle "\t\tMatched '$original_hit_name' to '$current_hit_name ($current_hit_ID, $current_hit_rank)'\n"; # TESTING
                my $current_hit_info = $current_hit_ID . "\t" . $current_hit_name . "\t" . $current_hit_rank;
                
                push (@blast_info, $current_hit_info); # Each element contains [ID\tname\trank] for one BLAST hit.
        } # End of BLASTHIT.
        
        #print "Number of unique BLAST taxa for this read: $number_of_blast_taxa\n";
        #print Dumper \@blast_info;

		# This section searched the names file DBM for information about each BLAST taxon.
		# However, it only considered up to $cap unique BLAST taxa. Duplicates and any additional unique taxa are ignored.
        # Names that could not be matched to the names file are carried forward with null values, so still count for diversity. Unfortunately, if they are the top, second-top or bottom hit the intersection cannot be calculated.
		
		
        # We have a list of unique BLAST taxa. How different are they?
        #-------------------------------------------------------------
        my $contrastinghit_ID = 0; my $contrastinghit_name = "none found"; # contrastinghit is second BLAST taxon: number 1 if you're counting from 0. Default the values to null.
		if ($number_of_blast_taxa > 1) { # If the BLAST taxa aren't all the same:
                my @contrastinghit_info = split ("\t", $blast_info[1]);
                $contrastinghit_ID = $contrastinghit_info[0]; # Fetch the contrasting hit ID and name from its info array.
                $contrastinghit_name = $contrastinghit_info[1]; 
		}
		# Calculate the taxonomic diversity score.
		my $tax_diversity_score = ($number_of_blast_taxa/$cap) - (1/$cap); # Remember,  $cap defaults to 100.
        #print "\t\tTaxa diversity = $tax_diversity\tScore = $tax_diversity_score\n";
        #print $log_filehandle "\t\tTaxa diversity = $tax_diversity\tScore = $tax_diversity_score\n";
		
        # So, are the top and second BLAST taxa in the same genus? Family? The search down from each taxon stops at class in most cases. If they aren't at least in the same class, it's not a worthy intersection? To do: sounds a little arbitrary.
        my $intersect_ID = 0; my $intersect_rank = "none"; my $intersect_name = "none found"; my $tophit_ID = 0; my $tophit_route = 0; # Default the intersect values to null.
        my @tophit_info = split ("\t", $blast_info[0]); # tophit is the first BLAST taxon. It's the top BLAST hit.

        if ($number_of_blast_taxa > 1) { # If the BLAST taxa are not all the same:
                $tophit_ID = $tophit_info[0];

				unless ($tophit_ID == 0) { # The BLAST taxa might not all be the same, but if the top taxon has ID 0, a proper ID was never found and we can't calculate an intersect.
                		$tophit_route = retrieve_taxonomic_structure ($tophit_ID, $nodesfileDBM, 1); # retrieve_taxonomic_structure() returns the route from this taxon down to either the root (if the taxon rank was higher than class) or down to its parent class (if the taxon rank was class or lower). To do: not sure why there's this difference.
	
						my $contrastinghit_route = retrieve_taxonomic_structure ($contrastinghit_ID, $nodesfileDBM, 1);
                        
						$intersect_ID = find_taxonomic_intersect ($tophit_route, $contrastinghit_route); # find_taxonomic_intersect() returns the lowest shared rank between the two routes. If there wasn't any shared rank or if one or more routes were undefined, it returns an ID of 0.
						if ($intersect_ID == 0) {
								$intersect_rank = "none"; $intersect_name = "none found";
						} else { # If there was an intersect, find its rank and name.
								$intersect_rank = retrieve_rank ($intersect_ID, $nodesfileDBM);
								$intersect_name = retrieve_name ($intersect_ID, $namesfileDBMids);
						}	
				}
		}
        
        my $tophit_name = 'noname';
        if ($tophit_info[1]) {
            $tophit_name = $tophit_info[1];
        }

        
		# How robust is the intersection?
		#--------------------------------
		# To help evaluate the robustness of the first intersect, ascertain the spread of blast hits phylogenetically by taking the bottom hit (within $cap) and getting that intersect: the top intersect. Note that these things aren't used to calculate or evaluate anything, but they are in the output.
		my $topintersect_ID = 0; my $topintersect_rank = "none found"; my $topintersect_name = "none found"; my $bottomhit_ID = 0; my $bottomhit_name = "none found"; # Default to null values.
		if ($number_of_blast_taxa == 2) { # If there were only two BLAST taxa, the top intersect is the same as the intersect and the bottom hit is the same as the contrasting hit.
            $topintersect_ID = $intersect_ID; $topintersect_rank = $intersect_rank; $topintersect_name = $intersect_name;
            $bottomhit_ID = $contrastinghit_ID; $bottomhit_name = $contrastinghit_name; 
		}
		if ($number_of_blast_taxa > 2) { # If there were more than two BLAST taxa, it gets more complicated.
				my $bottomhit = pop @blast_info; # $bottom is the final unique BLAST taxon ID. The least good BLAST match.
                my @bottomhit_info = split ("\t", $bottomhit);
                $bottomhit_ID = $bottomhit_info[0];
				$bottomhit_name = $bottomhit_info[1];

				unless ($bottomhit_ID == 0) { # If the bottom hit doesn't have an ID, we can't calculate a top intersect.
					my $bottomhit_route = retrieve_taxonomic_structure ($bottomhit_ID, $nodesfileDBM, 1); # retrieve_taxonomic_structure() returns the route from this taxon to either its parent class or the root.
					$topintersect_ID = find_taxonomic_intersect ($tophit_route, $bottomhit_route); # Returns the lowest common taxon for the top and bottom BLAST taxa.	
				}
                
				if ($topintersect_ID == 0) {			
                    $topintersect_name = "none found"; #$topintersect_rank = "none";
				} else {
                    $topintersect_name = retrieve_name ($topintersect_ID, $namesfileDBMids); #$topintersect_rank = retrieve_rank ($topintersect_ID, $nodesfileDBM);
				}
		}
        
        
		# Print all of this information to an intersects.txt file
		#--------------------------------------------------------
		open (INTERSECTS, ">>".$corename."/"."$corename".".intersects.txt") or die "Cannot write intersects file ".$corename.".intersects.txt\n$!\n"; # Open intersect file for appending.

		chomp $header;
		print INTERSECTS "Query: $header, first hit: $tophit_name ($tophit_ID), expect: $expect, identities: $identities, next hit: $contrastinghit_name ($contrastinghit_ID), last hit: $bottomhit_name ($bottomhit_ID), phylogenetic range of hits up to cap: $topintersect_name ($topintersect_ID), number of hits: $number_of_blast_hits, taxa diversity: $number_of_blast_taxa, taxa diversity score: $tax_diversity_score, classification intersect: $intersect_name ($intersect_ID)\n";
		# Note that there was originally "most distant classification intersect" ($topintersect_rank) after "last hit", but I don't see why the rank of $topintersect_name needs to go in the output, especially as all things topintersect are only there for interest. They aren't actually used to calculate or evaluate anything.
        # Similarly, I removed "id confidence class ($intersect_rank) from the end.
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
        tie (%namesfileDBMids, "DB_File", $namesfileDBMids, O_RDONLY, 0666, $DB_BTREE) or die "Can't open $namesfileDBMids: $!\n";
        
        if (exists $namesfileDBMids{$query_ID}) {
            $name = $namesfileDBMids{$query_ID};
        } else {
            print "\t\tERROR: ID $query_ID is not 0 but was not found in names file. Assuming no name.\n";
            print $log_filehandle "\t\tERROR: ID $query_ID is not 0 but was not found in names file. Assuming no name.\n";
        }
    }
    untie %namesfileDBMids;
	return $name; 
}


sub find_taxonomic_intersect {
##### Compare both taxonomic routes, find intersect
	my ($first_route, $second_route) = @_; # Routes are the taxonomic IDs of successive parent taxa.
	my $intersect = 0; # If either the first or second route aren't defined, return the intersect as ID 0.
	if (defined $first_route && defined $second_route){
		my @first_route = (); my @second_route = ();
		@first_route = split (/\t/, $first_route);
		@second_route = split (/\t/, $second_route);
		my $connected = 0; # The $connected flag activates when the first shared rank is stored as $intersect and prevents it from being overwritten. Because the routes go from lower to higher taxa, higher shared ranks are ignored.
		foreach my $rank (@first_route) { # Start with each rank in the first route.
			foreach my $otherrank (@second_route) { # Take the ranks in the second route.
				if ($rank == $otherrank) { # If the ranks match,
					unless ($connected == 1) { # (If the routes are already flagged as connected, move on to the next rank in the second route.)
						$connected = 1; # Flag them as connected.
						$intersect = $rank; # $intersect becomes this shared rank.
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
        my %nodesfileDBM = (); # Set up a fresh hash to hold the nodes DBM file.
        tie (%nodesfileDBM, "DB_File", $nodesfileDBM, O_RDONLY, 0666, $DB_BTREE) or die "Can't open $nodesfileDBM: $!\n";
        
        if (exists $nodesfileDBM{$query_ID}) {
            my @node_info = split("\t", $nodesfileDBM{$query_ID});
            $rank = $node_info[1]; # The 0th field is parent node. We don't need that here.
        } else {
            print "ERROR: ID $query_ID is not 0 but was not found in nodes file. Rank unassigned.\n";
        }
        untie %nodesfileDBM;
    }
	return $rank; 
}


sub retrieve_taxonomic_structure {
##### Get hierarchy from nodes.dmp file
	my ($query_ID, $nodesfileDBM, $stop_at_class) = @_; # $query_ID is a non-0 taxonomic ID. $nodesfileDBM is the path to the nodes index file. $stop_at_class is a 0/1 switch asking whether the route should stop at class (enough for the main PIA) or not (while checking names).
    my $route = undef; # If there's a problem, default the route to undef.
    
    unless ($query_ID == 0) { # ID 0 is not in the nodes file.
    
        my $exit = 0; # When to exit the do loop below.
        my $next_level_ID; # The ID of the parent node.
        my $rank; # The rank of the current ID.
        my @route = (); # @route is a list of tab-separated taxonomic IDs moving down from the current node. @route goes down to the root if the first taxon is ranked higher than class. If the first taxon is a class or lower, @route only goes as far as class. To do: why?
        my %nodesfileDBM = (); # Set up a fresh hash to hold the nodes DBM file.
        tie (%nodesfileDBM, "DB_File", $nodesfileDBM, O_RDONLY, 0666, $DB_BTREE) or die "Can't open $nodesfileDBM: $!\n";
        
        do {
            push (@route, $query_ID); # Add the current ID to @route.
    
            if (exists $nodesfileDBM{$query_ID}) {
                    my @node_info = split("\t", $nodesfileDBM{$query_ID});
                    $next_level_ID = $node_info[0];
                    $rank = $node_info[1];
            } 
            
            if ($query_ID == $next_level_ID) { # If the current node is its parent, we're at the root. We have a route.
                $exit = 1;
            } elsif ($stop_at_class == 1 and $rank =~/class/) { # If the current node is ranked class, that's far enough. We have a route.
                $exit = 1;
            }
            $query_ID = $next_level_ID; # If we're not yet at the root or a class, move on to the current parent node.
            
        } until ($exit == 1);
        
        untie %nodesfileDBM;
        $route = join ("\t", @route);
    }
    return $route;
}


sub search_names_file {
##### Try to match the a cleaned BLAST name to the names file
my ($namesfileDBMnames, $nodesfileDBM, $log_filehandle, $original_hit_name, $expected_phylogenetic_range, $current_hit_name_ref) = @_;
my @current_hit_name = @$current_hit_name_ref;

my %namesfileDBMnames = (); # Set up a hash to hold the names DBM file.
    tie (%namesfileDBMnames, "DB_File", $namesfileDBMnames, O_RDONLY, 0666, $DB_BTREE) or die "Can't open $namesfileDBMnames: $!\n";
    my $current_hit_ID = 0; # Like in the main code body, default to 0.
    my $not_in_expected_range = 0; # A flag for when the name does have at least one match in the names file, but none in the expected phylogenetic range.

    while (@current_hit_name) { # Do the following until last is called or you run out of words in $current_hit_name.
                    my $current_hit_name = join (' ', @current_hit_name);
                    #print "Match attempt: $current_hit_name\n";
                         
                    if (exists $namesfileDBMnames{$current_hit_name}) { # If the BLAST name has a match in the names file,
                        my @namesfile_matches = split ("\t", $namesfileDBMnames{$current_hit_name}); # See if there is more than one ID under this name.
                        if (exists $namesfile_matches[1]) { # If there's more than one ID, things gonna get complicated.
                                my @potential_taxa_in_right_range = (); # If there's more than one ID in the expected phylogenetic range, things gonna get even more complicated.
                                print "\t\tMultiple potential matches for '$original_hit_name':";
                                print $log_filehandle "\t\tMultiple potential matches for '$original_hit_name':";

                                foreach my $potential_taxon (@namesfile_matches) { # For each potential taxon, pull out its clarification and see whether it is within the expected phylogenetic range. For example, if an ID is in Insecta but the reads should all be in Viridiplantae, the ID is not the right one.
                                        print "\t$potential_taxon";
                                        print $log_filehandle "\t$potential_taxon";
   
                                        if ($potential_taxon == $expected_phylogenetic_range) { # If it's the same as the range, accept it.
                                                push (@potential_taxa_in_right_range, $potential_taxon);
                                        } else {
                                                my @potential_parent_taxa = split("\t", retrieve_taxonomic_structure($potential_taxon, $nodesfileDBM, 0) ); # If not, check the parent IDs.
                                                @potential_parent_taxa = reverse @potential_parent_taxa; # The expected phylogenetic range is probably a high taxon, so start with the higher parent taxa.
                                                foreach my $potential_parent_taxon (@potential_parent_taxa) {
                                                    if ($potential_parent_taxon == $expected_phylogenetic_range) {
                                                        push (@potential_taxa_in_right_range, $potential_taxon);
                                                        
                                                    }
                                                }
                                        }
                                }
                               
                                # If any potential taxa were in the expected phylogenetic range, choose the highest. We can't distinguish between them any further, so assign cautiously. Note that if multiple taxa have the joint highest rank, it will pick randomly because they are stored in a hash.
                                if (@potential_taxa_in_right_range) {
                                    my %potential_taxa_in_right_range_ranks = ();
                                    foreach my $potential_taxon (@potential_taxa_in_right_range) { # Find the taxonomic rank for each one.
                                        $potential_taxa_in_right_range_ranks{retrieve_rank($potential_taxon, $nodesfileDBM)} = $potential_taxon;
                                    }
                                    
                                    # Pick the taxon with the highest rank.
                                    my @ranks = ('superkingdom', 'kingdom', 'phylum', 'subphylum', 'class', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'family', 'genus', 'species', 'subspecies');
                                    foreach my $rank (@ranks) {
                                        if (exists $potential_taxa_in_right_range_ranks{$rank}) {
                                            $current_hit_ID = $potential_taxa_in_right_range_ranks{$rank};
                                            print "\n\t\t\tAssigned '$original_hit_name' to $current_hit_ID\n";
                                            return ($current_hit_ID, $not_in_expected_range); # We've assigned a taxon to the BLAST name, so exit the search.
                                        }
                                    }
                                }
                                
                                print "\n";
                                print $log_filehandle "\n";
                                
                                $not_in_expected_range = 1; # If there were no potential taxa in the expected phylogenetic range, flag and exit subroutine.
                                return ($current_hit_ID, $not_in_expected_range);

                        } else { # If there's just one ID for this name, things are a lot simpler.
                                $current_hit_ID = $namesfile_matches[0]; # Save the ID. We'll use this to look up more information shortly.
                                return ($current_hit_ID, $not_in_expected_range); # Exit the subroutine.
                        }
                    } else {
                        pop @current_hit_name; # If no match, take a word off the end of $current_hit_name and try again.
                    }
    } # End of while loop. If a name makes it this far, it wasn't found in the names file and its ID will remain at 0.
    return ($current_hit_ID, $not_in_expected_range);
}