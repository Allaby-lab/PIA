#!/usr/bin/perl

######################################################
##################### PIA_inner.pl ###################
######################################################
#########  Phylogenetic Intersection Analysis ########
############## Robin Allaby, UoW 2013 ################
######################################################
############## Version 4.7, 2019-08-23 ###############
######################################################

# NOTE: 4.7 is the same as 4.6 but without changes 2 and 3, which need to be discussed

# Edited by Roselyn Ware, UoW 2015
# A method of metagenomic phylogenetic assignation, designed to be robust to partial representation of organisms in the database		
		
# Further edited by Becky Cribdon, UoW 2019.
# - Names search only returns scientific names.
# - Added a log.
# - Made variable names more consistent and meaningful.
# - Removed $blastfile2 and functions that were not fully implemented.
# - Uses pre-made DBM index files for nodes.dmp and names.dmp
# - Names search accounts for duplicate names. It only takes potential names in the expected phylogenetic range, then chooses the one with the highest rank (to be conservative). If multiple names have the highest rank, it chooses from them randomly.
# - BLAST hits that can't be identified are ignored.
# - Only one BLAST hit is looked at per hit score. But instead of just taking the first one, it now takes an intersection of all hits with that hit score, and the taxonomic diversity of all of those hits still counts towards the diversity score. Note that if $cap was already being reached, $cap will now be reached earlier (with closer BLAST hits).
# - Intersections are taken wherever possible. No stopping at class.

# Issues:
# - Parent and child taxa, like "Zea" and "Zea mays", are taken as two different taxa. This isn't strictly right, but I'm not sure how to fix it.
# - Re-name "taxa diversity" and "taxa diversity score" in the intersects file to "taxon" or "taxonomic". But I'll have to update the FASTA making script too.
# - Incorporate the dodgy chimpanzee names and others.
# - Change the column order in Summary Basic to ID, name, count.
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
	use Data::Dumper qw(Dumper);
    use DB_File;
    use Fcntl;


######################################################										#####################
########### Check arguments and Input Data ###########										###### Modules ######
######################################################										#####################

##### Get arguments from command line #####
	my %options=();
	getopts('f:b:c:C:hm:p:s:', \%options); 														#Getopt::Std
    
	# If other text found on command line, do:
	print "Other things found on the command line:\n" if $ARGV[0];
	foreach (@ARGV)	{
        print "$_\n";
	}			
	
##### Check header filename and open file #####
	# The header file is the read names (headers) extracted from a FASTA file. PIA.pl makes the header file.
	my $header_filename = FileChecks::process_filename($options{f});									#FileChecks.pm
	
##### Check BLAST filename and open file #####
 	my $blast_filename = FileChecks::process_filename($options{b});									#FileChecks.pm

##### See if cap input #####
	# BLAST hits are all assigned to taxa. $cap is the maximum number of taxa to be looked at. If $cap is hit, no more BLAST hits will be considered. $cap is used to calculate the taxonomic diversity score, so affects the minimum score threshold (option s). Default is 100.	
	my $cap_opt = ($options{c});
	my $cap = 100;
	if ($cap_opt) { $cap = $cap_opt;} # If there is an option, overwrite the default.

##### See if min % coverage input #####
	# The minimum percentage coverage a top BLAST hit must have for a read to be taken forward. Default is 95.
	my $min_coverage_perc_opt = ($options{C});
	my $min_coverage_perc = 95;
	if ($min_coverage_perc_opt) { $min_coverage_perc = $min_coverage_perc_opt;} # If there is an option, overwrite the default.

##### Display help file and exit if -h flag called #####
	my $helpfile="Helpfile_PIA.txt";
	if ($options{h}){
		FileChecks::process_help($helpfile);														#FileChecks.pm	
	}
    
##### See if expected phylogenetic range input #####	
	my $expected_phylogenetic_range_opt = ($options{p});
	my $expected_phylogenetic_range = 1;
	if ($expected_phylogenetic_range_opt) {$expected_phylogenetic_range = $expected_phylogenetic_range_opt;}

##### See if min taxonomic diversity score input #####
	# The minimum taxonomic diversity score a read must have to make it to Summary_Basic.txt. Depends on $cap. Defaults to 0.1.
	my $min_taxdiv_score_opt = ($options{s});
	my $min_taxdiv_score = 0.1;
	if ($min_taxdiv_score_opt) { $min_taxdiv_score = $min_taxdiv_score_opt;} # If there is an option, overwrite the default.


######################################
########### Start log file ###########
######################################

    my $log_filename = $header_filename . '_PIA_inner_log.txt';
    open( my $log_filehandle, '>', $log_filename) or die "Cannot open $log_filename for writing.\n$!\n"; # This will overwrite old logs.
    use IO::Handle; # Enable autoflush.
    $log_filehandle -> autoflush(1); # Set autoflush to 1 for the log filehandle. This means that Perl won't buffer its output to the log, so the log will be updated in real time.
    
    # Print run parameters to log.
    print $log_filehandle "#####################################\nOther things found on the command line:\n" if $ARGV[0];
    foreach (@ARGV)	{
        print $log_filehandle "$_\n";
    }    
    print $log_filehandle "\n****Inputs****\nInput header:\t$header_filename\nInput BLAST:\t$blast_filename\nMinimum coverage for top BLAST hit:\t$min_coverage_perc %\nCap of BLAST taxa to examine:\t$cap\nMinimum taxonomic diversity score:\t$min_taxdiv_score\nExpected phylogenetic range (ID):\t$expected_phylogenetic_range\n\n";


######################################################
######### Make copies of the DBM index files #########
######################################################
    # TO DO: replace with file locking eventually
    my $nodesfileDBM = 'Reference_files/nodes.dmp.dbm' . "_$header_filename";
    system ("cp Reference_files/nodes.dmp.dbm $nodesfileDBM");
    my %nodesfileDBM = ();
    
    my $namesfileDBMnames = 'Reference_files/names.dmp_names.dbm' . "_$header_filename";
    system ("cp Reference_files/names.dmp_names.dbm $namesfileDBMnames");
    my %namesfileDBMnames = (); # And initialise a hash to hold the DBM later.
    
    my $namesfileDBMids = 'Reference_files/names.dmp_IDs.dbm' . "_$header_filename";
    system ("cp Reference_files/names.dmp_IDs.dbm $namesfileDBMids");
    my %namesfileDBMids = ();
    
    my $abbreviationsDBM = 'Reference_files/Species_abbreviations.txt.dbm' . "_$header_filename";
    system ("cp Reference_files/Species_abbreviations.txt.dbm $abbreviationsDBM");
    my %abbreviationsDBM = ();


######################################################
###################### Run PIA #######################
######################################################
    
	print "\n****Starting first round of PIA****\n";
	print $log_filehandle "\n\n****Starting first round of PIA****\n";
	
	my $corename = PIA($header_filename, $blast_filename, $cap, $min_coverage_perc); # PIA() returns a base name for this sample file. The base name is [header file]_out.
	print $log_filehandle "\nPIA() subroutine finished.\n\n";
	print "\nPIA() subroutine finished.\n\n";
    
	#my $corename = '50.header_out'; # IF NOT ACTUALLY RUNNING THE PIA; FOR TESTING


######################################################	
################## Summarise Data ####################
######################################################

##### Extract simple summary from the intersects file
	my $header_filename2 = "$corename"."/"."$corename.intersects.txt";
	my $simplesummary = simple_summary($header_filename2, $min_taxdiv_score);


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
    print "This run of PIA_inner.pl is finished.\n";
    print $log_filehandle "****This run of PIA_inner.pl is finished.****\n\n\n\n";
    close $log_filehandle;
    my $log_final_destination = $corename . '/' . $corename . '_PIA_inner_log.txt';
    system("mv $log_filename $log_final_destination");



######################################################	
######################################################
#################### SUBROUTINES #####################
######################################################
######################################################


sub simple_summary {
#### Create simple summary of output. The intersects file is more informative, but this is what most people will take as the output.
	my ($header_filename2, $min_taxdiv_score) = @_;
    
	my @intersects_file= FileChecks::open_file($header_filename2);										#FileChecks.pm					
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
            my $ID_and_name = $ID . "\t" . $name_field; # Join the ID and namewith a tab.
            
            if (exists $intersects{$ID_and_name}) {
                $intersects{$ID_and_name} = $intersects{$ID_and_name} + 1;
            } else {
                $intersects{$ID_and_name} = 1;
            }
        }
	}
    
	my @name=split ("\/",$header_filename2); # Pick out a sample name from $PIAinputname to use in the output file.
	my $name=$name[0];
    
    my $summary_basic_filename = $name."_Summary_Basic.txt";
	open (my $summary_basic_filehandle, ">", $header_filename2 . "_Summary_Basic.txt") or die "Cannot write output file\n$!\n";
	print $summary_basic_filehandle "#Series:\t$name\n"; # Output $name as a header.

    foreach my $intersect (keys %intersects) {
        unless ($intersect eq "0\tnone found") {
            print $summary_basic_filehandle $intersect . "\t" . $intersects{$intersect} . "\n";
        }
    }

	close $summary_basic_filehandle;
	return $summary_basic_filename;
}


sub PIA {
##### Run phylogenetic intersection analysis

	# Summary
	#--------
	# List the header of each read
    # Look at the BLAST file entry by entry. If an entry matches a header:
    #   Skip if no BLAST hits
    #   Skip if insufficient coverage
    #   List BLAST hits
    #   Match BLAST hits to taxa
    #       Discard hits to taxa already seen
    #       Only consider up to $cap taxa
    #       Discard hits not matched to a taxon
    #       If hits share the same hit score, "average" them to their intersection
    #   Calculate the taxonomic diversity score for this read
    #   Find the intersection between the top and contrasting hits
    #   Find the intersection between the top and bottom hits
    #   Print all of this information to the intersects.txt file
	# Return $corename
    
	my ($header_filename, $blast_filename, $cap, $min_coverage_perc) = @_;
	
	my $corename = $header_filename . "_out"; # # Generate a core name: a base name for naming things linked to this sample. For example, "test_sample.header_out".
	`mkdir $corename`; # Create a separate directory for outputs (makes folder organisation tidier).
    
    
	# List the header of each read
	#-----------------------------
	open (my $header_filehandle, $header_filename) or die "Cannot open $header_filename\n!$\n"; # Extract headers from the header file. A header is the ID number given to a read, like a name, so a header represents a read. "Header" and "read" are used somewhat interchangably in this code, but it never deals directly with a read. It always works via the header.

	my %headers_all = (); # Will store all headers for this analysis.
	
	while (my $header = <$header_filehandle>) {
		substr ($header, 0, 1) = ""; # Remove the first character, which is the ">" symbol. These separate FASTA records. We don't need them.
		chomp $header; # Remove any trailing newlines.
        $headers_all{$header} = undef;		
	}
	close $header_filehandle;
	
	my $number_of_headers = scalar(keys %headers_all);
	print $log_filehandle "$number_of_headers reads to process.\n\n";
	

	# Search the BLAST file for entries for reads in @headers
	#--------------------------------------------------------
    # Entries for reads we want will then be analysed.    
    open (my $blast_filehandle, $blast_filename) or die "Cannot open $blast_filename\n"; # Open the BLAST file.
    $/ = 'Effective search space used: '; # Set the read separator to 'Effective search space used: ', which comes right near the end of every BLAST entry.
    
    #print Dumper @headers_all; print "\n\n";
    my $entry_count = 0;
    my $number_of_headers_to_find = $number_of_headers;
    my $current_header_number = 0;
    
    BLASTENTRY: while (my $entry = <$blast_filehandle>) { # Look at every entry in $blast_filename.
        $entry_count ++;
        
        # Isolate the query name (equivalent to header name).
        my @split_on_query = split(/\nQuery= |\n\nLength=/, $entry); # Split on the query field first followed by the length field. This is not an 'or'. It's one after the other, chopping off text from the left and right sides to leave just the query name in the middle.
        #print Dumper @split_on_query; print "\n\n";
        my $query;   
        $query = $split_on_query[1];

        if (exists $headers_all{$query}) { # If we've found the entry for a header, 
            delete($headers_all{$query}); # Delete the header from %headers_all
            $number_of_headers_to_find --;
            $current_header_number ++;
            print "\t$current_header_number of $number_of_headers: $query\n";
            print $log_filehandle "\t$current_header_number of $number_of_headers: $query\n";
            
            
            # First check: are there any BLAST hits at all?
            #----------------------------------------------
            if (index($split_on_query[2], "\n\n\n***** No hits found *****") != -1) { print "No hits for this entry\n"; next BLASTENTRY; } # $split_on_query[2] is the second and final part of the entry.
            
            
            # Second check: coverage
            #-----------------------
            # We calculate it ourselves because the percentage in the BLAST output is rounded. We need the read length and top-BLAST-hit match length. The latter is the second number in the Identities line.
            my @entry_remaining = split("\n", $split_on_query[2]); # We'll use this array later for listing BLAST hits too.
            my $read_length = $entry_remaining[0];
            
            my @split_on_Identities = split('Identities = ', $split_on_query[2]); # Split the partial BLAST entry at 'Identities ='.
            my @split_on_Gaps = split(', Gaps', $split_on_Identities[1]); # Split that to isolate just the Identities field.
            my $Identities = $split_on_Gaps[0]; # While we're at it, note the Identities field. This will be exported in the intersects file right at the end of PIA().

            my @Identities = split(/\/| \(/, $Identities); # The match length is the / value in Identities.
            my $match_length = $Identities[1];
            
            my $coverage = $match_length / $read_length;
  			my $min_coverage = $min_coverage_perc / 100;
			if ($coverage < $min_coverage) { next; } # If the top BLAST hit doesn't have at least $min_coverage, skip this read.
            
            
            # Note the E value of the top hit
            #--------------------------------
            my @split_on_E_value = split(",  Expect = ", $split_on_Identities[0]); # Split that to isolate just the E value.
            my $E_value = $split_on_E_value[1];
            $E_value =~ s/\s+$//g; # Trim whitespace off the  end. This will be exported in the intersects file right at the end of PIA().
            
            
            # Collect all of the BLAST hits
            #------------------------------
            my @blast_hits;
            foreach my $line (@entry_remaining[4..$#entry_remaining]) { # $entry_remaining[4] is the first BLAST hit.
                if ($line eq '') { last; } # The block finishes with an empty line.

                my @line = split (/\s\s+/, $line); # Otherwise, split $line on two or more spaces. This separates the columns.
				pop @line; # Discard the last element, the E value.
				push (@blast_hits, join(" ", $line[1], $line[-1])); # The last element is now the score. Save this and the description ([1]) in @blast_hits.
            }
            
            
            # Find the taxonomic ID for each BLAST hit using the names index file
            #====================================================================	
            my $number_of_blast_hits = scalar @blast_hits; # This is NOT the number of unique BLAST taxa. It's the number of hits to process for this read.
            my $number_of_blast_taxa = 0; # *This* is the number of unique BLAST taxa. It's a running count independent of %blast_taxa.
            my $lastname = ();
            my %blast_taxa = (); # A list of the unique BLAST taxa. Only used to check whether hits are to taxa that have already been seen. Later hit processing uses the full $blast_hit.
            my %scores = (); # A list of the unique hit scores.


    		BLASTHIT:    foreach my $blast_hit (@blast_hits) {

                if ($number_of_blast_taxa == $cap) { last BLASTHIT; } # But stop if hits to $cap taxa have already been found! $cap is 100 by default.
				
				my @blast_hit = split(/ /, $blast_hit); # Split the BLAST hit line on single spaces.
                my $current_hit_ID = 0;
                my $current_score = pop @blast_hit; # Since we processed the hit line above, the hit score is the last element.
                
                if (exists $scores{$current_score}) { next BLASTHIT; } # If we've already processed a hit with this score, move on to the next hit. (Edited in change 2)
                
                my $original_hit = join (' ', @blast_hit); # Every other element is a word in the description field. That's what we'll try to match to a name.
               
               
                # Clean up the name to make it more closely resemble names in the names file
                #---------------------------------------------------------------------------
                my $clean_hit = $original_hit;
                if (index ($clean_hit, ',') != -1) { # If $hit contains a comma,
                    my @clean_hit = split(',', $clean_hit); # Remove anything after a comma. A tiny number of names do contain commas, but they are low taxa.
                    # TO DO: Otherwise we'd have to remove commas after removing a word. Try this!
                    $clean_hit = $clean_hit[0];
                }
                
                my @clean_hit1 = split(' ', $clean_hit); # Split into words.
                
                $clean_hit1[-1] =~ s/\.{3}//; # Remove the ellipsis (I used to think sometimes just "..", but actually I've never seen that) from the last element, if there is one.
                $clean_hit1[-1] =~ s/\s$//g; # Also remove any trailing space that the last step might have left.
                
                my @clean_hit2 = (); # An intermediate for further cleaning.
                foreach my $word (@clean_hit1) {
                    if ($word eq 'UNVERIFIED:' or $word eq 'PREDICTED:' or $word eq 'TPA:' or $word eq 'TPA_asm:') { # Remove these words. They are about the sequence, not the organism. TPA stands for third party annotation, by the way.
                        next;
                    }
                    push (@clean_hit2, $word);
                }
                
                if (@clean_hit2[0..1]) { # If there are at least two words in $clean_hit
                    my $first_two_words = join (' ', @clean_hit2[0..1]);
                    if ($first_two_words eq 'Genomic sequence') {
                        @clean_hit2 = @clean_hit2[3..$#clean_hit2]; # If the first two words are "Genomic sequence", the hit will go "Genomic sequence for Zea mays" or "Genomic sequence to Zea mays" or something like that. So remove the first three words.
                    }
                }
                
                $clean_hit = join (' ', @clean_hit2);                
                
                
                # Search for a matching taxon name
                #---------------------------------
				my $found = 0;
                
                # First, some hard-coded difficult names that won't be matched otherwise.
                if (index ($clean_hit, 'MACACA MULATTA BAC clone') != -1 || index ($clean_hit, 'Rhesus Macaque BAC') != -1 || index ($clean_hit, 'Rhesus Macaque Centromeric ') != -1 ) { # There are 400 NCBI sequences for the Rhesus macaque where the organism name is in all-caps and an unknown number where the species epithet is capitalised.
                    $current_hit_ID = 9544;
                    $found = 1;
                }
                
                if ($found == 0 && $original_hit eq 'Wild cabbage satellite DNA') { # Brassica oleraceae. It's the single "wild cabbage" sequence on the NCBI (X07519.1).
                    $current_hit_ID =  3712;
                    $found = 1;
                }
                if ($found == 0 && $original_hit eq 'complete chromosome Acholeplasma palmae') { # Acholeplasma palmae
                    $current_hit_ID = 38986;
                    $found = 1;
                }
                if ( ($found == 0) && ($clean_hit eq 'strain 231/09 Stolbur phytoplasma draft' || $clean_hit eq 'strain 284/09 Stolbur phytoplasma draft') ) { # Candidatus Phytoplasma solani. These are the only two NCBI sequences in this format.
                    $current_hit_ID = 69896;
                    $found = 1;
                }
                if ($found == 0 && $original_hit eq 'Uncultured Aigarchaeota clone MB-2008-349-18 16S ribo...') { # 'uncultured thaumarchaeote'.
                    $current_hit_ID = 651141;
                    $found = 1;
                }
                if ($found == 0 && $original_hit eq 'Uncultured chloroflexi bacterium partial 16S rRNA gen...') { # 'uncultured Chloroflexi bacterium' with the genera capitalised.
                    $current_hit_ID = 166587;
                    $found = 1;
                }
                if ($found == 0 && $original_hit eq 'F.odoratum large subunit ribosomal RNA') { # This exact hit is to Flavobacterium odoratum (a synonym). Could be confused with Funastrum odoratum. But since it's a synonym, it's unlikely to be submitted again.
                    $current_hit_ID =  256;
                    $found = 1;
                }
                
                my $current_hit = $clean_hit; # If $clean_hit wasn't any of those, save a working copy of it as $current_hit.
                my @current_hit = split(' ', $current_hit); # Split the name into words again.

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
    
                    if ($not_in_expected_range == 1) { # If there were matches but not in the right range, move on to the next BLAST hit.
                        print "\t\tBLAST hit '$original_hit' matched names file, but no matches in the expected phylogenetic range.\n";
                        print $log_filehandle "\t\tBLAST hit '$original_hit' matched names file, but no matches in the expected phylogenetic range.\n";
                        next BLASTHIT; # Unidentified hits do not count towards anything.
                    }
                    
                    if ($current_hit_ID != 0) { # If there was a good match, mark as found.
                        $found = 1;
                    }
                }
                
                
                if ($found == 0) { # We've tried all hard-coded options and we've tried iterating word by word. Now try making the first word lower-case and iterating again.
                    my $first_word = $current_hit[0];
                    my @first_word = split ('', $first_word);
                    my $first_letter = $first_word[0];

                    if ($first_word[0] =~ /^\p{Uppercase}/) { # If the first letter of the first word is upper-case,
                        $first_word[0] = lc $first_word[0]; # Make it lower case.
                        $first_word = join ('', @first_word); # Join the word back together.
                        $current_hit[0] = $first_word; # Overwrite the first word of @current_hit.
                        
                        my $current_hit_temp = join (' ', @current_hit);
                        
                        ($current_hit_ID, $not_in_expected_range) = search_names_file($namesfileDBMnames, $nodesfileDBM, $log_filehandle, $original_hit, $expected_phylogenetic_range, \@current_hit); # And run the iterative names file search again.
                        
                        if ($not_in_expected_range == 1) { # If there were matches but not in the right range, move on to the next BLAST hit.
                            print "\t\tBLAST hit '$original_hit' matched names file, but no matches in the expected phylogenetic range.\n";
                            print $log_filehandle "\t\tBLAST hit '$original_hit' matched names file, but no matches in the expected phylogenetic range.\n";
                            next BLASTHIT; # Unidentified hits do not count towards anything.
                        }
                        
                        if ($current_hit_ID != 0) { # If there was a good match, mark as found.
                            $found = 1;
                        }
                    
                    } # If the first word was lower case already, never mind.
                }
                
                untie %namesfileDBMnames;
                
                
                if ($found == 0) { # If the name still hasn't been matched to an ID, move on to the next BLAST hit.
                    print "\t\tFailed to match BLAST hit '$original_hit' to names file!\n";
                    print $log_filehandle "\t\tFailed to match BLAST hit '$original_hit' to names file!\n";
                    next BLASTHIT; # Unidentified hits do not count towards anything.
                }
                
				if ( exists $blast_taxa{$current_hit_ID} ) { # Now we've identified the hit, have we already seen this taxon? If so, move on to the next hit.
					next BLASTHIT; 
				}
                $number_of_blast_taxa ++; # Not already seen? Add one unique BLAST taxon to the count.
				$blast_taxa{$current_hit_ID} = undef; # Note the hit name in %blast_taxa.
                

                $scores{$current_score} = $current_hit_ID; # If we haven't already seen this score, make a note.

        } # End of BLASTHIT.
        

        # %scores is a list of scores paired with a single taxonomic ID. Each ID represents either a real BLAST hit or an 'average' for hits with the same score. To recreate the list of BLAST hits, sort by score (descending).
        # This is a roundabout way of doing it, but I intend to implement change 2, so it's temporary.
        my @finished_hit_IDs = ();
        foreach my $finished_hit_ID (sort {$b <=> $a} keys %scores) {
            push (@finished_hit_IDs, $scores{$finished_hit_ID});
        }
        
        
        # Look up the more information about the finished hits using their IDs
        #---------------------------------------------------------------------
        my @all_hit_info;
        foreach my $finished_hit_ID (@finished_hit_IDs) {
            my $finished_hit_name = retrieve_name ($finished_hit_ID, $namesfileDBMids);
            my $finished_hit_rank = retrieve_rank ($finished_hit_ID, $nodesfileDBM);
            my $finished_hit_info = $finished_hit_ID . "\t" . $finished_hit_name . "\t" . $finished_hit_rank;
            push (@all_hit_info, $finished_hit_info); # Add the [ID\tname\trank] for this finished BLAST hit to the end of @all_hit_info.
        }
        
        my $number_of_finished_blast_hits = @all_hit_info;
        if ($number_of_finished_blast_hits == 0) {
            print "\t\tNo hits identified for this read. Something might have gone wrong.\n";
            print $log_filehandle "\t\tNo hits identified for this read. Something might have gone wrong.\n";
            next BLASTENTRY; # Move on to the next read (for which there is a BLAST entry)
        }
        #print "\tNumber of BLAST taxa identified for this read: $number_of_blast_taxa\n";
        #print "\tNumber of BLAST hits for this read after filtering by taxon and score: $number_of_finished_blast_hits\n";
        #print Dumper \@all_hit_info; print "\n\n******\n\n";
		
		
        # Calculate the taxonomic diversity score
        #----------------------------------------
        # Note that this is based on the number of taxa in the BLAST hits, not the number of hits after filtering by score.
		my $tax_diversity_score = ($number_of_blast_taxa/$cap) - (1/$cap); # Remember, $cap defaults to 100.
        
        my $contrastinghit_ID = 0; my $contrastinghit_name = "none found"; # contrastinghit is second BLAST hit: number 1 if you're counting from 0. Default the values to null.
		if ($number_of_finished_blast_hits > 1) { # If there is more than one finished hit:
                my @contrastinghit_info = split ("\t", $all_hit_info[1]);
                $contrastinghit_ID = $contrastinghit_info[0]; # Fetch the contrasting hit ID and name from its info array.
                $contrastinghit_name = $contrastinghit_info[1]; 
		}
		
        # Find the intersection between the top and contrasting hits
        #-----------------------------------------------------------
        # This is the taxon the read will be assigned to.
        my $intersect_ID = 0; my $intersect_rank = "none"; my $intersect_name = "none found"; my $tophit_ID = 0; my $tophit_route = 0; # Default the intersect values to null.
        my @tophit_info = split ("\t", $all_hit_info[0]); # tophit is the top BLAST hit. Unless the top hit wasn't identified, it's also the first BLAST taxon. We ignore unidentified hits.

        $tophit_ID = $tophit_info[0];
        
        if ($number_of_finished_blast_hits > 1) { # If there is more than one BLAST hit:
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

        
		# Find the intersection between the top and bottom hits
		#------------------------------------------------------
		# We also find the intersection between the top hit and the bottom (within $cap): bottom_intersect. This isn't used in any calculations, but it is printed in the intersects file and might be useful one day.
		my $bottom_intersect_ID = 0; my $bottom_intersect_rank = "none found"; my $bottom_intersect_name = "none found"; my $bottomhit_ID = 0; my $bottomhit_name = "none found"; # Default to null values.
		if ($number_of_finished_blast_hits == 2) { # If there were only two hits in the end, the top intersect is the same as the intersect and the bottom hit is the same as the contrasting hit.
            $bottom_intersect_ID = $intersect_ID; $bottom_intersect_rank = $intersect_rank; $bottom_intersect_name = $intersect_name;
            $bottomhit_ID = $contrastinghit_ID; $bottomhit_name = $contrastinghit_name; 
		}
		if ($number_of_finished_blast_hits > 2) { # If there were more than two hits in the end, it gets more complicated.
				my $bottomhit = pop @all_hit_info; # $bottom is the final BLAST hit after filtering. The least good BLAST match (within our standards).
                my @bottomhit_info = split ("\t", $bottomhit);
                $bottomhit_ID = $bottomhit_info[0];
				$bottomhit_name = $bottomhit_info[1];

				unless ($bottomhit_ID == 0) { # If the bottom hit doesn't have an ID, we can't calculate a top intersection.
					my $bottomhit_route = retrieve_taxonomic_structure ($bottomhit_ID, $nodesfileDBM, 1);
					$bottom_intersect_ID = find_taxonomic_intersect ($tophit_route, $bottomhit_route);
				}
                
				if ($bottom_intersect_ID == 0) {			
                    $bottom_intersect_name = "none found";
				} else {
                    $bottom_intersect_name = retrieve_name ($bottom_intersect_ID, $namesfileDBMids);
				}
		}
        
        
		# Print all of this information to the intersects.txt file
		#---------------------------------------------------------
		open (my $intersects_filehandle, ">>".$corename."/"."$corename".".intersects.txt") or die "Cannot write intersects file ".$corename.".intersects.txt\n$!\n"; # Open intersect file for appending.

		print $intersects_filehandle "Query: $query, first hit: $tophit_name ($tophit_ID), expect: $E_value, identities: $Identities, next hit: $contrastinghit_name ($contrastinghit_ID), last hit up to cap: $bottomhit_name ($bottomhit_ID), phylogenetic range of hits up to cap: $bottom_intersect_name ($bottom_intersect_ID), number of hits: $number_of_blast_hits, taxa diversity: $number_of_blast_taxa, taxa diversity score: $tax_diversity_score, classification intersect: $intersect_name ($intersect_ID)\n";
        
        
        close $intersects_filehandle;


        } else { print "No header matches this entry.\n"; } # The end of BLASTENTRY
        
    if ($number_of_headers_to_find == 0) { last; } # We've found the last header. No point searching the read of the BLAST file.
    }
    close $blast_filehandle;
    $/ = "\n"; # Set the read separator back to the default (newline).
    
	return $corename; # Return the $corename path for use in other subroutines.
}


sub retrieve_name {
##### Use taxonomic ID to get name
	my ($query_ID, $namesfileDBMids) = @_; # $query_ID is the taxonomic ID in question. We want to find its name.
	my $name = undef; # Default the name to undef.
    
    unless ($query_ID == 0) { # ID 0 is null. It's not in the names file.
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
##### Find intersection of two taxa by comparing their routes down the phlyogenetic tree
	my ($first_route, $second_route) = @_; # Routes are the taxonomic IDs of successive parent taxa.
	my $intersect_ID = 0; # If either the first or second route aren't defined, return the intersect ID 0.
	if (defined $first_route && defined $second_route){
		my @first_route = split (/\t/, $first_route);
		my @second_route = split (/\t/, $second_route);
		my $connected = 0; # The $connected flag activates when the first shared ID is stored as $intersect_ID and prevents it from being overwritten. Because the routes go from lower to higher taxa, higher shared IDs are ignored.
		foreach my $first_ID (@first_route) { # Start with each ID in the first route.
			foreach my $second_ID (@second_route) { # Take the IDs in the second route.
				if ($first_ID == $second_ID) { # If the IDs match,
					unless ($connected == 1) { # And if they're not yet flagged as connected,
						$connected = 1; # Flag them as connected.
						$intersect_ID = $first_ID; # $intersect_ID becomes this shared rank.
					} # If the routes are already flagged as connected, move on to the next ID in the second route. TO DO: stop processing when a connection is found?
				}
			}
		}
	}
	return $intersect_ID;	
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
            $rank = $node_info[1]; # Rank is the 1st field.
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
    my $route = undef; # Default the route to undef.
    
    unless ($query_ID == 0) {
    
        my $exit = 0; # When to exit the do loop below.
        my $next_level_ID; # The ID of the parent node.
        my $rank; # The rank of the current ID.
        my @route = (); # @route is a list of tab-separated taxonomic IDs moving down from the current node.
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
    
        while (@current_hit_name) { # Do the following until last is called or you run out of words in @current_hit_name.
                        my $current_hit_name = join (' ', @current_hit_name);
                        #print "\tAttempting '$current_hit_name'\n";
                             
                        if (exists $namesfileDBMnames{$current_hit_name}) { # If the BLAST name has a match in the names file,
                            my @namesfile_matches = split ("\t", $namesfileDBMnames{$current_hit_name}); # See if there is more than one ID under this name.
                            if (exists $namesfile_matches[1]) { # If there's more than one ID, things gonna get complicated.
                                    my @potential_taxa_in_right_range = (); # If there's more than one ID in the expected phylogenetic range, things gonna get even more complicated.
    
                                    foreach my $potential_taxon (@namesfile_matches) { # For each potential taxon, pull out its clarification and see whether it is within the expected phylogenetic range. For example, if an ID is in Insecta but the reads should all be in Viridiplantae, the ID is not the right one.       
                                            if ($potential_taxon == $expected_phylogenetic_range) { # If it's the same as the range, accept it.
                                                    push (@potential_taxa_in_right_range, $potential_taxon);
                                            } else {
                                                    my @potential_parent_taxa = split("\t", retrieve_taxonomic_structure ($potential_taxon, $nodesfileDBM, 0) ); # If not, check the parent IDs.
                                                    @potential_parent_taxa = reverse @potential_parent_taxa; # The expected phylogenetic range is probably a high taxon, so start with the higher parent taxa.
                                                    foreach my $potential_parent_taxon (@potential_parent_taxa) {
                                                        if ($potential_parent_taxon == $expected_phylogenetic_range) {
                                                            push (@potential_taxa_in_right_range, $potential_taxon);
                                                            
                                                        }
                                                    }
                                            }
                                    }
                                    
                                    my $number_of_potential_taxa_in_right_range = @potential_taxa_in_right_range;
                                    if ($number_of_potential_taxa_in_right_range == 1) { # TO DO: test
                                        $current_hit_ID = $potential_taxa_in_right_range[0];
                                        return ($current_hit_ID, $not_in_expected_range); # We've assigned a taxon to the BLAST name, so exit the search.
                                    }
                                    
                                    # If multiple potential taxa were in the expected phylogenetic range, choose the highest. We can't distinguish between them any further, so assign cautiously.
                                    # If multiple taxa have the joint highest rank, it will pick randomly because they are stored in a hash.
                                    my %potential_taxa_in_right_range_ranks = ();
                                    print "\t\tMultiple matches for '$original_hit_name':";
                                    print $log_filehandle "\t\tMultiple matches for '$original_hit_name':"; 
                                    foreach my $potential_taxon (@potential_taxa_in_right_range) { # Find the taxonomic rank for each one. Print while we're at it.
                                        print "\t$potential_taxon";
                                        print $log_filehandle "\t$potential_taxon";
                                        my $rank = retrieve_rank($potential_taxon, $nodesfileDBM);
                                        if (exists $potential_taxa_in_right_range_ranks{$rank}) {
                                            $potential_taxa_in_right_range_ranks{$rank} = $potential_taxa_in_right_range_ranks{$rank} . "\t" . $potential_taxon;
                                        } else {
                                            $potential_taxa_in_right_range_ranks{$rank} = $potential_taxon;
                                        }
                                    }
                                    
                                    my @ranks = ('superkingdom', 'kingdom', 'phylum', 'subphylum', 'class', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'family', 'genus', 'species', 'subspecies', 'no rank'); # The ranks are ordered highest to lowest. It checks for the highest first.
                                    foreach my $rank (@ranks) {
                                        if (exists $potential_taxa_in_right_range_ranks{$rank}) { # If at least one taxon matches this rank:
                                            my @potential_taxa_right_rank = split ("\t", $potential_taxa_in_right_range_ranks{$rank}); # Save them in an array.
                                            if ($potential_taxa_right_rank[1]) { # If there are multiple taxa, choose between them randomly.
                                                $current_hit_ID = $potential_taxa_right_rank[rand @potential_taxa_right_rank];
                                                print "\n\t\t\tMultiple matches at highest rank ($rank). Randomly assigned to $current_hit_ID\n";
                                                print $log_filehandle "\n\t\t\tMultiple matches at highest rank ($rank). Randomly assigned to $current_hit_ID\n";
                                                return ($current_hit_ID, $not_in_expected_range); # We've assigned a taxon to the BLAST name, so exit the search.
                                            } else {
                                                $current_hit_ID = $potential_taxa_in_right_range_ranks{$rank}; # If there's only one taxon at this rank, assign to it.
                                                print "\n\t\t\tAssigned to the highest-ranked match: $current_hit_ID ($rank)\n";
                                                print $log_filehandle "\n\t\t\tAssigned to the highest-ranked match: $current_hit_ID ($rank)\n";
                                                return ($current_hit_ID, $not_in_expected_range); # We've assigned a taxon to the BLAST name, so exit the search.
                                            }
                                        }
                                    }

                                    $not_in_expected_range = 1; # If there were no potential taxa in the expected phylogenetic range, flag and exit subroutine.
                                    return ($current_hit_ID, $not_in_expected_range);
    
                            } else { # If there's just one ID for this name, things are a lot simpler.
                                    $current_hit_ID = $namesfile_matches[0];
                                    return ($current_hit_ID, $not_in_expected_range); # Exit the subroutine.
                            }
                        } else {
                            pop @current_hit_name; # If no match, take a word off the end of $current_hit_name and try again.
                        }
        } # End of while loop. If a name makes it this far, it wasn't found in the names file and its ID will remain at 0.
        return ($current_hit_ID, $not_in_expected_range);
}