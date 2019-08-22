#!/usr/bin/perl 

######################################################
################# PIA_index_maker.pl #################
######################################################
####### Creates DBM index files for later use ########
############## Becky Cribdon, UoW 2019################
######################################################
############## Version 1.0, 2019-08-22 ###############
######################################################

	use strict;
	use warnings;
    use Getopt::Std;
    use DB_File;
    use Fcntl;
	use Data::Dumper qw(Dumper);
    
# Run as follows:
# perl PIA_index_maker.pl [optional reference file locations]
# It will take a few minutes.

# Summary:
# - Makes an index of the names file where keys are names and values are any matching taxonomic IDs (tab-delimited).
# - Makes an index of the names file where keys are IDs and values are the matching scientific names.
# - Makes an index of the nodes file where keys are IDs and values are parent nodes and taxonomic ranks.

my %options=();
getopts('n:N:', \%options);
    
##### Locate nodes.dmp and names.dmp files #####
	# If nodes.dmp and names.dmp aren't in their default location (Reference_files/), use -n and -N to specify where they are.	
	my $nodesfile="Reference_files/nodes.dmp";
	if ($options{n}){
		$nodesfile=($options{n});
	}
	my $namesfile="Reference_files/names.dmp";   
	if ($options{N}){
		$namesfile=($options{N});
	}


##############################################		
########### Set up DBM index files ###########
##############################################
print "\n\n****Setting up DBM index files****\n";


# Nodes index where keys are IDs
#-------------------------------
my $nodesfileDBM = $nodesfile . '.dbm';

if (-e $nodesfileDBM) {
    print "$nodesfileDBM already exists. Overwriting.\n";
    unlink $nodesfileDBM;
}
   
my %nodesfileDBM = (); # Keys are taxonomic IDs and values are parent nodes and taxonomic ranks.
tie (%nodesfileDBM, "DB_File", $nodesfileDBM, O_RDWR|O_CREAT, 0666, $DB_BTREE) or die "Can't open $nodesfileDBM.\n$!\n";

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
close $nodes_filehandle;
untie %nodesfileDBM;


# Names index where keys are names
#---------------------------------
# This will be used to look up the taxonomic IDs of BLAST hits.
my $namesfileDBMnames = $namesfile . '_names.dbm'; # Including $header_filename in the DBM filename ensures that every thread makes its own index file, preventing read/write conflict. Otherwise I'd have to get into file locks and oh my goodness not now.
        
if (-e $namesfileDBMnames) {
    print "$namesfileDBMnames already exists. Overwriting.\n";
    unlink $namesfileDBMnames;
}
        
my %namesfileDBMnames = (); # Keys are names (scientific and not, because BLAST contains all sorts) and values are the corresponding taxonomic ID. It is built as a hash and then written to the DBM file.
tie (%namesfileDBMnames, "DB_File", $namesfileDBMnames, O_RDWR|O_CREAT, 0666, $DB_BTREE) or die "Can't open $namesfileDBMnames.\n$!\n"; # Open the DBM file (in read-write mode; O_RDWR). If it doesn't exist, create one (O_CREAT). Note that these options use the Fcntl module.

# Now populate the hash:
open (my $names_filehandle, $namesfile) or die "Could not open $namesfile.\n$!\n";
while (1) { # Run this loop until "last" is called.
        my $line = <$names_filehandle>; # Read the next line from the names file.
        if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.

        my @line = split(/\|/, $line); # Split the line by | characters.
        my $name_full = $line[1]; # Pick out the name.
        $name_full =~ s/^\s+|\s+$//g; # Trim whitespace (tabs) off the start and end.
        my $name; # We may have to shorten it.

        if (length ($name_full) > 52) {
            $name = substr ($name_full, 0, 52); # We only want the first 52 characters of the name (counting from 0) because that's the maximum length for most BLAST hits. Hits can be up to 55 characters intact, but 56 or higher get truncated so an ellipsis (...) can go on the end.
            # A few taxa have very long names that only differ at the very end; so far, it seems these don't match and are discarded.
            my @name = split (' ', $name); # Now split on spaces. This also removes any trailing space left by the character trim.
            if ($name[1]) {
                pop @name; # Remove the last word, which since $name_full was more than 52 characters, will have been truncated.
            }
            $name = join (' ', @name);
        } else {
            $name = $name_full;
        }

        my $ID = $line[0]; # Also get the ID.
        $ID =~ s/^\s+|\s+$//g;

        if (exists $line[2]) { # Examine the second names field. This is usually empty, but otherwise gives clarification if the name is also used for another taxon. Or other extra information that isn't relevant to us. But a second names field suggests that this $name may be a duplicate, so append this ID rather than potentially overwriting a previous one.
            if (exists $namesfileDBMnames{$name}) {
                $namesfileDBMnames{$name} = $namesfileDBMnames{$name} . "\t$ID";
            } else {
                $namesfileDBMnames{$name} = "$ID";
            }
        } # Only doing a hash lookup if there's a second names field saves us doing a hash lookup every time we add to the hash.
}
close $names_filehandle;
untie %namesfileDBMnames;

    
# Names index where keys are IDs
#-------------------------------
# This will be used to look up the scientific names of taxa using their IDs.
my $namesfileDBMids = $namesfile . '_IDs.dbm';

if (-e $namesfileDBMids) {
    print "$namesfileDBMids already exists. Overwriting.\n";
    unlink $namesfileDBMids;
}
   
my %namesfileDBMids = (); # Keys are taxonomic IDs and values are scientific names.
tie (%namesfileDBMids, "DB_File", $namesfileDBMids, O_RDWR|O_CREAT, 0666, $DB_BTREE) or die "Can't open $namesfileDBMids.\n$!\n";

# Now populate the hash:
open ($names_filehandle, $namesfile) or die "Could not open $namesfile.\n$!\n";
while (1) { # Run this loop until "last" is called.
        my $line = <$names_filehandle>; # Read the next line from the names file.
        if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.

        my @line = split(/\|/, $line); # Split the line by | characters.

        if ($line[3] eq "\tscientific name\t") {
            my $name = $line[1];
            $name =~ s/^\s+|\s+$//g; # Trim whitespace off the start and end.
            my $ID = $line[0];
            $ID =~ s/^\s+|\s+$//g;
            $namesfileDBMids{$ID} = $name; # Store in the hash.
        } # If that line is not for a scientific name, ignore it and move on.
}
close $names_filehandle;
untie %namesfileDBMids;


print "Finished indexing.\n\n";  
