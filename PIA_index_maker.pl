#!/usr/bin/perl 

######################################################
################# PIA_index_maker.pl #################
######################################################
####### Creates DBM index files for later use ########
############## Becky Cribdon, UoW 2019################
######################################################
############## Version 4.8, 2019-09-04 ###############
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
# - Makes an index of the names file where keys are IDs and values are the matching scientific names.
# - Makes an index of the nodes file where keys are IDs and values are parent nodes and taxonomic ranks.

my %options=();
getopts('n:N:"', \%options);
    
##### Locate reference files #####
	# If they aren't in their default location (Reference_files/), use the options to specify where they are.
	my $nodes_filename = 'Reference_files/nodes.dmp';
	if ($options{n}){
		$nodes_filename = ($options{n});
	}
	
    my $names_filename = 'Reference_files/names.dmp';   
	if ($options{N}){
		$names_filename = ($options{N});
	}


##############################################		
########### Set up DBM index files ###########
##############################################

print "\n\n****Setting up DBM index files****\n";


# Nodes index where keys are IDs
#-------------------------------
my $nodesfileDBM = $nodes_filename . '.dbm';

if (-e $nodesfileDBM) {
    print "$nodesfileDBM already exists. Overwriting.\n";
    unlink $nodesfileDBM;
}
   
my %nodesfileDBM = (); # Keys are taxonomic IDs and values are parent nodes and taxonomic ranks.
tie (%nodesfileDBM, "DB_File", $nodesfileDBM, O_RDWR|O_CREAT, 0666, $DB_BTREE) or die "Can't open $nodesfileDBM.\n$!\n";

# Now populate the hash:
open (my $nodes_filehandle, $nodes_filename) or die "Could not open $nodes_filename.\n$!\n";
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

    
# Names index where keys are IDs
#-------------------------------
# This will be used to look up the scientific names of taxa using their IDs.
my $namesfileDBM = $names_filename . '.dbm';

if (-e $namesfileDBM) {
    print "$namesfileDBM already exists. Overwriting.\n";
    unlink $namesfileDBM;
}
   
my %namesfileDBM = (); # Keys are taxonomic IDs and values are scientific names.
tie (%namesfileDBM, "DB_File", $namesfileDBM, O_RDWR|O_CREAT, 0666, $DB_BTREE) or die "Can't open $namesfileDBM.\n$!\n";

# Now populate the hash:
open (my $names_filehandle, $names_filename) or die "Could not open $names_filename.\n$!\n";
while (1) { # Run this loop until "last" is called.
        my $line = <$names_filehandle>; # Read the next line from the names file.
        if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.

        my @line = split(/\|/, $line); # Split the line by | characters.

        if ($line[3] eq "\tscientific name\t") {
            my $name = $line[1];
            $name =~ s/^\s+|\s+$//g; # Trim whitespace off the start and end.
            my $ID = $line[0];
            $ID =~ s/^\s+|\s+$//g;
            $namesfileDBM{$ID} = $name; # Store in the hash.
        } # If that line is not for a scientific name, ignore it and move on.
}
close $names_filehandle;
untie %namesfileDBM;


print "Finished indexing.\n\n";  
