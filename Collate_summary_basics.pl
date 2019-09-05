#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper qw(Dumper); # Only used for commented-out test prints.

# Version 1.1
# Updated 2019-09-04

# This script collates Summary_Basic.txt files. It can optionally include a pre-PIA spreadsheet from MEGAN. The MEGAN spreadsheet must be in "taxonID_to_count" format.
# It uses fullnamelineage.dmp from the NCBI: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
#
# Run as follows:
# >perl Collate_spreadsheets.pl -f [optional path to fullnamelineage.dmp] -m [optional MEGAN spreadsheet] -o [optional output name] [at least one summary basic]

# Summary
#--------
# If there is a MEGAN file (slightly different format to summary basics),
#   Edit the header for outputting later.
#   For each other line,
#       If no hits, skip.
#       Otherwise, save in the master hash. Keys are ID, values are rest of line.
#   Find full taxonomies of all IDs and add to hash values.
#
# For each summary basic,
#   Add header section to master header.
#   For each other line,
#       Save in the incoming hash. Keys are ID, values are hit number.
#   For each element in the master hash,
#       If also in the incoming hash, add new hit number to the end of the value. Delete from incoming hash.
#       If not also in incoming hash, add a new hit number of 0.
#   For elements remaining in the incoming hash,
#       Find full taxonomies and add to master hash. Delete from incoming hash.
#
# Print master hash to output.


#######################################################################################################################################################################

# Setting up
#===========
print "\n\n##### Welcome to the PIA output collator #####\n\n";

my %options = ();
getopts('f:m:o:h', \%options);

if ($options{h}) { # If the help option is called, print the help text and exit.
    print "Usage: perl Collate_summary_basics.pl [-fmoh] [summary basics]

Option	Description			Explanation
-f	Full taxonomy file name		Path to fullnamelineage.dmp if it is not in the current directory. 
-m	MEGAN file name			Pre-PIA MEGAN output to collate the summary basics with. Must be in the format taxonID_to_count.
-o	Output file name                Name of output file. Defaults to 'Collapsed.txt'.
-h	Help				Print this message.

Other arguments	    Description
[summary basics]    Output files from the PIA in the format ID, name, hit count.

";
	exit;
}

my $full_taxonomy_filename = 'fullnamelineage.dmp'; # Default path to full taxonomy file (i.e. in this directory).
if ($options{f}) { # If the full taxonomy option is called, overwrite the default with the option input.
    $full_taxonomy_filename = $options{f};
    print "-f: Using full taxonomy file '$full_taxonomy_filename'\n";
}

my $megan_filename = ($options{m}); # If the master file option is called, take the input as the master file name.
if ($megan_filename) {
    print "-m: Using pre-PIA MEGAN file '$megan_filename'\n";
}

my $output_filename = 'Collated.txt'; # Default output file name.
if ($options{o}) { # If the output name option is called, overwrite the default with the option input.
    $output_filename = $options{o};
    print "-o: Using output file name '$output_filename'\n";
}


my @summary_basic_filenames = @ARGV; # Any other inputs are summary basics.

print "\nWorking on:\n";

my $header_master;
my $number_of_sample_columns = 0; # Start at 0.
my %data_by_ID = (); # Keys are taxonomic IDs and values are everything else to be outputted for that taxon.


#######################################################################################################################################################################

# Process the MEGAN file (optional)
#==================================
if ($megan_filename) { # If there is a MEGAN file to be processed.

    print "\t$megan_filename\n";
    open (my $megan_filehandle, $megan_filename) or die "Could not open file $megan_filename.\n$!\n";
    
    my %data_by_ID_partial = (); # Keys are taxonomic IDs and values will eventually be everything else to be outputted on that line APART FROM full taxonomies.
    
    $header_master = readline($megan_filehandle); # Save the header. We'll add to this later.
    $header_master =~ s/\s+$//; # Remove trailing whitespace of any kind. Files edited in Excel can have carriage returns as well as newlines.
    my @header_master = split ("\t", $header_master);
    $number_of_sample_columns = scalar @header_master - 1; # Note the number of sample columns.
    $header_master[0] = "Full_taxonomy\tName\tTaxonomic_ID"; # Change the 0th element ("#Datasets") to the final non-sample column titles.
    $header_master = join ("\t", @header_master); # Join the line back into a scalar.
    
    # For every remaining line,
    MASTERLINES: foreach my $line (readline($megan_filehandle)) {
        my $line_sum = 0; # Reset $line_sum. Otherwise it's cumulative across lines. 
        my @line_elements = split("\t", $line); # Split the line on tabs.
        
        foreach my $element (@line_elements[1 .. $#line_elements]) { # Sum all elements but the 0th. That's the taxonomic ID.
            $line_sum += $element;
        }
        if ($line_sum == 0) { next MASTERLINES;}  # If there aren't any hits to this taxon, skip it.
        
        # Otherwise, save the line in %data_by_ID.
        my $ID = $line_elements[0];
        my $rest_of_line = join ("\t", @line_elements[1 .. $#line_elements]);
        $rest_of_line =~ s/\s+$//; # Remove trailing whitespace of any kind.
        $data_by_ID_partial{$ID} = $rest_of_line;
    }
    close $megan_filehandle;
    
    if (exists $data_by_ID_partial{131567}) { # If one of the taxa is ID 131567, it's "cellular organisms" and doesn't have a full taxonomy. Fill its information for $data_by_ID here to avoid formatting problems.
        $data_by_ID{131567} = "cellular organisms\tcellular organisms\t" . $data_by_ID_partial{131567};
        delete $data_by_ID_partial{131567};
    }
    
    # Look up the IDs in fullnamelineage.dmp
    open (my $full_taxonomy_filehandle, $full_taxonomy_filename) or die "Could not open $full_taxonomy_filename.\n$!\n";
    foreach my $line (readline($full_taxonomy_filehandle)) {
        my @line = split(/\|/, $line); # Split the line by | characters.
        my $ID = $line[0]; # Element [0] is the ID.
        $ID =~ s/\s+//g; # Remove any whitespace.
        
        if (exists $data_by_ID_partial{$ID}) { # If this is an ID we're looking for:
            
            $line[1] =~ s/^\s+|\s+$//g; # $line[1] is the scientific name. Remove any leading or trailing whitespace.
            $line[2] =~ s/^\s+|\s+$//g; # $line[2] is the taxonomy. Remove any leading or trailing whitespace.
            my $full_taxonomy = $line[2] . ' ' . $line[1]; # Add scientific name onto the end of the taxonomy to get the full taxonomy. Apart from this, the taxonomy is still in NCBI format.
            
            $data_by_ID{$ID} = "$full_taxonomy\t$line[1]\t" . $data_by_ID_partial{$ID}; # Save the new information in %data_by_ID.
            delete $data_by_ID_partial{$ID};
        }
    }
    
    foreach my $remaining_ID (keys %data_by_ID_partial) {
        print "\t\tWarning: failed to find $remaining_ID in fullnamelineage.dmp. Defaulting name and taxonomy to 'none found'.\n";
        $data_by_ID{$remaining_ID} =  "none found\tnone found\t" . $data_by_ID_partial{$remaining_ID}; # The full taxonomy and name are 'none found'.
    }
    
    close $full_taxonomy_filehandle;

    } else {
    $header_master = "Full_taxonomy\tName\tTaxonomic_ID"; # If there wasn't a MEGAN file to process, set up these variables ready for the first summary basic.
}


#######################################################################################################################################################################

# Process the summary basics
#===========================
foreach my $summary_basic_filename (@summary_basic_filenames) {
    print "\t$summary_basic_filename\n";
    open (my $summary_basic_filehandle, $summary_basic_filename) or die "Could not open $summary_basic_filename.\n$!\n";

    # Save every other line in %incoming_by_ID. This hash resets with each new sample.
    my %incoming_by_ID = ();
    SAMPLELINE: foreach my $line (readline($summary_basic_filehandle)) {
        #print "|$line|";
        if (index ($line, '#') != -1) { # If the line contains a hash symbol, which indicates the header line,
            if (index ($line, 'Input FASTA:') != -1) { # If it's the FASTA line, extract the FASTA name.
                my @fasta_field = split ('Input FASTA:', $line);
                my $fasta_name = $fasta_field[1];
                $fasta_name =~ s/\s+$//; # Remove trailing whitespace of any kind, just in case.
                $header_master = $header_master . "$fasta_name"; # Add this sample to the master header. Note that it already has a preceding tab.
            }
            next; # Then move on.
        }
        
        my @line = split ("\t", $line);
        $line[2] =~ s/\s+$//; # Remove trailing whitespace of any kind, just in case.
        $incoming_by_ID{$line[0]} = $line[2]; # Key is ID; value is hit count.
    }

    # Combine common taxa in %data_by_ID and %incoming_by_ID.
    foreach my $existing_taxon (keys %data_by_ID) {
        if (exists $incoming_by_ID{$existing_taxon}) { # If this taxon is also in the incoming summary basic,
            $data_by_ID{$existing_taxon} = $data_by_ID{$existing_taxon} . "\t" . $incoming_by_ID{$existing_taxon}; # Add the incoming hit count on the end of the value.
            delete $incoming_by_ID{$existing_taxon};
        } else {
            $data_by_ID{$existing_taxon} = $data_by_ID{$existing_taxon} . "\t0"; # Otherwise, add a 0.
        }
    }

    
    # Add any new taxa to %data_by_ID after looking up their full taxonomies.
    my $preceding_empty_columns = "\t0" x $number_of_sample_columns;
    
    if (exists $incoming_by_ID{131567}) { # If one of the taxa is ID 131567, it's "cellular organisms" and doesn't have a full taxonomy. Fill its information for $data_by_ID here to avoid formatting problems.
        $data_by_ID{131567} = "cellular organisms\tcellular organisms$preceding_empty_columns\t" . $incoming_by_ID{131567};
        delete $incoming_by_ID{131567};
    }
        
    open (my $full_taxonomy_filehandle, $full_taxonomy_filename) or die "Could not open $full_taxonomy_filename.\n$!\n";
    foreach my $line (readline($full_taxonomy_filehandle)) {
        my @line = split(/\|/, $line);   # Split the line by | characters.
        my $ID = $line[0];     # Element [1] is the ID.
        $ID =~ s/\s+//g; # Remove any whitespace.
        
        if (exists $incoming_by_ID{$ID}) { # If this is an ID we're looking for:
    
            $line[1] =~ s/^\s+|\s+$//g; # $line[1] is the scientific name. Remove any leading or trailing whitespace.
            $line[2] =~ s/^\s+|\s+$//g; # $line[2] is the taxonomy. Remove any leading or trailing whitespace.
                 
            my $full_taxonomy = $line[2] . ' ' . $line[1]; # Add scientific name onto the end of the taxonomy to get the full taxonomy. Apart from this, the taxonomy is still in NCBI format.
            
            $data_by_ID{$ID} = "$full_taxonomy\t$line[1]$preceding_empty_columns\t" . $incoming_by_ID{$ID}; # Save the information from here and %incoming_by_ID in %data_by_ID.
            delete $incoming_by_ID{$ID};
        }
    }
    close $full_taxonomy_filehandle;
    
    foreach my $remaining_ID (keys %incoming_by_ID) {
        print "\t\tWarning: failed to find $remaining_ID in fullnamelineage.dmp. Defaulting name and taxonomy to 'none found'.\n";
        $data_by_ID{$remaining_ID} = "none found\tnone found$preceding_empty_columns\t" . $incoming_by_ID{$remaining_ID}; # The full taxonomy and name are 'none found'.
    }


    $number_of_sample_columns ++;
    close $summary_basic_filehandle;
}


#######################################################################################################################################################################

# Print output
#=============
open( my $output_filehandle, '>', $output_filename) or die "Cannot open $output_filename for writing.\n$!\n";

print $output_filehandle "$header_master\n"; # Print the header first.

foreach my $ID (keys %data_by_ID) {
    my @entry = split ("\t", $data_by_ID{$ID});
    my $sample_columns = join ("\t", @entry[2..$#entry]); # Extract the sample columns - that's column 3 onwards.
    print $output_filehandle "$entry[0]\t$entry[1]\t$ID\t$sample_columns\n"; # Print the full taxonomy, name, ID, then sample columns.
}

if ($output_filename eq 'Collated.txt') {
    print "\nExported collated spreadsheet as 'Collated.txt'\n";
}

print "\nDone.\n\n";