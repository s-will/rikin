package RRIkin;

use 5.008003;
use strict;
use warnings;


# set the version for version checking
our $VERSION     = 1.00;

our @ISA         = qw(Exporter);
our @EXPORT      = qw(
  rev_compl
  read_table
  lookup_table
);


## ------------------------------------------------------------
## subs

## reverse complement
sub rev_compl {
    my $s=shift;
    $s =~ tr/ACGT/TGCA/;
    return reverse($s);
}

############################################################
## @brief read table from file with header
## @param filename name of the input file
## @return list of hash; one hash per row; names are indices
############################################################
sub read_table {
    my $tab_filename=shift;

    my @tab_header;
    my @tab=();
    
    open(IN,$tab_filename) || die "Cannot read input file $tab_filename.\n";
    while(my $line=<IN>) {
	next if ($line =~ /^\s*#/); ## skip comments
	
	if (@tab_header==0) {
	    @tab_header=split /\s+/, $line;
	    next;
	}
	my @list=split /\s+/, $line;
	my %linehash;
	for my $i (0..@tab_header-1) {
	    $linehash{$tab_header[$i]} = $list[$i];
	}
	push @tab, { %linehash };
    }
    close IN;
    
    return \@tab;
}

############################################################
## @brief lookup entry in table
## @param tab
## @param key
## @param value
## @return row hash
############################################################
sub lookup_table {
    my $tab=shift;
    my $key=shift;
    my $value=shift;
    
    for my $row (@$tab) {
	if ($row->{$key} eq $value) {
	    return $row;
	}
    }
    return;
}

## ------------------------------------------------------------


1;
