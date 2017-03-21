#!/usr/bin/perl

use strict;
use warnings;

my $input=shift;
my $output=shift;

open FH, "<$input";
open OUT, ">$output";

while (<FH>) {
    if (/(Library_ID\S+)/) {
        print OUT "$1\n";
    }
    if (/^\S+?,(\S+.*)/) {
        print OUT "$1\n";
    }
}
