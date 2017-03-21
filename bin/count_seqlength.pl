#!/usr/bin/env perl

use strict;
use warnings;

my $file=shift;

open FH, "<$file";

while (<FH>) {
	if (/^>/) {
		print;
	}
	else {
		my $seq=$_;
		chomp $seq;
		my $len=(length($seq));
		$seq=~ s/N//g;
		my $newlen=length($seq);
		my $div = $newlen/3;
		print "full length $len\t length without Ns $newlen AA len $div\n";
	}
}
