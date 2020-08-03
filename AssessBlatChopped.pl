#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f );

# Usage
my $usage = "
AssessBlatChopped.pl - outputs matches from a blat file where choppsed sequences are compared

Copyright (C) 2020 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl AssessBlatChopped.pl options
 required:
  -f  (path to) psl file
";

#############

getopts('f:');
die $usage unless ($opt_f);

my $blat = $opt_f if $opt_f;

my @out;

my $count = 0;

open(IN, $blat) || die "can't open $blat\n";

while (<IN>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ((defined $data[16])&&($data[16] =~ /\d/)) {
    $count +=1;
    my @name1 = split /\./, $data[9];
    my @name2 = split /\./, $data[13];
    push @out, "$count\t$data[0]\t$name1[-1]\t$name2[-1]";
  }
}

close (IN);

my $result = join "\n", @out;

my $outfile = "$blat.matches.txt";

unless ( open(META, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print META "Count\tMatch\tPos1\tPos2\n$result";
close (META);
