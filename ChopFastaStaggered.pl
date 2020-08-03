#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_w);

# Usage
my $usage = "
ChopFastaStaggered.pl
Breaks up a single-sequence fasta into 600bp chunks, in staggered windows

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

Usage: perl ChopFastaStaggered.pl options
 required:
  -f	(path to) fasta
 optional:
  -w  width of each chopped segment [default = 600]
";

#############

getopts('f:w:');
die $usage unless ($opt_f);

my $fasta = $opt_f if $opt_f;

my $subsetlength;

if (defined $opt_w) {
  $subsetlength = $opt_w;
} else {
  $subsetlength = 600;
}

open(FAS, $fasta) || die "can't open $fasta\n";

my $name;

my $subseqcount = 1;

my $nuccount1 = 0;

my $nuccount2 = 0;

my @out;

my @window1;

my @window2;

my $pause = 0;

while (<FAS>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  if ($line =~ /^>/) {
    $name = $line;
  } else {
    if ($nuccount1 >= $subsetlength) {
      my $seq = join "\n", @window1;
      push @out, "$name.$subseqcount\n$seq";
      $subseqcount +=1;
      $nuccount1 = 0;
      @window1 = ();
    }
    if ($nuccount2 >= $subsetlength) {
      my $seq = join "\n", @window2;
      push @out, "$name.$subseqcount\n$seq";
      $subseqcount +=1;
      $nuccount2 = 0;
      @window2 = ();
    }
    push @window1, $line;
    $nuccount1 += length ($line);
    if ($pause >=5) {
      push @window2, $line;
      $nuccount2 += length ($line);
    }
    $pause +=1;
  }
}

close (FAS);

my $result = join "\n", @out;

my $outfile = "$fasta.chopped$subsetlength"."staggered.fas";

unless ( open(META, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print META $result;
close (META);
