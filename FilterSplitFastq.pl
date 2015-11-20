use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_n $opt_q $opt_g $opt_o);

# Usage
my $usage = "
FilterSplitFastqs.pl - reads fastq files that consists of multiple samples identified in read title. Splits samples into individual files and renames them.
Copyright (C) 2015 by Jacob A Tennessen

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

Usage: perl FilterSplitFastqs.pl options
 required:
  -f	a comma-delimited list of FASTQ files.
  -n	a tab-delimited table with two columns. First column is final samples name (will be file name), second column is identified name in the FASTQ read title
 optional:
  -q	ASCII character for miniumum quality score for a base to be retained (otherwise connverted to N) [default = 4 (qual of 20)]
  -g    mimimum number of good sites in a read for the read to be retained [default = 36]
  -o    output folder [default is same folder as first FASTQ file]
";

#############

# command line processing.
getopts('f:n:q:g:o:');
die $usage unless ($opt_f);
die $usage unless ($opt_n);

my ($fastqs, $namefile, $qualthresh, $goodthresh, $outfolder);

$fastqs = $opt_f if $opt_f;
$namefile = $opt_n if $opt_n;
$qualthresh = $opt_q ? $opt_q : 4;
$goodthresh = $opt_g ? $opt_g : 36;
$outfolder = $opt_o if $opt_o;

my @fastqs = split ",", $fastqs;

my %fastqs;

foreach my $f (@fastqs) {
    $fastqs{$f} = 1;
    unless (defined $outfolder) {
	my @filedata = split /\//, $f;
	my $finalfile = pop @filedata;
	if (defined $filedata[0]) {
	    push @filedata, "/";
	    $outfolder = join "/", @filedata;
	}
    }
}

open(EXP, "$namefile") || die "can't open $namefile\n";

my %inds;

while (<EXP>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    $inds{$data[1]} = $data[0];
}

close (EXP);

my $outputsize = 500; #how often to output data

my %HoAgoods;

my %cleared;

my $seenlines = 0;

my $goodlines = 0;

foreach my $fq (keys %fastqs) {
    
    unless (open(FASTQ, "$fq")) {
        print "Cannot open file \"$fq\"\n\n";
        next;
    }
    
    my $linetype = 1;
    
    my @fastqdata;
    
    while (<FASTQ>) {
        my $line = $_;
        $line =~ s/\r|\n//g;
        push @fastqdata, $line;
        if ($linetype == 4) {
            my @linenamedata = split ":", $fastqdata[0];
            my $lastbit = pop @linenamedata;
            my @linenamedata2 = split /\s+/, $lastbit;
            unless (defined $inds{$linenamedata2[0]}) {
                print "Name missing for $linenamedata2[0].\n";
                exit;
            }
            my $finalname = $inds{$linenamedata2[0]};
            $seenlines +=1;
            my @seq = split "", $fastqdata[1];
            my @quals = split "", $fastqdata[3];
            my $qualcount = 0;
            my $reject = 0;
            my $firstgood;
            my $lastgood;
            my $goodcount = 0;
            foreach my $q (@quals) {
                if (($q cmp $qualthresh) == -1) {
                    $seq[$qualcount] = "N";
                } else {
                    $goodcount +=1;
                    unless (defined $firstgood) {
                        $firstgood = $qualcount;
                    }
                    $lastgood = $qualcount;
                }
                $qualcount +=1;   
            }
            if ($goodcount < $goodthresh) {
                $reject = 1;
            }
            unless ($reject == 1) {
                $fastqdata[1] = join "", @seq[$firstgood..$lastgood];
                $fastqdata[3] = join "", @quals[$firstgood..$lastgood];
                my $goodline = join "\n", @fastqdata;
                push @{$HoAgoods{$finalname}}, $goodline;
                $goodlines +=1;
            }
            @fastqdata = ();
            $linetype = 0;
            foreach my $fn (keys %HoAgoods) {
                if (scalar (@{$HoAgoods{$fn}}) >= $outputsize) {
                    my $result = join "\n", (@{$HoAgoods{$fn}});
		    my $out = "$fn.fq";
		    if (defined $outfolder) {
			$out = "$outfolder$fn.fq";
		    }
                    if (defined $cleared{$fn}) {
                        unless ( open(META, ">>$out") ) {
                            print "Cannot open file \"$out\" to write to!!\n\n";
                            exit;
                        }
                        print META "$result\n";
                        close (META);                        
                    } else {
                        unless ( open(META, ">$out") ) {
                            print "Cannot open file \"$out\" to write to!!\n\n";
                            exit;
                        }
                        print META "$result\n";
                        close (META);
                        $cleared{$fn} = 1;                        
                    }
                    (@{$HoAgoods{$fn}}) = ();
                }
            }
        }
        $linetype +=1;
    
    }
    
    close (FASTQ);
   
}

foreach my $fn (keys %HoAgoods) {
    if (scalar (@{$HoAgoods{$fn}}) > 0) {
        my $result = join "\n", (@{$HoAgoods{$fn}});
	my $out = "$fn.fq";
	if (defined $outfolder) {
	    $out = "$outfolder$fn.fq";
	}
        if (defined $cleared{$fn}) {
            unless ( open(META, ">>$out") ) {
                print "Cannot open file \"$out\" to write to!!\n\n";
                exit;
            }
            print META "$result\n";
            close (META);                        
        } else {
            unless ( open(META, ">$out") ) {
                print "Cannot open file \"$out\" to write to!!\n\n";
                exit;
            }
            print META "$result\n";
            close (META);
            $cleared{$fn} = 1;                        
        }
        (@{$HoAgoods{$fn}}) = ();
    }
}

print "Out of $seenlines lines, $goodlines (with at least $goodthresh good sites) were retained.\n";
