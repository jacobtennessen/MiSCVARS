use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_g $opt_q $opt_c $opt_m $opt_o);

# Usage
my $usage = "
MakeOnemapFormatFromVcfNoParentsTwoGenos.pl - reads a vcf file and converts it to OneMap format, only for sites inferred to be heterozygous in one parent
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

Usage: perl MakeOnemapFormatFromVcfNoParentsTwoGenos.pl options
 required:
  -g	a vcf file. All individuals are assumed to be offspring
 optional:
  -q	miniumum quality score for a genotype to be called [default = 15]
  -c    minimum number of times each of at least two genotypes must be seen (e.g. at least this many reference homozygotes and heterozygotes) [default = 6]
  -m    maximum proportion of sites with missing data allowed [default = 0.3]
  -o    output folder [default is same folder as vcf file]

";

#############

# command line processing.
getopts('g:q:c:m:o:');
die $usage unless ($opt_g);

my ($genotypes, $lowqual, $commonthresh, $missingthresh, $outfolder);

$genotypes = $opt_g if $opt_g;
$lowqual = $opt_q ? $opt_q : 15;
$commonthresh = $opt_c ? $opt_c : 6;
$missingthresh = $opt_m ? $opt_m : 0.3;
$outfolder = $opt_o if $opt_o;

my @vcffiledata = split /\//, $genotypes;

my $vcf = pop @vcffiledata;

my @vcfnamedata = split /\./, $vcf;

my $outputfile = "OneMap_no_parents_qual$lowqual"."_$vcfnamedata[0].txt";

if (defined $outfolder) {
    unless ($outfolder =~ /\/$/) {
        $outfolder = "$outfolder/";
    }
    $outputfile = "$outfolder$outputfile";
} else {
    if (defined $vcffiledata[0]) {
        my $dir = join "/", @vcffiledata;
        $outputfile = "$dir/$outputfile";
    }
}

my $model = "D1.10";

my @plcode = ("AA","AB","BB","AC","BC","CC","AD","BD","CD","DD","AE","BE","CE","DE","EE");

my $goodsites = 0;

my $namecount = 0;

my @resultarray;

open(IN, $genotypes) || die "can't open $genotypes\n";

while (<IN>) {
	my $line = $_;
	$line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    my $limit = scalar(@data) - 1;
    if ($data[0] =~ /CHROM/gi) {;
        my @fullnames = @data[9..$limit];
        foreach my $fn (@fullnames) {
            $namecount +=1;
        }
        next;
    } elsif ($line=~ /^#/) {
        next;
    }
    my $chromsite = "$data[0]_$data[1]";
    my @inds = @data[9..$limit];
    my %seenF1genos;
    my @tempgenos;
    my $indcount = 0;
    my $missingcount = 0;
    my %hets;
    my %homos;
    foreach my $i (@inds) {
	my @idata = split ":", $i;
	if (defined $idata[2]) { #SNP called
	    if ($idata[2] == 0) {
            push @tempgenos, "-";
            $missingcount +=1;
	    } else {
            my @pldata = split ",", $idata[1];
            my $low;
            for (my $pl = 0; $pl < (scalar(@pldata)); $pl++) {
                if (($pldata[$pl] > 0)&&($pldata[$pl] < $lowqual)) { #if there are ambiguous possible genotypes, it's missing data.
                $low = "-";
                } elsif ((!defined $low)&&($pldata[$pl] == 0)) {
                $low = $plcode[$pl];
                } elsif ((defined $low)&&($pldata[$pl] == 0)) { #if more than one genotype match perfectly, missing data
                $low = "-";
                }
            }
            unless (defined $low) {
                $low = "-";
            }
            my @sa = split "", $low;
            if (defined $sa[1]) {
                if (defined $seenF1genos{$low}) {
                    $seenF1genos{$low} += 1;
                } else {
                   $seenF1genos{$low} = 1; 
                }
                if ($sa[0] =~ /$sa[1]/) {
                    $homos{$low} = 1;
                 } else {
                    $hets{$low} = 1;
                }
                
            }
            if ($low =~ /-/) {
                $missingcount +=1;
            }
            push @tempgenos, $low;
	    }
	}
	$indcount +=1;
    }
    if (scalar(keys %seenF1genos) <= 1) { #if not more than one genotype seen in offspring, useless.
        next;
    }
    unless (scalar (keys %hets) == 1) { #reject more complex inheritance patterns, only allow one type of heterozygote
        next;
    }
    unless (scalar (keys %homos) == 1) { #reject more complex inheritance patterns, only allow one type of homozygote
        next;
    }
    if ($missingcount/$namecount > $missingthresh) { #too many offspring data are missing
        next;
    }
    my $common = 0;
    foreach my $ksfg (keys %seenF1genos) {
        if ($seenF1genos{$ksfg} >= $commonthresh) {
            $common +=1;
        }
    }
    unless ($common == 2) {
        next;
    }
    my $hetfail = 0;
    foreach my $h (keys %hets) {
        my @alleles = split "", $h;
        unless ((defined $homos{"$alleles[0]$alleles[0]"})||(defined $homos{"$alleles[1]$alleles[1]"})) {
            $hetfail = 1;
        }
    }
    if ($hetfail == 1) {
        next;
    }
    
    my @finalgenos;

    for (my $tg = 0; $tg < (scalar(@tempgenos)); $tg++) {
	    if (defined $homos{$tempgenos[$tg]}) {
            $tempgenos[$tg] = "a";
	    } elsif (defined $hets{$tempgenos[$tg]}) {
            $tempgenos[$tg] = "ab";
	    }
	    push @finalgenos, $tempgenos[$tg];
    }
    my $finalgenoline = join ",", @finalgenos;
    push @resultarray, "$chromsite\t$model\t$finalgenoline";
}

close (IN);

my @finalresultarray;

foreach my $mate (@resultarray) {
    push @finalresultarray, "*$mate";
    $goodsites +=1;
}

my $result = join "\n", @finalresultarray;

unless ( open(META, ">$outputfile") ) {
    print "Cannot open file \"$outputfile\" to write to!!\n\n";
    exit;
}

print META "$namecount $goodsites\n$result";

close (META);
