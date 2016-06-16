use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_g $opt_q $opt_c $opt_m $opt_o $opt_h);

# Usage
my $usage = "
MakeOnemapFormatFromVcfNoParentsGeneralFormat.pl - reads a vcf file with any genotype format and converts it to OneMap format, all sites designted as heterozygous in the same parent, though other inheritance patterns may be present
Copyright (C) 2016 by Jacob A Tennessen

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

Usage: perl MakeOnemapFormatFromVcfNoParentsGeneralFormat.pl options
 required:
  -g	a vcf file. All individuals are assumed to be offspring
 optional:
  -q	miniumum quality score for a genotype to be called [default = 15]
  -c    minimum number of times each of at least two genotypes must be seen (e.g. at least this many reference homozygotes and heterozygotes) [default = 6]
  -m    maximum proportion of sites with missing data allowed [default = 0.3]
  -o    output folder [default is same folder as vcf file]
  -h    maximum number of different homozygous genotypes (e.g. use 1 to exclude markers heterozygous in both parents) [default = 2]

";

#############

# command line processing.
getopts('g:q:c:m:o:h:');
die $usage unless ($opt_g);

my ($genotypes, $lowqual, $commonthresh, $missingthresh, $outfolder, $maxhomos);

$genotypes = $opt_g if $opt_g;
$lowqual = $opt_q ? $opt_q : 15;
$commonthresh = $opt_c ? $opt_c : 6;
$missingthresh = $opt_m ? $opt_m : 0.3;
$outfolder = $opt_o if $opt_o;
$maxhomos = $opt_h ? $opt_h : 2;

my @vcffiledata = split /\//, $genotypes;

my $vcf = pop @vcffiledata;

my @vcfnamedata = split /\./, $vcf;

my $suffix = pop @vcfnamedata;

my $shortname = join "\.", @vcfnamedata;

my $outputfile = "OneMap_no_parents_qual$lowqual"."_$shortname.txt";

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
    my $dppos;
    my $plpos;
    my @genoformat = split ":", $data[8];
    for (my $gf = 0; $gf < (scalar(@genoformat)); $gf++) {
        if ($genoformat[$gf] =~ /DP/) {
            $dppos = $gf;
        } elsif ($genoformat[$gf] =~ /PL/) {
            $plpos = $gf;
        }
    }
    unless ((defined $dppos)&&(defined $plpos)) {
        next;
    }
    foreach my $i (@inds) {
	my @idata = split ":", $i;
	if (defined $idata[$dppos]) { #SNP called
	    if (($idata[$dppos] !~ /\d/)||($idata[$dppos] == 0)||($idata[$plpos] !~ /\d/)) {
            push @tempgenos, "-";
            $missingcount +=1;
	    } else {
            my @pldata = split ",", $idata[$plpos];
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
                unless ($sa[0] =~ /$sa[1]/) {
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
    if ($missingcount/$namecount > $missingthresh) { #too many offspring data are missing
        next;
    }
    my $common = 0;
    foreach my $ksfg (keys %seenF1genos) {
        if ($seenF1genos{$ksfg} >= $commonthresh) {
            $common +=1;
        }
    }
    unless ($common >= 2) {
        next;
    }
    my $homo1;
    my $homo2;
    foreach my $h (keys %hets) {
        my @alleles = split "", $h;
        $homo1 = "$alleles[0]$alleles[0]";
        $homo2 = "$alleles[1]$alleles[1]";
    }
    my @finalgenos;
    
    my %seenhomos;

    for (my $tg = 0; $tg < (scalar(@tempgenos)); $tg++) {
	    if ($tempgenos[$tg] =~ /$homo1/) {
            $tempgenos[$tg] = "a";
            $seenhomos{"a"} = 1;
	    } elsif (defined $hets{$tempgenos[$tg]}) {
            $tempgenos[$tg] = "ab";
	    } elsif ($tempgenos[$tg] =~ /$homo2/) {
            $tempgenos[$tg] = "b";
            $seenhomos{"b"} = 1;
	    }
	    push @finalgenos, $tempgenos[$tg];
    }
    if ((scalar(keys %seenhomos)) > $maxhomos) {
        next;
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
