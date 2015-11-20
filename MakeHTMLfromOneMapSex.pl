use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_s $opt_p $opt_r);

# Usage
my $usage = "
MakeHTMLfromOneMapSex.pl- reads OneMap format file and coverts to color-coded HTML for easy viewing
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

Usage: perl MakeHTMLfromOneMapSex.pl options
 required:
  -f	a OneMap file (or equivalent - needs three columns, with comma-delimited genotypes in thrid column)
  -s	a file listing all samples names in order as they appear in onemap file (whitespace delimited)
 optional:
  -p 	phenotype file
  -r	reference position from which colors are determined (should be a non-recombinant) [default = 5]
";

#############

# command line processing.
getopts('f:s:p:r:');
die $usage unless ($opt_f);
die $usage unless ($opt_s);

my ($bigfile, $samplenames, $phenos, $refpos);

$bigfile = $opt_f if $opt_f;
$samplenames = $opt_s if $opt_s;
$phenos = $opt_p if $opt_p;
$refpos = $opt_r ? $opt_r : 5;

my $usemp = 0;

unless (open(SAMP, $samplenames)) {
  print "Cannot open file \"$samplenames\"\n\n";
  exit;
}

my @samples;

while (<SAMP>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split /\s+/, $line;
  push @samples, @data;
}

close (SAMP);

my @malesexes;
my @femalesexes;

if ((defined $phenos)&&(open(SEX, $phenos))) {

  my %malesexes;
  my %HoAfemalesexes;
  
  while (<SEX>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    if ($line =~ /Plant/) {
      next;
    }
    my @data = split "\t", $line;
    if (defined $malesexes{$data[0]}) {
      if ($data[1] =~ /N/) {
	if ($malesexes{$data[0]} =~ /U/) {
	  $malesexes{$data[0]} = "S";
	} else {
	  unless ($malesexes{$data[0]} =~ /S/) {
	    print "Male sex conflict for $data[0]";
	    exit;
	  }
	}
      } elsif ($data[1] =~ /Y/) {
	if ($malesexes{$data[0]} =~ /U/) {
	  $malesexes{$data[0]} = "M";
	} else {
	  unless ($malesexes{$data[0]} =~ /M/) {
	    print "Male sex conflict for $data[0]";
	    exit;
	  }
	}
      }
    } else {
      if ($data[1] =~ /N/) {
	$malesexes{$data[0]} = "S";
      } elsif ($data[1] =~ /Y/) {
	$malesexes{$data[0]} = "M";
      } else {
	$malesexes{$data[0]} = "U";
      }
    }
    push @{$HoAfemalesexes{$data[0]}}, $data[2];
  }
  
  close (SEX);
  
  my %femalesexes;
  
  foreach my $f (keys %HoAfemalesexes) {
    my $total = 0;
    my $count = 0;
    foreach my $v (@{$HoAfemalesexes{$f}}) {
      if ($v =~ /\d/) {
	$total += $v;
	$count +=1;
      }
    }
    if ($count > 0) {
      my $mean = $total/$count;
      if ($mean >= 50) {
	$femalesexes{$f} = "F";
      } else {
	$femalesexes{$f} = "S";
      }
    } else {
      $femalesexes{$f} = "U";
    }
  }
  
  foreach my $s (@samples) {
    unless (defined $malesexes{$s}) {
      $malesexes{$s} = "-";
    }
    unless (defined $femalesexes{$s}) {
      $femalesexes{$s} = "-";
    }
    push @malesexes, $malesexes{$s};
    push @femalesexes, $femalesexes{$s};
  }

}

unless (open(BIG, $bigfile)) {
  print "Cannot open file \"$bigfile\"\n\n";
  exit;
}

my $samplesize;

my @out;

push @out, "<html>";
push @out, "<head>";
push @out, '  <meta content="text/html; charset=ISO-8859-1';
push @out, ' http-equiv="content-type">';
push @out, "  <title>$bigfile</title>";
push @out, "</head>";
push @out, "<body>";
push @out, "<span style=\"font-family: Courier New;\">";

my @blues;
my @yellows;
my @hets;
my @homos;
my @incompatible;
my $opentable = 0;
my $linecount = 0;

while (<BIG>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  unless (($line =~ /Fvb/gi)||($line =~ /BAC/gi)) {
    @incompatible = ();
    if ($opentable == 1) {
      for (my $z = 0; $z < $samplesize; $z++) {
	if (($hets[$z] > 0)&&($homos[$z] > 0)) {
	  if (($yellows[$z] > 0)&&($blues[$z] > 0)) {
	    if (($yellows[$z] > 1)&&($blues[$z] > 1)) {
	      $incompatible[$z] = "*";
	    } else {
	      $incompatible[$z] = "e";
	    }
	  } else {
	    $incompatible[$z] = "&nbsp";
	  }
	} else {
	  $incompatible[$z] = "i";
	}
	$blues[$z] = 0;
	$yellows[$z] = 0;
	$hets[$z] = 0;
	$homos[$z] = 0;
      }
      my $ilist = join "", @incompatible;
      push @out, "<tr>\n<td></td>\n<td></td>\n<td>$ilist</td>\n</tr>";
      push @out, "</table><br>$line<br>".'<table style="width:100%">';
    } else {
      push @out, "$line<br>".'<table style="width:100%">';
    }
    $opentable = 1;
    $linecount = 0;
    next;
  }
  my @gdata = split "\t", $line;
  if (($usemp == 0)&&($gdata[0] =~ /mp$/)) {
    next;
  }
  $gdata[2] =~ s/ab/h/g;
  my @genos = split ",", $gdata[2];
  unless (defined $samplesize) {
    $samplesize = scalar(@genos);
    unless ((scalar(@samples)) == $samplesize) {
      my $wrong = (scalar(@samples));
      print "$wrong doesn't match $samplesize.\n";
      exit;
    }
    my %HoAoutnumbers;
    for (my $x = 1; $x<= $samplesize; $x++) {
      my @digits = split "", $x;
      for (my $d = 3; $d >= 1; $d--) {
	if (defined $digits[$d-1]) {
	  push @{$HoAoutnumbers{$d}}, $digits[$d-1];
	} else {
	  push @{$HoAoutnumbers{$d}}, "-"
	}
      }
      push @yellows, 0;
      push @blues, 0;
      push @hets, 0;
      push @homos, 0;
    }
    for (my $d = 1; $d <= 3; $d++) {
      my $numblist = join "", @{$HoAoutnumbers{$d}};
      push @out, "<tr>\n<td></td>\n<td></td>\n<td>$numblist</td>\n</tr>";
    }
    if (defined $femalesexes[0]) {
      my $femalesexlist = join "", @femalesexes;
      push @out, "<tr>\n<td>Female_Fert</td>\n<td></td>\n<td>$femalesexlist</td>\n</tr>";
    }
    if (defined $malesexes[0]) {
      my $malesexlist = join "", @malesexes;
      push @out, "<tr>\n<td>Male_Fert</td>\n<td></td>\n<td>$malesexlist</td>\n</tr>";
    }
  }
  my $refgeno;
  my $localrefpos = $refpos;
  until (defined $refgeno) {
    unless (defined $genos[$localrefpos]) {
      $localrefpos = 2;
    }
    if ((($genos[$localrefpos] =~ /a/)||($genos[$localrefpos] =~ /h/))&&($blues[$localrefpos] == 0)&&(($localrefpos == $refpos)||($linecount == 0)||($yellows[$localrefpos] > 0))) {
      $refgeno = $genos[$localrefpos];
    } else {
      $localrefpos +=1;
    }
    if (($localrefpos + 1) == $refpos) {
      $refgeno = $genos[$refpos];
    }
  }
  my $gcount = 0;
  my $lastcolor = "n";
  my $spanopen = 0;
  foreach my $g (@genos) {
    if ($g =~ /h/) {
      $hets[$gcount] +=1;
    } elsif ($g =~ /a/i) {
      $homos[$gcount] +=1;
    }
    if ($g =~ /$refgeno/i) {
      if ($gdata[0] =~ /mp$/) {
	unless ($lastcolor =~ /y/) {
	  if ($spanopen == 1) {
	    $g = '</span><span style="background-color: rgb(140, 198, 63);">'."$g";
	  } else {
	    $g = '<span style="background-color: rgb(140, 198, 63);">'."$g";
	  }
	}
      } else {
	unless ($lastcolor =~ /y/) {
	  if ($spanopen == 1) {
	    $g = '</span><span style="background-color: rgb(255, 255, 0);">'."$g";
	  } else {
	    $g = '<span style="background-color: rgb(255, 255, 0);">'."$g";
	  }
	}
	$yellows[$gcount] += 1;
      }
      $spanopen = 1;
      $lastcolor = "y";
    } elsif ($g =~ /-/) {
      if ($spanopen == 1) {
	$g = "</span>m";
      } else {
	$g = "m";
      }
      $lastcolor = "n";
      $spanopen = 0;
    } else {
      if ($gdata[0] =~ /mp$/) {
	unless ($lastcolor =~ /b/) {
	  if ($spanopen == 1) {
	    $g = '</span><span style="background-color: rgb(145, 39, 143);">'."$g";
	  } else {
	    $g = '<span style="background-color: rgb(145, 39, 143);">'."$g";
	  }
	}
      } else {
	unless ($lastcolor =~ /b/) {
	  if ($spanopen == 1) {
	    $g = '</span><span style="background-color: rgb(153, 255, 255);">'."$g";
	  } else {
	    $g = '<span style="background-color: rgb(153, 255, 255);">'."$g";
	  }
	}
	$blues[$gcount] += 1;
      }
      $spanopen = 1;
      $lastcolor = "b";
    }
    $gcount +=1;
  }
  my $newgenos = join "", @genos;
  if ($spanopen == 1) {
    $newgenos = "$newgenos</span>";
  }
  my $newline = "<tr>\n<td>$gdata[0]</td>\n<td>$gdata[1]</td>\n<td>$newgenos</td>\n</tr>";
  push @out, $newline;
  $linecount +=1;
}

close (BIG);

@incompatible = ();

for (my $z = 0; $z < $samplesize; $z++) {
  if (($hets[$z] > 0)&&($homos[$z] > 0)) {
    if (($yellows[$z] > 0)&&($blues[$z] > 0)) {
      if (($yellows[$z] > 1)&&($blues[$z] > 1)) {
	$incompatible[$z] = "*";
      } else {
	$incompatible[$z] = "e";
      }
    } else {
      $incompatible[$z] = "&nbsp";
    }
  } else {
    $incompatible[$z] = "i";
  }
}

my $ilist = join "", @incompatible;

push @out, "<tr>\n<td></td>\n<td></td>\n<td>$ilist</td>\n</tr>";

push @out, "</table>";
push @out, "</span>";
push @out, "</body>";
push @out, "</html>";

my $result = join "\n", @out;

my $outfile = "$bigfile.html";

unless ( open(HWV, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}

print HWV $result;

close (HWV);
