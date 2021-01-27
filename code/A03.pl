#!/usr/bin/perl
use Data::Dumper;
use Getopt::Std;
use Scalar::Util qw/looks_like_number/;
use List::Util qw/shuffle sum min max/;
use POSIX qw/ceil floor/;
use File::Basename;
use warnings;
use strict;
# --------------------------------------------------------------------
#  ______.pl
#  Author: James Robert White, PhD
#  Email: jwhite@respherabio.com
# --------------------------------------------------------------------
#  description of program
# --------------------------------------------------------------------
use vars qw/$opt_i/;
getopts("i");
my $usage =
".USAGE.
./_______ -i < input >

.DESCRIPTION.

.OPTIONS.
  -i input\
\n";

die $usage unless defined $opt_i;
# --------------------------------------------------------------------
my $outDir = "../analysis/".basename($0); $outDir =~ s/\.pl$//g;
if (-e $outDir){
  `rm -r $outDir`;
}
`mkdir $outDir`;
# --------------------------------------------------------------------
