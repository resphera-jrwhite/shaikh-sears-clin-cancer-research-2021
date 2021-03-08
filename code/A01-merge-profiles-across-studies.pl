#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
# ---------------------------------------------------------
# Create a single merged dataset across all cohorts to be
# used for signal recapture and new analysis
# ---------------------------------------------------------
# GLOBAL VARS
# ---------------------------------------------------------
my %data    = ();
my %samples = ();
# ---------------------------------------------------------
# SET OUT DIR
# ---------------------------------------------------------
my $outDir = "../analysis/A01-merge-profiles-across-studies";
if (-e $outDir){
  `rm -r $outDir`;
}
`mkdir $outDir`;
# ---------------------------------------------------------
# FILE PATHS
# ---------------------------------------------------------
my $TOPDIR         = "../data/processed";
my $ROUTYSPECIES   = "$TOPDIR/Routy.2017/WGS/taxonomic-assignments/metaphlan2.species.txt";
my $FRANKELSPECIES = "$TOPDIR/Frankel.2017/WGS/taxonomic-assignments/metaphlan2.species.txt";
my $MATSONSPECIES  = "$TOPDIR/Matson.2018/16SrRNA/taxonomic-assignments/otu_table_even10000.txt";
my $CHAPUTSPECIES  = "$TOPDIR/Chaput.2017/16SrRNA/taxonomic-assignments/otu_table.species.txt";
my $GOPALAKRISHNANSPECIES = "$TOPDIR/Gopalakrishnan.2017/16SrRNA/taxonomic-assignments/otu_table_species.txt";
# ---------------------------------------------------------
# METAPHLAN2 WGS load
# ---------------------------------------------------------
my @header = ();
open IN, "$ROUTYSPECIES" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }
  # only bacteria / archaea
  next if ($A[0] =~ /(k__Eukaryota|k__Viruses)/);
  my @a0 = split /\|/, $A[0];
  my $species = $a0[$#a0];
  $species =~ s/^s__//g;
  for my $i (1 .. $#A){
    $data{$species}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$FRANKELSPECIES" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  # only bacteria / archaea
  next if ($A[0] =~ /(k__Eukaryota|k__Viruses)/);
  my @a0 = split /\|/, $A[0];
  my $species = $a0[$#a0];
  $species =~ s/^s__//g;
  for my $i (1 .. $#A){
    $data{$species}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$MATSONSPECIES" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  next if ($_ =~ /\#\ Constructed\ from\ biom\ file/);
  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $species = $A[0];
  for my $i (1 .. ($#A-1)){
    $data{$species}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$CHAPUTSPECIES" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  next if ($_ =~ /\#\ Constructed\ from\ biom\ file/);

  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $species = $A[0];
  for my $i (1 .. ($#A-1)){
    $data{$species}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$GOPALAKRISHNANSPECIES" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  next if ($_ =~ /\#\ Constructed\ from\ biom\ file/);
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  my $species = $A[0];

  for my $i (1 .. ($#A-1)){
    $data{$species}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

# zero out all unobserved entries
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    if (!defined($data{$k1}{$k2})){
      $data{$k1}{$k2} = 0;
    }
  }
}

# normalize
my %totals = ();
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $totals{$k2} += $data{$k1}{$k2};
  }
}

foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $data{$k1}{$k2} = 100*$data{$k1}{$k2}/$totals{$k2};
  }
}

open OUT, ">$outDir/out.A01-species.txt" or die;
foreach my $k2 (sort keys %samples){
  print OUT "\t$k2";
}
print OUT "\n";

foreach my $k1 (sort keys %data){
  print OUT "$k1";
  foreach my $k2 (sort keys %samples){
    print OUT "\t$data{$k1}{$k2}";
  }
  print OUT "\n";
}
close OUT;


# now for genus level ----
my $ROUTYGEN   = "$TOPDIR/Routy.2017/WGS/taxonomic-assignments/metaphlan2.genus.txt";
my $FRANKELGEN = "$TOPDIR/Frankel.2017/WGS/taxonomic-assignments/metaphlan2.genus.txt";
my $MATSONGEN  = "$TOPDIR/Matson.2018/16SrRNA/taxonomic-assignments/summarize_taxa/otu_table_even10000_L6.txt";
my $CHAPUTGEN  = "$TOPDIR/Chaput.2017/16SrRNA/taxonomic-assignments/summarize_taxa/otu_table_L6.txt";
my $GOPALAKRISHNANGEN = "$TOPDIR/Gopalakrishnan.2017/16SrRNA/taxonomic-assignments/summarize_taxa/otu_table_L6.txt";
# reset variables ---
%data    = ();
%samples = ();
# ---------------------------------------------------------
# METAPHLAN2 WGS load
# ---------------------------------------------------------
@header = ();
open IN, "$ROUTYGEN" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }
  # only bacteria / archaea
  next if ($A[0] =~ /(k__Eukaryota|k__Viruses)/);
  my @a0 = split /\|/, $A[0];
  my $GEN = $a0[$#a0];
  $GEN =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$GEN}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$FRANKELGEN" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  # only bacteria / archaea
  next if ($A[0] =~ /(k__Eukaryota|k__Viruses)/);
  my @a0 = split /\|/, $A[0];
  my $GEN = $a0[$#a0];
  $GEN =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$GEN}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$MATSONGEN" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $GEN = $A[0];
  my @GEN = split /\;/, $GEN;
  $GEN    = $GEN[$#GEN];
  $GEN =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$GEN}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$CHAPUTGEN" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $GEN = $A[0];
  my @GEN = split /\;/, $GEN;
  $GEN    = $GEN[$#GEN];
  $GEN =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$GEN}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$GOPALAKRISHNANGEN" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    for my $j (1..$#header){
      $header[$j] =~ s/\.R1//g;
    }
    next;
  }
  my $GEN = $A[0];
  my @GEN = split /\;/, $GEN;
  $GEN    = $GEN[$#GEN];
  $GEN =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$GEN}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

# zero out all unobserved entries
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    if (!defined($data{$k1}{$k2})){
      $data{$k1}{$k2} = 0;
    }
  }
}

# normalize
%totals = ();
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $totals{$k2} += $data{$k1}{$k2};
  }
}

foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $data{$k1}{$k2} = 100*$data{$k1}{$k2}/$totals{$k2};
  }
}

open OUT, ">$outDir/out.A01-genus.txt" or die;
foreach my $k2 (sort keys %samples){
  print OUT "\t$k2";
}
print OUT "\n";

foreach my $k1 (sort keys %data){
  print OUT "$k1";
  foreach my $k2 (sort keys %samples){
    print OUT "\t$data{$k1}{$k2}";
  }
  print OUT "\n";
}
close OUT;


# now for family level ----
my $ROUTYFAM   = "$TOPDIR/Routy.2017/WGS/taxonomic-assignments/metaphlan2.family.txt";
my $FRANKELFAM = "$TOPDIR/Frankel.2017/WGS/taxonomic-assignments/metaphlan2.family.txt";
my $MATSONFAM  = "$TOPDIR/Matson.2018/16SrRNA/taxonomic-assignments/summarize_taxa/otu_table_even10000_L5.txt";
my $CHAPUTFAM  = "$TOPDIR/Chaput.2017/16SrRNA/taxonomic-assignments/summarize_taxa/otu_table_L5.txt";
my $GOPALAKRISHNANFAM = "$TOPDIR/Gopalakrishnan.2017/16SrRNA/taxonomic-assignments/summarize_taxa/otu_table_L5.txt";
# reset variables ---
%data    = ();
%samples = ();
# ---------------------------------------------------------
# METAPHLAN2 WGS load
# ---------------------------------------------------------
@header = ();
open IN, "$ROUTYFAM" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }
  # only bacteria / archaea
  next if ($A[0] =~ /(k__Eukaryota|k__Viruses)/);
  my @a0 = split /\|/, $A[0];
  my $FAM = $a0[$#a0];
  $FAM =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$FAM}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$FRANKELFAM" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  # only bacteria / archaea
  next if ($A[0] =~ /(k__Eukaryota|k__Viruses)/);
  my @a0 = split /\|/, $A[0];
  my $FAM = $a0[$#a0];
  $FAM =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$FAM}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$MATSONFAM" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $FAM = $A[0];
  my @FAM = split /\;/, $FAM;
  $FAM    = $FAM[$#FAM];
  $FAM =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$FAM}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$CHAPUTFAM" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $FAM = $A[0];
  my @FAM = split /\;/, $FAM;
  $FAM    = $FAM[$#FAM];
  $FAM =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$FAM}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$GOPALAKRISHNANFAM" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    for my $j (1..$#header){
      $header[$j] =~ s/\.R1//g;
    }
    next;
  }
  my $FAM = $A[0];
  my @FAM = split /\;/, $FAM;
  $FAM    = $FAM[$#FAM];
  $FAM =~ s/^g__/g\./g;
  for my $i (1 .. $#A){
    $data{$FAM}{$header[$i]} += $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

# zero out all unobserved entries
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    if (!defined($data{$k1}{$k2})){
      $data{$k1}{$k2} = 0;
    }
  }
}

# normalize
%totals = ();
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $totals{$k2} += $data{$k1}{$k2};
  }
}

foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $data{$k1}{$k2} = 100*$data{$k1}{$k2}/$totals{$k2};
  }
}

open OUT, ">$outDir/out.A01-family.txt" or die;
foreach my $k2 (sort keys %samples){
  print OUT "\t$k2";
}
print OUT "\n";

foreach my $k1 (sort keys %data){
  print OUT "$k1";
  foreach my $k2 (sort keys %samples){
    print OUT "\t$data{$k1}{$k2}";
  }
  print OUT "\n";
}
close OUT;


# ---------------------------------------------------------
# add functional summaries
# ---------------------------------------------------------
# FILE PATHS
# ---------------------------------------------------------
my $ROUTYFUNCTIONAL          = "$TOPDIR/Routy.2017/16SrRNA/functional-inference/ko_kegglevel3.txt";
my $FRANKELFUNCTIONAL        = "$TOPDIR/Frankel.2017/16SrRNA/functional-inference/ko_kegglevel3.txt";
my $MATSONFUNCTIONAL         = "$TOPDIR/Matson.2018/16SrRNA/functional-inference/ko_kegglevel3.txt";
my $CHAPUTFUNCTIONAL         = "$TOPDIR/Chaput.2017/16SrRNA/functional-inference/ko_kegglevel3.txt";
my $GOPALAKRISHNANFUNCTIONAL = "$TOPDIR/Gopalakrishnan.2017/16SrRNA/functional-inference/ko_kegglevel3.txt";
# reset variables ---
%data    = ();
%samples = ();
# ---------------------------------------------------------
# PICRUSt level 3 kegg assignments
# ---------------------------------------------------------
@header = ();
open IN, "$ROUTYFUNCTIONAL" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }
  my $FUNCTIONAL = $A[0];
  for my $i (1 .. $#A){
    $data{$FUNCTIONAL}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$FRANKELFUNCTIONAL" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    for my $i (1 .. $#header){
      $header[$i] =~ s/\.R1//g;
    }
    next;
  }
  my $FUNCTIONAL = $A[0];
  for my $i (1 .. $#A){
    $data{$FUNCTIONAL}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$MATSONFUNCTIONAL" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  my $FUNCTIONAL = $A[0];
  for my $i (1 .. $#A){
    $data{$FUNCTIONAL}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$CHAPUTFUNCTIONAL" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;

  if (!defined($header[1])){
    @header = @A;
    next;
  }

  my $FUNCTIONAL = $A[0];
  for my $i (1 .. $#A){
    $data{$FUNCTIONAL}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

@header = ();
open IN, "$GOPALAKRISHNANFUNCTIONAL" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    for my $i (1 .. $#header){
      $header[$i] =~ s/\.R1//g;
    }
    next;
  }
  my $FUNCTIONAL = $A[0];
  for my $i (1 .. $#A){
    $data{$FUNCTIONAL}{$header[$i]} = $A[$i];
    $samples{$header[$i]} = 1;
  }
}
close IN;

# zero out all unobserved entries
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    if (!defined($data{$k1}{$k2})){
      $data{$k1}{$k2} = 0;
    }
  }
}

# normalize
%totals = ();
foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $totals{$k2} += $data{$k1}{$k2};
  }
}

foreach my $k1 (sort keys %data){
  foreach my $k2 (sort keys %samples){
    $data{$k1}{$k2} = 100*$data{$k1}{$k2}/$totals{$k2};
  }
}

open OUT, ">$outDir/out.A01-picrust.txt" or die;
foreach my $k2 (sort keys %samples){
  print OUT "\t$k2";
}
print OUT "\n";

foreach my $k1 (sort keys %data){
  print OUT "$k1";
  foreach my $k2 (sort keys %samples){
    print OUT "\t$data{$k1}{$k2}";
  }
  print OUT "\n";
}
close OUT;
