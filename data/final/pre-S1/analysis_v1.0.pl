#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# create a single merged dataset of consensus species level
# calls
my %wgsdata       = ();
my %wgssamples    = ();
my %amp16Sdata    = ();
my %amp16Ssamples = ();
my $ROOTDIR       = "....";

my @studyDirs = qw/Gopalakrishnan.2017 Matson.2018 Frankel.2017/;

foreach my $studyDir (@studyDirs){
  # load metaphlan species level
  my $f      = "$ROOTDIR/$studyDir/WGS/taxonomic-assignments/metaphlan2.species.txt";
  my @header = ();
  my %wgscountdata = ();
  my %wgstotals    = ();
  open IN, "$f" or die "Error cannot find $f\n";
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
      next;
    }

    if ($A[0] !~ /k__(Archaea|Bacteria)/){
      next;
    }
    my @spID      = split /\|/, $A[0];
    my $speciesID = $spID[6];
    for my $i (1 .. $#A){
      $wgscountdata{$speciesID}{$header[$i]} = $A[$i];
      $wgstotals{$header[$i]} += $A[$i];
      $wgssamples{$header[$i]} = 1;
    }
  }
  close IN;
  # normalize to prokaryotic % abundances
  foreach my $k1 (sort keys %wgscountdata){
    foreach my $k2 (sort keys %wgstotals){
      $wgsdata{$k1}{$k2} = 100*$wgscountdata{$k1}{$k2}/$wgstotals{$k2};
    }
  }

  # load 16S species assignments
  my $f2 = "$ROOTDIR/$studyDir/16SrRNA/taxonomic-assignments/otu_table_even10000.txt";
  if ($studyDir eq "Gopalakrishnan.2017"){
     $f2 = "$ROOTDIR/$studyDir/16SrRNA/taxonomic-assignments/otu_table_even8000.txt";
  }
  my @header2         = ();
  my %amp16Scountdata = ();
  my %amp16Stotals    = ();
  open IN, "$f2" or die "Error cannot find $f2\n";
  while(<IN>){
    chomp($_);
    next if ($_ =~ /\#\ Constructed\ from\ biom\ file/);
    my @A = split "\t", $_;
    if (!defined($header2[1])){
      if ($studyDir !~ /Frankel/){
        for my $i (0 .. $#A){
          $A[$i] =~ s/\.R1//g;
        }
      }
      @header2 = @A;
      next;
    }

    for my $i (1 .. ($#A-1)){
      $amp16Scountdata{$A[0]}{$header2[$i]} = $A[$i];
      $amp16Stotals{$header2[$i]} += $A[$i];
      $amp16Ssamples{$header2[$i]} = 1;
    }
  }
  close IN;

  foreach my $k1 (sort keys %amp16Scountdata){
    foreach my $k2 (sort keys %amp16Stotals){
      $amp16Sdata{$k1}{$k2} = 100*$amp16Scountdata{$k1}{$k2}/$amp16Stotals{$k2};
    }
  }
} # end of studyDirs loop

open OUT, ">wgs.prct.txt" or die;
print OUT "species";
foreach my $k2 (sort keys %wgssamples){
  print OUT "\t$k2";
}
print OUT "\n";
foreach my $k1 (sort keys %wgsdata){
  print OUT "$k1";
  foreach my $k2 (sort keys %wgssamples){
    if (!defined($wgsdata{$k1}{$k2})){
      $wgsdata{$k1}{$k2} = 0;
    }
    print OUT "\t$wgsdata{$k1}{$k2}";
  }
  print OUT "\n";
}
close OUT;

open OUT, ">16SrRNA.prct.txt" or die;
print OUT "species";
foreach my $k2 (sort keys %amp16Ssamples){
  print OUT "\t$k2";
}
print OUT "\n";
foreach my $k1 (sort keys %amp16Sdata){
  print OUT "$k1";
  foreach my $k2 (sort keys %amp16Ssamples){
    if (!defined($amp16Sdata{$k1}{$k2})){
      $amp16Sdata{$k1}{$k2} = 0;
    }
    print OUT "\t$amp16Sdata{$k1}{$k2}";
  }
  print OUT "\n";
}
close OUT;

my %pairmeta = ();
open IN, "metadata.txt" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  next if ($A[0] eq "PID");
  $pairmeta{$A[0]}{"16SrRNA"} = $A[1];
  $pairmeta{$A[0]}{"WGS"}     = $A[2];
  $pairmeta{$A[0]}{"Study"}   = $A[3];
}
close IN;

my %targets = ();
open IN, "targets.txt" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  $targets{$A[0]} = $A[1];
}
close IN;

open OUT, ">results.txt" or die;
print OUT "PID\tStudy\tSpecies\tWGSID\tabWGS\t16SID\tab16SrRNA\n";
foreach my $k1 (sort keys %pairmeta){
  my $amp16Sid = $pairmeta{$k1}{"16SrRNA"};
  my $wgsid    = $pairmeta{$k1}{"WGS"};
  my $studyid  = $pairmeta{$k1}{"Study"};

  foreach my $t1 (sort keys %targets){
    print OUT "$k1\t$studyid\t$t1\t$wgsid\t$wgsdata{$t1}{$wgsid}\t$amp16Sid\t$amp16Sdata{$targets{$t1}}{$amp16Sid}\n";
  }
}
close OUT;
