#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
# ---------------------------------------------------------------
# Author: James Robert White PhD
# Email: jwhite@respherabio.com
# ---------------------------------------------------------------
my %data     = ();
my %samples  = ();
# ---------------------------------------------------------------
my %YEAR     = ("Chaput"         => 2017,
                "Gopalakrishnan" => 2017,
                "Matson"         => 2018,
                "Frankel"        => 2017,
                "Routy NSCLC"    => 2017,
                "Routy RCC"      => 2017);
# ---------------------------------------------------------------
# load merged file across all studies
my @header   = ();
my %features = ();
my %meta     = ();
open IN, "../analysis/figure04/index-results.txt" or die;
while(<IN>){
  chomp($_);
  $_    =~ s/\"//g;
  my @A = split "\t", $_;

  # "SampleID"	"Study"	"CancerType"	"Response"	"Type"	"Value"
  # "ERR2213660"	"Routy NSCLC"	"NSCLC"	"R"	"F.prausnitzii"	16.2267166105393
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  my $sampid   = $A[0];
  my $study    = $A[1];
  my $response = $A[3];
  my $index    = $A[4];
  #next if ($index =~ /\.no\./);
  my $val      = $A[5];
  $data{$index}{$sampid} = $val;
  if ($index =~ /nr/){
    $data{$index}{$sampid} = -1*$val;
  }
  $meta{$sampid}{"Response"} = $response;
  $meta{$sampid}{"Study"}    = $study;
}
close IN;

# ------------------------------------------------
# for each index -- create input table
# and generate a forest plot ---
my $outDir = "../analysis/figure03";
`rm -r $outDir`;
`mkdir $outDir`;
`mkdir $outDir/tables`;
`mkdir $outDir/raw-tables`;
`mkdir $outDir/pdfs`;
`mkdir $outDir/rcodes`;

open REPORT, ">$outDir/meta-analysis-output.txt" or die;
print REPORT "Feature\tFixed SMD\tRand SMD\tFixed P\tRand P\tI2\tTau\n";
foreach my $f (sort keys %data){
  print "$f...\n";

  my %studies = ();
  # apply filtering criteria here
  foreach my $k1 (sort keys %meta){ # sample ID
    push @{$studies{$meta{$k1}{"Study"}}}, $k1;
  }

  open OUT2, ">$outDir/raw-tables/$f\.txt" or die;
  open OUT, ">$outDir/tables/$f\.txt" or die;
  print OUT "study\tyear\tn.e\tmean.e\tsd.e\tn.c\tmean.c\tsd.c\n";
  foreach my $study (sort keys %studies){ # for this study
    my @samples        = @{$studies{$study}}; # get the valid samples
    my @responder      = ();
    my @nonresponder   = ();
    foreach my $s (@samples){
      if ($meta{$s}{"Response"} eq "R"){
        push @responder, $data{$f}{$s}; # add the feature value
      }elsif($meta{$s}{"Response"} eq "NR"){
        push @nonresponder, $data{$f}{$s};
      }
    }
    next if (!defined($responder[1]) or !defined($nonresponder[1]));

    my ($responderavg,$respondersamplestddev,$respondersamplestderr,$respondern)             = meanstde(@responder);
    my ($nonresponderavg,$nonrespondersamplestddev,$nonrespondersamplestderr,$nonrespondern) = meanstde(@nonresponder);

    foreach my $it (@responder){
      print OUT2 "$study\tR\t$it\n";
    }
    foreach my $it (@nonresponder){
      print OUT2 "$study\tNR\t$it\n";
    }

    print OUT "$study\t$YEAR{$study}\t$respondern\t$responderavg\t$respondersamplestddev\t$nonrespondern\t$nonresponderavg\t$nonrespondersamplestddev\n";
  } # end of study loop
  close OUT;
  close OUT2;

  # plotting
  open  RIN, ">$outDir/rcodes/$f\.rcode.in" or die;
  print RIN
  "library(meta)
t1<-read.table(\"$outDir/tables/$f\.txt\",sep=\"\t\", header=TRUE)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=t1, sm=\"SMD\", studlab=paste(study,year))
summary(meta1)
print(c(\"efixedp\",meta1\$pval.fixed))
print(c(\"erandp\",meta1\$pval.random))

out1 <- paste(\"$outDir/pdfs\/\", \"$f\.pdf\", sep=\"\")

pdf(out1,paper=\"special\",width=16,height=5,bg=\"white\")
forest(meta1, col.diamond.fixed=\"#2980B9\", col.diamond.random=\"#C0392B\", lab.e=\"Responder\", lab.c=\"Non-responder\", leftcols=c(\"studlab\", \"n.c\", \"n.e\"), leftlabs=c(\"Study\", \"Non-responder\", \"Responder\"), xlim=c(-2,2))
dev.off()

\n";
  close RIN;
  system("R CMD BATCH $outDir/rcodes/$f\.rcode.in $outDir/rcodes/$f\.rcode.out");

  # pullElements
  # SMD             95%-CI     z p-value
  # Fixed effect model   -0.2687 [-0.5053; -0.0321] -2.23  0.0260
  # Random effects model -0.2748 [-0.5575;  0.0079] -1.91  0.0568
  #
  # Quantifying heterogeneity:
  # tau^2 = 0.0328; H = 1.17 [1.00; 1.81]; I^2 = 26.5% [0.0%; 69.4%]
  #
  # Test of heterogeneity:
  # Q d.f. p-value
  # 6.80  5  0.2357
  my $fixedp   = "NA";
  my $fixedsmd = "NA";
  my $randp    = "NA";
  my $randsmd  = "NA";
  my $I2       = "NA";
  my $tausq    = "NA";

  open IN2, "$outDir/rcodes/$f\.rcode.out";
  while(<IN2>){
    chomp($_);
    my @A = split " ", $_;
    if ($_ =~ /Fixed effect model/){
      $fixedsmd = $A[3];
      $fixedp   = $A[7];
    }
    if ($_ =~ /Random effects model/){
      $randsmd = $A[3];
      $randp   = $A[7];
    }
    if ($_ =~ /^tau/){
      $I2    = $A[10];
      $tausq = $A[2];
      $tausq =~ s/\;//g;
    }
    if ($_ =~ /\[1\]\ \"efixedp\"/){
      $fixedp   = $A[2];
    }
    if ($_ =~ /\[1\]\ \"erandp\"/){
      $randp    = $A[2];
    }
  }
  close IN2;

  print REPORT "$f\t$fixedsmd\t$randsmd\t$fixedp\t$randp\t$I2\t$tausq\n";

} # end of features
close REPORT;


# ---------------------------------------------------------
# SUBROUTINES
# ---------------------------------------------------------
sub meanstde
{
  my (@B)       = @_;

  if (!defined($B[1])){
    return("NA","NA","NA");
  }

  # BEGIN TRANSFORMATION
  my @tmpB = @B;
  for my $b (0..$#B){
  	my $p     = $B[$b];
  	# my $asinsqrt = asin(sqrt($p));
    # $tmpB[$b] = $asinsqrt;
    $tmpB[$b] = $p;
  }
  @B = @tmpB;
  # END OF TRANSFORMATION


  my $sum       = 0;
  my $sumofsqrs = 0;
  foreach my $s (@B){
    $sum       += $s;
    $sumofsqrs += ($s**2);
  }
  my $favg          = $sum/($#B+1);
  my $variance      = ($sumofsqrs - ($sum**2/($#B+1)))/$#B;
  my $samplestddev  = sqrt($variance);
  my $samplestderr  = $samplestddev/sqrt($#B+1);
  my $n             = ($#B+1);

  # formatting
  $favg         = sprintf("%3.8f", $favg);
  $samplestddev = sprintf("%3.8f", $samplestddev);
  $samplestderr = sprintf("%3.8f", $samplestderr);

  if ($samplestddev ne "NA" and $samplestddev == 0.000){
    $samplestddev = 0.001;
  }
  if ($samplestderr ne "NA" and $samplestderr == 0.000){
    $samplestderr = 0.001;
  }

  my @returnarray   = ($favg,$samplestddev,$samplestderr,$n);

  return(@returnarray);
}
