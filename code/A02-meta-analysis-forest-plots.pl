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

# ---------------------------------------------------------
# SET OUT DIR
# ---------------------------------------------------------
my $outDir = "../analysis/A02-meta-analysis-forest-plots";
if (-e $outDir){
  `rm -r $outDir`;
}
`mkdir $outDir`;

# expanding meta-analysis results to other feature types:
my @featureTypes = qw/picrust species/;

foreach my $ft (@featureTypes){
  `mkdir $outDir/$ft`;
  # ---------------------------------------------------------------
  # load merged file across all studies
  my @header   = ();
  my %features = ();
  open IN, "../analysis/A01-merge-profiles-across-studies/out.A01-$ft.txt" or die "Cannot locate out.A01-$ft.txt\n";
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;

    if (!defined($header[1])){
      @header = @A;
      next;
    }

    my $nonZeros = 0;
    for my $i (1 .. $#A){
      $data{$A[0]}{$header[$i]} = $A[$i];
      $samples{$header[$i]} = 1;
      if ($A[$i] > 0){
        $nonZeros++;
      }
    }
    if ($nonZeros >= 10){
      $features{$A[0]} = 1;
    }
  }
  close IN;

  my %meta = ();
  @header  = ();
  open IN, "../data/final/metadata.2019-04-22.txt" or die;
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
      next;
    }

    for my $i (1 .. $#A){
      $meta{$A[0]}{$header[$i]} = $A[$i];
    }
  }
  close IN;
  # ---------------------------------------------------
  # for each taxon -- create input table
  # and generate a forest plot
  `rm -r $outDir/$ft/tables $outDir/$ft/pdfs $outDir/$ft/rcodes`;
  `mkdir $outDir/$ft/tables`;
  `mkdir $outDir/$ft/raw-tables`;
  `mkdir $outDir/$ft/pdfs`;
  `mkdir $outDir/$ft/rcodes`;

  my @studies = ("Chaput", "Frankel", "Gopalakrishnan", "Matson",
  "Routy NSCLC", "Routy RCC");

  open REPORT, ">$outDir/$ft/out.$ft.meta-analysis.report.txt" or die;
  print REPORT "Feature\tInclusion\tComparator\tFixed SMD\tRand SMD\tFixed P\tRand P\tI2\tTau";
  foreach my $study (@studies){
    print REPORT "\t$study";
  }
  print REPORT "\n";
  foreach my $keepi (qw/UseforMeta/){
  foreach my $responsei (qw/Response/){
    my $fdex = 1;
    foreach my $f (sort keys %features){
      if ($f =~ /^otu/){
        next;
      }

      my $numColon = $f =~ tr/\://;
      if ($numColon > 2){
        next;
      }

      my $fname = $f;
      if (length($f) > 75){
        $fname = substr($f, 0, 75);
      }
      $fname =~ s/^f\__/f\./g;

      if ($ft eq "picrust"){
        $fname =~ s/\ /\./g;
        $fname =~ s/[\;\:\)\(\,\'\"\-\/\_\]\[]/\./g;
        $fname =~ s/\.+/\./g;
      }
      $fname .= "\-$fdex";
      print "$keepi\t$responsei\t$f\t$fname\t$fdex\n";
      $fdex++;

      my %studies = ();
      # apply filtering criteria here
      foreach my $k1 (sort keys %meta){ # sample ID
        # is the keepi correct?
        if ($meta{$k1}{$keepi} ne "Yes"){
          next;
        }else{ #keep it
          push @{$studies{$meta{$k1}{"Study"}}}, $k1;
        }
      }

      my %nonZeroPrct = ();
      open OUT2, ">$outDir/$ft/raw-tables/$fname\.$keepi\.$responsei\.txt" or die;
      open OUT, ">$outDir/$ft/tables/$fname\.$keepi\.$responsei\.txt" or die;
      print OUT "study\tyear\tn.e\tmean.e\tsd.e\tn.c\tmean.c\tsd.c\n";
      foreach my $study (sort keys %studies){ # for this study
        my @samples        = @{$studies{$study}}; # get the valid samples
        my @responder      = ();
        my @nonresponder   = ();
        foreach my $s (@samples){
          if ($meta{$s}{$responsei} eq "R"){
            push @responder, $data{$f}{$s}; # add the feature value
          }elsif($meta{$s}{$responsei} eq "NR"){
            push @nonresponder, $data{$f}{$s};
          }
        }
        next if (!defined($responder[1]) or !defined($nonresponder[1]));

        my ($responderavg,$respondersamplestddev,$respondersamplestderr,$respondern)             = meanstde(@responder);
        my ($nonresponderavg,$nonrespondersamplestddev,$nonrespondersamplestderr,$nonrespondern) = meanstde(@nonresponder);

        my $nCnt = 0;
        my $gt0  = 0;
        foreach my $it (@responder){
          print OUT2 "$study\tR\t$it\n";
          if ($it > 0){
            $gt0++;
          }
          $nCnt++;
        }
        foreach my $it (@nonresponder){
          print OUT2 "$study\tNR\t$it\n";
          if ($it > 0){
            $gt0++;
          }
          $nCnt++;
        }
        $nonZeroPrct{$study} = 100*$gt0/$nCnt;

        print OUT "$study\t$YEAR{$study}\t$respondern\t$responderavg\t$respondersamplestddev\t$nonrespondern\t$nonresponderavg\t$nonrespondersamplestddev\n";
      } # end of study loop
      close OUT;
      close OUT2;

      # plotting
      open  RIN, ">$outDir/$ft/rcodes/$fname\.$keepi\.$responsei\.rcode.in" or die;
      print RIN
      "library(meta)
    t1<-read.table(\"$outDir/$ft/tables/$fname\.$keepi\.$responsei\.txt\",sep=\"\t\", header=TRUE)
    meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=t1, sm=\"SMD\", studlab=paste(study,year))
    summary(meta1)

      out1 <- paste(\"$outDir/$ft/pdfs\/\", \"$fname\.$keepi\.$responsei\.pdf\", sep=\"\")
      pdf(out1,paper=\"special\",width=16,height=5,bg=\"white\")
      forest(meta1, col.diamond.fixed=\"#2980B9\", col.diamond.random=\"#C0392B\", lab.e=\"Responder\", lab.c=\"Non-responder\", leftcols=c(\"studlab\", \"n.c\", \"n.e\"), leftlabs=c(\"Study\", \"Non-responder\", \"Responder\"), xlim=c(-2,2))
      dev.off()

    \n";
      close RIN;
      system("R CMD BATCH $outDir/$ft/rcodes/$fname\.$keepi\.$responsei\.rcode.in $outDir/$ft/rcodes/$fname\.$keepi\.$responsei\.rcode.out");

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

      open IN2, "$outDir/$ft/rcodes/$fname\.$keepi\.$responsei\.rcode.out";
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
          last;
        }
      }
      close IN2;

      print REPORT "$fname\t$keepi\t$responsei\t$fixedsmd\t$randsmd\t$fixedp\t$randp\t$I2\t$tausq";
      foreach my $study (@studies){
        print REPORT "\t$nonZeroPrct{$study}";
      }
      print REPORT "\n";

    } # end of features
  } # end of response
  } # end of keep
  close REPORT;
} # end of feature Types

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
