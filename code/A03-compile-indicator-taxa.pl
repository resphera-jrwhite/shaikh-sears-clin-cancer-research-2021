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
# ---------------------------------------------------------
# SET OUT DIR
# ---------------------------------------------------------
my $outDir = "../analysis/A03-compile-indicator-taxa";
if (-e $outDir){
  `rm -r $outDir`;
}
`mkdir $outDir`;
# ---------------------------------------------------------------
my %YEAR     = ("Chaput"         => 2017,
                "Gopalakrishnan" => 2017,
                "Matson"         => 2018,
                "Frankel"        => 2017,
                "Routy NSCLC"    => 2017,
                "Routy RCC"      => 2017);
# ---------------------------------------------------------------
# load merged file across all studies
my %features = ("Akkermansia_muciniphila"=>1,
"Alistipes_indistinctus"=>1,
"Anaerostipes_hadrus"=>1,
"Anaerotruncus_colihominis"=>1,
"Bacteroides_caccae"=>1,
"Bacteroides_coprocola"=>1,
"Bacteroides_fragilis"=>1,
"Bacteroides_thetaiotaomicron"=>1,
"Bacteroides_uniformis"=>1,
"Bifidobacterium_adolescentis"=>1,
"Bifidobacterium_longum"=>1,
"Clostridium_hathewayi"=>1,
"Clostridium_hylemonae"=>1,
"Clostridium_methylpentosum"=>1,
"Collinsella_aerofaciens"=>1,
"Dorea_formicigenerans"=>1,
"Enterococcus_faecium"=>1,
"Faecalibacterium_prausnitzii"=>1,
"g.Akkermansia"=>1,
"g.Bifidobacterium"=>1,
"g.Enterococcus"=>1,
"g.Eubacterium"=>1,
"g.Faecalibacterium"=>1,
"g.Lactobacillus"=>1,
"g.Parabacteroides"=>1,
"g.Ruminococcus"=>1,
"Gemmiger_formicilis"=>1,
"Holdemania_filiformis"=>1,
"Klebsiella_pneumoniae"=>1,
"Megasphaera_micronuciformis"=>1,
"Methanobrevibacter_smithii"=>1,
"Oribacterium_sinus"=>1,
"Oxalobacter_formigenes"=>1,
"Parabacteroides_merdae"=>1,
"Parasutterella_excrementihominis"=>1,
"Prevotella_buccalis"=>1,
"Roseburia_hominis"=>1,
"Roseburia_intestinalis"=>1,
"Ruminococcus_obeum"=>1,
"Scardovia_wiggsiae"=>1,
"Streptococcus_parasanguinis"=>1,
"Veillonella_parvula"=>1);

my %hasdata = ();
my @files = qw/out.A01-species.txt out.A01-genus.txt/;
for my $f (@files){
  my @header = ();
  open IN, "../analysis/A01-merge-profiles-across-studies/$f" or die;
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
      next;
    }
    for my $i (1 .. $#A){
      $data{$A[0]}{$header[$i]} = $A[$i];
      $hasdata{$header[$i]} = 1;
    }
  }
  close IN;
}


my %meta   = ();
my @header = ();
open IN, "../data/final/metadata.2019-04-22.txt" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($header[1])){
    @header = @A;
    next;
  }
  next if (!defined($hasdata{$A[0]}));
  for my $i (1 .. $#A){
    $meta{$A[0]}{$header[$i]} = $A[$i];
  }
}
close IN;

my @mcols = ("Study", "PatientID", "CancerType", "Timept", "DataType",
"Response", "Clean Seqs", "UseforOrig", "UseforMeta");

open OUT, ">$outDir/out.A03.indicator-taxa.txt" or die;
print OUT "SampleID";
for my $j (0 .. $#mcols){
  print OUT "\t$mcols[$j]";
}
foreach my $f (sort keys %features){
  print OUT "\t$f";
}
print OUT "\n";

foreach my $id (sort keys %meta){
  print OUT "$id";
  for my $j (0 .. $#mcols){
    print OUT "\t$meta{$id}{$mcols[$j]}";
  }
  foreach my $f (sort keys %features){
    print OUT "\t$data{$f}{$id}";
  }
  print OUT "\n";
}
close OUT;


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
