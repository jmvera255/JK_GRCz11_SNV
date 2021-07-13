use strict;
use Getopt::Std;

my(%opts);
getopts("i:k:s:", \%opts);

my $usage = "Usage: perl chrom_blocks.pl -i <genome.fa.fai> -k <int>

Where -k sets the # of kb per window.

optional: -s sets the step size of windows if a sliding window is desired.";

if(not defined($opts{k})){
  die $usage;
}

# convert kb to bp
$opts{k} = ($opts{k} * 1000) - 1;

# define slide/tiling length
my $slide = $opts{k} + 1;
if(defined($opts{s})){
   $slide = ($opts{s} * 1000);
}

open(FILE, "< $opts{i}") || die "Input file not found, stopped";
while(my $line =<FILE>){
  chomp($line);
  my @tabs = split("\t", $line);

  # print windows to stdout
  for(my $i = 1; $i < $tabs[1] - $opts{k}; $i = $i + $slide){
    my $stop = $i + $opts{k};
    print("$tabs[0]:$i-$stop\n");
  }
  # add equal length window at end of sequence/chrom
  my $start = $tabs[1] - $opts{k};
  print("$tabs[0]:$start-$tabs[1]\n");
}

