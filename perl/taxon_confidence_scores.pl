#!/usr/bin/env perl

=head1 NAME

  taxon_confidence_scores.pl

=head1 DESCRIPTION

  Provided a model tree and the clScores.txt file for that dataset, 
  produce the average +/- standard error for each taxon and 
  taxonomic level in the model tree. 

=head1 SYNOPSIS

  taxon_confidence_scores.pl -i <input dir> -o <output dir>

  **Currently, looks for model.tree and clScores.txt as input from the in directory.

=head1 OPTIONS

=over

=item B<--input-dir, -i>
  Input directory containing all.fa, spp.tx, and spp.lineage files. 

=item B<--output-dir, -o>
  Output directory.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  taxon_confidence_scores.pl -i ~/Desktop/PECAN/FFT_NS_2/FFT_NS_2/V3V4_models_bac_only -o ~/Desktop/PECAN/FFT_NS_2/FFT_NS_2

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::MoreUtils qw( part );
use List::Util qw( sum );
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(uniq);
use Cwd;
use Statistics::Basic qw(:all);


## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-dir|i=s"       => \my $inDir,
  "output-dir|o=s"      => \my $outDir,
 
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if ( !$inDir )
{
  warn "\n\n\tERROR: No input directory provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if ( $inDir && !$outDir )
{
  print "\n\n\tMissing out directory. Using input directory for output.";
  print "\n\n";
  $outDir = $inDir;
}

elsif ( !$inDir && !$outDir )
{
  warn "\n\n\tERROR: No input or output directories provided.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

####################################################################
##                               MAIN
####################################################################

my $modelTreeFile = "$inDir/model.tree";
print "---Using $modelTreeFile for model tree file\n";
my $clScoreFile = "$inDir/clScore_new.txt";
print "---Using $clScoreFile for classification score file\n";
my $sppSeqIDFile = "$inDir/sppSeqID.lineage";
print "---Using $sppSeqIDFile for classification score file\n";

my %clscore = readTbl($clScoreFile);
my %sppSeqID = read2colTbl($sppSeqIDFile);

print "---Obtaining all taxa (nodes & leaves) from $modelTreeFile\n";
my $modelTreeLabels = "$outDir/modelTree.labels";
my $cmd = "nw_labels $modelTreeFile > $modelTreeLabels";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

open IN, "<$modelTreeLabels" or die "Cannot open $modelTreeLabels for appending: $OS_ERROR";
my $taxonConf = "$outDir/taxon_confidence.txt";
print "---Writing taxon confidence values to $taxonConf...\n";
open OUT, ">$taxonConf" or die "Cannot open $taxonConf for appending: $taxonConf";
print OUT "taxon\tnSeqs\tavg±stderr\n";
## Calculation loop
my $spp;
foreach my $taxon ( <IN>)
  {
    chomp $taxon;
    my $sum;
    my $avg;
    my $stddev;
    my $stderr;
    if ($debug)
      {
        print "\nSearching for $taxon in lineages\n";
      }
    my @seqIDs = grep { $sppSeqID{$_} =~ $taxon } keys %sppSeqID;
    if ($debug)
      {
        print scalar @seqIDs . " sequences found for $taxon.\n";
      }
    my @clscores;
    foreach my $id ( @seqIDs ) 
    { 
      if ( exists $clscore{$id} )
      {
        my $score = $clscore{$id};
        push @clscores, $score;
      }
    }
    if (scalar @clscores > 0)
    {
      $sum = sum @clscores;
      $avg = $sum / scalar @clscores;
      $stddev = stddev \@clscores;
      $stderr = $stddev / scalar @clscores;
      print OUT "$taxon\t". scalar @clscores . "\t$avg±$stderr\n";
      if ($debug)
      {
        print "Score for $taxon is $avg ± $stderr\n";
      }
    }
    else
    {
      print "clScores not found for $taxon\n";
    }
   $spp++; 
  }
close IN;
close OUT;
print "---Confidence value calculations completed for $spp taxa.\n";

####################################################################
##                               SUBS
####################################################################
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub read2colTbl{

  my $file = shift;
  my %tbl;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

exit;