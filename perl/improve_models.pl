#!/usr/bin/env perl

=head1 NAME

  improve_models.pl

=head1 DESCRIPTION

  Using the classification scores produced by score_classification.pl, determine which PECAN 
  models can be improved by removing sequences with scores of < 1. The script proceeds as follows: 

    1. Identify models with poorly-classified reference sequences (sub poor_classification)
        - All correct (sum of scores per taxon/nSeqs of taxon = 1)
        - Partial (0 < sum of scores per taxon/nSeqs of taxon < 1)
        - All incorrect (sum of scores per taxon/nSeqs of taxon = 0)
    2. Remove misclassified reference sequences from models if they have a score < 1
    3. Produce a whole new set of models from new .tx and .lineage files with above seq's removed
    3. Test new model classification abilities using score_classification subroutine.
    4. Test whether the resulting avg classification score = 100% or not. If not, then returns to
       step 1. using newly produced models, classification scores, lineages, taxonomy files, etc.

=head1 SYNOPSIS

  improve_models.pl -i <input dir> -o <output dir>

=head1 OPTIONS

=over

=item B<--input-dir, -i>
  Input directory containing all.fa, spp.tx, and spp.lineage, clScore.txt, and MC_order7_results.txt files. 

=item B<--output-dir, -o>
  Output directory.

=item B<--lineage-file, -l>
  If a lineage file other than spp.lineage is to be used, specify here.

=item B<--taxonomy-file, -t>
  If a taxonomy file other than spp.tx is to be used, specify here.

=item B<--classification-score-file, -c>
  If a score file other than clScore.txt is to be used, specify here.

=item B<--pecan-output-file, -p>
  If a pecan classification file other than MC_order7_results.txt is to be used, specify here.

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

  improve_models.pl -i ~/Desktop/PECAN/FFT_NS_2/FFT_NS_2/V3V4_models_bac_only -o ~/Desktop/PECAN/FFT_NS_2/FFT_NS_2

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

## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################
GetOptions(
  "input-dir|i=s"       => \my $inDir,
  "output-dir|o=s"      => \my $outDir1,
  "lineage-file|l=s"    => \my $LineageFile,
  "taxonomy-file|t=s"   => \my $TxFile,
  "class-scores|c=s"    => \my $clScoreFile,
  "pecan-output|p=s"    => \my $pecanFile,
  "fasta-file|f=s"      => \my $faFile,
  "offset-coef|z=s"     => \my $offsetCoef,
  "tx-size-thld|n=s"    => \my $txSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
  "skip-err-thld"       => \my $skipErrThld,
  "quiet"               => \my $quiet,
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

if ( $inDir && !$LineageFile && !$TxFile && !$clScoreFile && !$pecanFile)
{
  $LineageFile = "$inDir/spp.lineage";
  $TxFile = "$inDir/spp.tx";
  $clScoreFile = "$inDir/clScore.txt";
  $pecanFile = "$inDir/MC_order7_results.txt";
}
elsif ( !$LineageFile && !$TxFile && !$clScoreFile && !$pecanFile && !$inDir && !$outDir1 )
{
  warn "\n\n\tERROR: No input, output directories or files provided.";
  exit 0;
}
elsif ( !$inDir )
{
  warn "\n\n\tERROR: No input directory provided.";
  $outDir1 = $inDir;
}
elsif ( !$outDir1 )
{
  warn "\n\n\tERROR: No input directory provided.";
  $outDir1 = "merge_model_out";
}

if ( $inDir && !$LineageFile )
{
  $LineageFile = "$inDir/spp.lineage";
}

if ( $inDir && !$TxFile )
{
  $TxFile = "$inDir/spp.tx";
}

if ( $inDir && !$clScoreFile )
{
  $clScoreFile = "$inDir/clScore.txt";
}

if ( $inDir && !$pecanFile)
{
  $pecanFile = "$inDir/MC_order7_results.txt";
}

if ( $inDir && !$faFile)
{
  $faFile = "$inDir/all.fa";
}

my $quietStr = "";
if ( $quiet  )
{
  $quietStr = "--quiet";
}

my $debugStr = "";
if ($debug)
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ($verbose)
{
  $verboseStr = "--verbose";
}

my $skipErrThldStr = "";
if ($skipErrThld)
{
  $skipErrThldStr = "--skip-err-thld";
}

my $ppEmbeddingStr = "";
if ($ppEmbedding)
{
  $ppEmbeddingStr = "--pp-embedding"
}

if (!$offsetCoef)
{
  $offsetCoef = 0.9;
}

if (!$txSizeThld)
{
  $txSizeThld = 10;
}
####################################################################
##                               MAIN
####################################################################

my $run = 0;

START:

my $outDir = $outDir1 . "_" . $run;
print "\n---Using $LineageFile for lineage file\n";
print "---Using $TxFile for taxonomy file\n";
print "---Using $clScoreFile for classification score file\n";
print "---Using $outDir for writing new files\n";

print "---Printing results to $outDir\n";
my $cmd = "rm -rf $outDir; mkdir $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

## Read in scores, taxonomy, and PECAN output file
my %clscore = readTbl($clScoreFile);
my %tx = readTbl($TxFile);
my %pecan = readTbl($pecanFile); 

## Make a unique array of taxa
my @uniqtx = uniq (values %tx);

print "---Determining models with poorly-classified sequences\n";
my ($rmisclassified) = poor_classification(\@uniqtx);

## Push to an single array all of the seq's that were misclassified (score < 1)
## For removal from dataset.

my $rmFile = "$outDir/seqs_to_remove.seqID";
print "---Writing IDs of misclassified sequences to $rmFile\n";
open OUT, ">$rmFile" or die "Cannot open $rmFile for writing: $OS_ERROR";
print OUT @$rmisclassified;
close OUT;

print "---Removing misclassified sequences (classification score < 1) from taxonomy file\n";
my $newtxFile = "$outDir/spp_new.tx";
$cmd = "select_tx.pl -i $TxFile -e $rmFile -o $newtxFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Removing misclassified sequences (classification score < 1) from lineage file\n";
my $newLinFile = "$outDir/spp_new.lineage";
$cmd = "select_fullTx.pl -t $newtxFile -f $LineageFile -o $newLinFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Building new fasta file with misclassified sequences removed\n";
my $reducefaFile = $outDir . "/select.fa";
$cmd = "select_seqs.pl -i $faFile -s $newtxFile -o $reducefaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Building model tree and creating taxon's reference fasta files\n";
$cmd = "buildModelTree $quietStr -l $newLinFile -i $reducefaFile -t $newtxFile -o $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Building MC models\n";
$cmd = "buildMC -t $outDir/spp_paths.txt -k 8 -d $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Estimating error thresholds\n";
$cmd = "est_error_thlds --offset-coef $offsetCoef --tx-size-thld $txSizeThld -d $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Classifying $reducefaFile using new models\n";
$cmd = "classify $skipErrThldStr -d $outDir -i $reducefaFile -o $outDir"; # $ppEmbeddingStr
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Comparing ref seq's taxonomy with the classification results\n";
my $sppSeqID = "$outDir/sppSeqID.lineage";
$cmd = "get_lineage.pl -a $newtxFile -l $newLinFile -o $sppSeqID";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "---Using $sppSeqID for lineage file\n";
my %newlineage = read2colTbl($sppSeqID);

my $newpecanFile = "$outDir/MC_order7_results.txt";
my %newpecan = readTbl($newpecanFile); 
print "---Using $newpecanFile for lineage file\n";

my @clIDs = keys %newpecan;

print "---Removed ". scalar @$rmisclassified . " misclassified sequences from taxonomy & lineage files.\n";
print "---Models built from $newLinFile, $newtxFile, and $reducefaFile.\n";
print "---Scoring ref seq's taxonomy with the new classification results\n";
my $rclScore = score_classification( \@clIDs, \%newlineage, \%newpecan);
my %clScore = %$rclScore;
my $clScore_sum = sum values %clScore;
my $clScore_perc = 100.0 * $clScore_sum /  scalar @clIDs; 

if ($skipErrThldStr)
{
  print "---Without error thresholds, the new classification score for ". scalar @clIDs . " is: " . $clScore_perc . "\n";
}
else
{
  print "---Using error thresholds, the new classification score for ". scalar @clIDs . " is: " . $clScore_perc . "\n";
}

my $newclScoreFile = "$outDir/clScore_new.txt";
print "---New classification scores written to $newclScoreFile\n";
open OUT, ">$newclScoreFile" or die "Cannot open $newclScoreFile for appending: $OS_ERROR";
print OUT "seqID\tclassification_score\n";
foreach my $key (keys %clScore)
{
  print OUT $key . "\t". $clScore{$key} . "\n";
}
close OUT;

if ($clScore_perc < 100)
{
  $LineageFile = $newLinFile;
  $TxFile = $newtxFile;
  $clScoreFile = $newclScoreFile;
  $inDir = $outDir;
  $run++;
  print "Classification score average = $clScore_perc. Model improvement continuing to iteration $run.\n";
  goto START;
}
else 
{
  print "Classification score average = $clScore_perc. Model improvement complete.\n";
}


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

sub poor_classification{
  my $runiqtx = shift;
  my @uniqtx = @$runiqtx;
  my $spp=0;
  my @seqIDs;
  my %sppSeqID;
  my %spsums;
  my $score;
  my $sum;
  my $avg;
  my $id;
  my @seqs_to_remove;
  my @seqs_greater_than_50;
  my @seqs_equal_to_50;
  my @misclassified;
  my @poorclass;


  ## For each unique taxon of the reference data 
  ## Obtain all seqIDs of that taxon...
  ## and build array of corresponding classification scores
  foreach my $key ( @uniqtx) 
  {
    @seqIDs = grep { $tx{$_} eq $key } keys %tx; 
    my @clscores;
    foreach $id ( @seqIDs ) 
    {
      if ( exists $clscore{$id} ) 
      {
        $score = $clscore{$id}; 
        push @clscores, $score; 
        if ($debug)
        {
          print "The score for $id was $clscore{$id}\n";
        }
      }
    }

    ##Place the clScore sum into a hash for that species.
    $sum = sum @clscores;
    $avg = $sum / scalar @clscores;
    $spsums{$key} = $avg;
    $spp++;

    ## If the mean classification score for a taxon is eq 1, 
    ## Go to the next taxon
    if ( $avg == 1 )
    {
      if ($debug)
      {
        print "All sequences of $key are correctly classified\n";
      }
    }
    ## If the average classification score of a taxon is less than 1
    ## loop through the sequences of that taxon and if their score
    ## is 0, push the seqID to the misclassified array
    elsif ( $avg < 1 )
    {
      foreach $id ( @seqIDs ) 
      {
        $score = $clscore{$id};
        if ($score < 1)
          {
          push @misclassified, $id."\n";
          }
      }
    } 
  }
return (\@misclassified)
}

sub score_classification{ 
  my $rclIDs = shift;
  my $rtestLineage = shift;
  my $rclTx = shift;

  my %testLineage = %$rtestLineage;
  my %clTx = %$rclTx;

  my %clScore;

  my @clIDs = @$rclIDs;
  foreach my $s (@clIDs) ## For each classified sequence
  {
    chomp $s;
    ##print "\nHere is the lineage for $s: ". $testLineage{$s} . "\n";

    ## t is equal to the classified taxa
    my $t = $clTx{$s};
    if (!$t)
    {
      next;
    } 
    else
    {
      ##print "The value of $s is $clTx{$s}.\n";
    }
    ## a is equal to the reference lineage
    my $a = $testLineage{$s};
    if (!$a)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in testLineage hash table\n;"
      }
      next;
    }
    #print "\nThe value of $s is $testLineage{$s}\n";
    ## if the classified taxa is found in the reference lineage
    my $value = 0;
    if ( grep /$t/, $a) 
    {   
      #print "\n $t is part of $a.\n";
      $value++;
      if ($debug)
      {
        print "\n $t is part of $a.\n";
      }
      if ( $value > 0) ##If the classifed taxonomy is found in the array
      {
        if ($clTx{$s} =~ /^d_/) 
        {
          $clScore{$s} = 0;
        }
        elsif ($clTx{$s} =~ /^p_/) 
        {
          $clScore{$s} = 0.25;
        }
        elsif ($clTx{$s} =~ /^c_/) 
        {
          $clScore{$s} = 0.325;
        }
        elsif ($clTx{$s} =~ /^o_/) 
        {
          $clScore{$s} = 0.5;
        }
        elsif ($clTx{$s} =~ /^f_/) 
        {
          $clScore{$s} = 0.625;
        }
        elsif ($clTx{$s} =~ /^g_/) 
        {
          $clScore{$s} = 0.75;
        }          
        elsif ($clTx{$s} =~ /^sg_/) 
        {
          $clScore{$s} = 0.825;
        }
        else 
        {
          $clScore{$s} = 1;
        }
      }
      elsif ($value = 0) ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
      {
        $clScore{$s} = 0;
      }   
      #print "\nThe score for $s is $clScore{$s}\n";
    }
    else ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
    {
      $clScore{$s} = 0;
    } 
  }
return \%clScore;
}

exit;