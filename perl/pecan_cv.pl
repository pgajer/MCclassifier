#!/usr/bin/env perl

=head1 NAME

  pecan_cv.pl

=head1 DESCRIPTION

  Perform cross-validation on PECAN classifier.

=head1 SYNOPSIS

  pecan_cv.pl -i <phylo group>

=head1 OPTIONS

=over

=item B<--input-group-name, -i>
  Prefix of input group. For example, if Firmicutes_group_6_dir is a directory
  Firmicutes_group_6 group, then the input group name is "Firmicutes_group_6".

=item B<--num-folds, -n>
  Number of folds. That is the number of parts into which the data is split.

=item B<--offset-coef, -o>
  offset = log10(offsetCoef) is the amount of log10 posterior probability below the
  min(pp) that the error threshold is set.  Default value of offset coef: 0.99.

=item B<--tx-size-thld, -t>
  Size of taxon (number of its ref seq's) below which error threshold is computed
  differently than for the taxons of size at or above the threshold.

=item B<--skip-err-thld>
  Apply --skip-err-thld to the classifier

=item B<--cv-sp-size-thld, -s>
  Cross-validation species size threshold. Only species of that size or higher
  will have their ref seq's partitioned. Default value: 100.

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

  pecan_cv.pl --skip-err-thld --offset-coef 0.9 -n 10 -i Firmicutes_group_6_V3V4

  pecan_cv.pl --cv-sp-size-thld 100 --debug -n 10 -i Firmicutes_group_5_V3V4

  pecan_cv.pl --skip-err-thld --cv-sp-size-thld 50 -i Firmicutes_group_6_V3V4

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::MoreUtils qw( part );
use List::Util qw( sum );

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $nFolds       = 10;
my $offsetCoef   = 0.99;
my $txSizeThld   = 10;
my $cvSpSizeThld = 100;

GetOptions(
  "input-group|i=s"     => \my $grPrefix,
  "num-folds|n=i"       => \$nFolds,
  "offset-coef|o=f"     => \$offsetCoef,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "cv-sp-size-thld|s=i" => \$cvSpSizeThld,
  "skip-err-thld"       => \my $skipErrThld,
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

if ( !$grPrefix )
{
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
# elsif ( !$nFolds )
# {
#   print "\n\n\t10-fold CV\n\n";
# }

####################################################################
##                               MAIN
####################################################################

my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

my $debugStr = "";
my $quietStr = "--quiet";
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

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "\n\n\tERROR: $grDir does not exist";
  print "\n\n";
  exit 1;
}

chdir $grDir;

my $faFile      = $grPrefix . "_final.fa";
my $lineageFile = $grPrefix . "_final.lineage";
my $treeFile	= $grPrefix . "_final.tree";
my $txFile      = $grPrefix . "_final.tx";

if ( ! -e $lineageFile )
{
  warn "\n\n\tERROR: $lineageFile does not exist";
  print "\n\n";
  exit 1;
}
elsif ( ! -e $faFile )
{
  warn "\n\n\tERROR: $faFile does not exist";
  print "\n\n";
  exit 1;
}
elsif ( ! -e $treeFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n";
  exit 1;
}
elsif ( ! -e $txFile )
{
  warn "\n\n\tERROR: $txFile does not exist";
  print "\n";
  exit 1;
}

if ( -l $treeFile )
{
  $treeFile = readlink($treeFile);
}

print "\n";

print "\r--- Creating cross-validation reports directory";
my $cvReportsDir = "cv_reports_dir";
my $cmd = "mkdir -p $cvReportsDir";
print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "\r--- Extracting seq IDs from trimmed alignment fasta file              ";
my @seqIDs = get_seqIDs_from_fa($faFile);

if ($debug)
{
  print "\nNo. of seq's in the phylo group: " . @seqIDs . "\n";
}

print "\r--- Parsing lineage table                                             ";
my %lineageTbl = readLineageTbl($lineageFile);


## testing if lineage and fa files has the same seq IDs
print "\r--- Checking if seqIDs of $faFile and $lineageFile are the same       ";
my @lSeqIDs = keys %lineageTbl;
my @commSeqIDs = comm(\@seqIDs, \@lSeqIDs);
if (@commSeqIDs != @seqIDs || @commSeqIDs != @lSeqIDs)
{
  warn "WARNING: seq IDs of the fasta file and lineage file do not match";
  print "Number of elements in the fasta file: " . @seqIDs . "\n";
  print "Number of elements in the lineage file: " . @lSeqIDs . "\n";
  print "Number of elements common to the fasta and lineage files: " . @commSeqIDs . "\n";
  exit 1;
}

print "\r--- Testing if treeFile is consistent with faFile                      ";
## extracting leaves' IDs
my $treeLeavesFile = "$grPrefix" . "_tree.leaves";
$cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

## looking at the difference between leaf IDs and newTx keys
my @treeLeaves = readArray($treeLeavesFile);

my @commTL = comm(\@treeLeaves, \@seqIDs);
if (@commTL != @treeLeaves || @commTL != @seqIDs)
{
  warn "\n\n\tERROR: seq IDs of $treeFile and $faFile do not match";
  print "\tNumber of seq IDs of $treeFile: " . @treeLeaves . "\n";
  print "\tNumber of seq IDs of $faFile:   " . @seqIDs . "\n";
  print "\tNumber of common seq IDs elements: " . @commTL . "\n\n";
  exit 1;
}


my $offsetStr = sprintf("o%d", int(100*$offsetCoef));
my $cvReportFile = "$cvReportsDir/accuracy_report_$nFolds" . "FCV_$offsetStr" . ".txt";
if ($skipErrThld)
{
  $cvReportFile = "$cvReportsDir/accuracy_report_$nFolds" . "FCV_minSpSize$cvSpSizeThld" . "_skip_error_thld.txt";
}

open ROUT, ">$cvReportFile" or die "Cannot open $cvReportFile for writing: $OS_ERROR";
print ROUT "phGr\titr\ttestSize\tnKnownSpp\tnNovelSpp\tpTPcommSpp\tnTPcommSpp\tnIDsCommSpp\tpFPnovelSpp\tnFPnovelSpp\tnIDsNovelSpp\n";

print "\r                                                                                                                     ";


##
## Performing stratified (per species) partition of seq IDs
##

my %tx = readTbl($txFile);
my %spTbl;
for (keys %tx)
{
  push @{$spTbl{$tx{$_}}}, $_;
}

my @base;  # This is array of seq IDs of species with less than
           # $cvSpSizeThld. They will always be in the training set.
my @parts; # Array with nFolds elements, each a ref to an array with seq IDs of
           # one of the nFolds parts
my $nBigSpp = 0; # number of species with at least $cvSpSizeThld elements
my %bigSpp;      # $bigSpp{$sp} = size($sp)
for my $sp (keys %spTbl)
{
  if ( @{$spTbl{$sp}} >= $cvSpSizeThld )
  {
    my @a = @{$spTbl{$sp}};
    $bigSpp{$sp} = @a;
    $nBigSpp++;
    fisher_yates_shuffle(\@a);

    my $i = 0;
    my @spParts = part { $i++ % $nFolds } @a;

    for my $i (0..($nFolds-1))
    {
      push @{$parts[$i]}, @{$spParts[$i]};
    }
  }
  else
  {
    push @base, @{$spTbl{$sp}};
  }
}

if (1)
{
  my $nSpp = keys %spTbl;
  print "\n\nNo. all spp:                    $nSpp\n";
  print     "No. spp with at least $cvSpSizeThld seq's: $nBigSpp\n\n";

  my @spp = sort {$bigSpp{$b} <=> $bigSpp{$a}} keys %bigSpp;
  printFormatedTbl(\%bigSpp, \@spp);
  my $sum = 0;
  for (keys %bigSpp)
  {
    $sum += $bigSpp{$_};
  }

  print "Total number of seq's:      " . ($sum + @base) . "\n";
  print "No. of seq's in big species: $sum\n";

  print "\nNo. parts: " . @parts . "\n";
  for my $p (@parts)
  {
    print "Part size: " . @{$p} . "\n";
  }
  print "\n";
  # print "Seq's in part 1: @{$parts[0]}\n\n";
  # print "size(base): " . @base . "\n";
  # exit 0;
}

##
## cross-validation loop
##
my @pCorrectClKnownSpp;
foreach my $i (0..($nFolds-1))
{
  print "\r[$i] Creating cross-validation directory";
  my $cvDir = "cvDir_$i";
  my $cmd = "rm -rf $cvDir; mkdir $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my @testIDs = @{$parts[$i]};
  my @trainIDs;
  foreach my $j (0..($nFolds-1))
  {
    if ( $j != $i )
    {
      push @trainIDs, @{$parts[$j]};
    }
  }

  push @trainIDs, @base;

  my $nTestIDs = @testIDs;
  my $nTrainIDs = @trainIDs;

  ##print "i: $i\ttest size: " . @testIDs . "\ttrain size: " . @trainIDs . "\n";

  my %trainLineageTbl;
  @trainLineageTbl{@trainIDs} = @lineageTbl{@trainIDs};

  my %spLineage;
  for my $id ( keys %trainLineageTbl )
  {
    my $lineage = $trainLineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $subGe = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    $subGe = "sg_$subGe";
    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $spLineage{$sp} = "$subGe\t$ge\t$fa\t$or\t$cl\t$ph\td_Bacteria";
  }

  print "\r[$i] Creating training set species lineage file                 ";
  my $trSpLineageFile = "$cvDir/train.spLineage";
  open OUT, ">$trSpLineageFile" or die "Cannot open $trSpLineageFile for writing: $OS_ERROR";
  for my $sp (keys %spLineage)
  {
    print OUT "$sp\t" . $spLineage{$sp} . "\n";
  }
  close OUT;

  print "\r[$i] Writing training set seqIDs to a file                        ";
  my $trIDsFile = "$cvDir/train.seqIDs";
  open OUT, ">$trIDsFile" or die "Cannot open $trIDsFile for writing: $OS_ERROR";
  for my $id (@trainIDs)
  {
    print OUT "$id\n";
  }
  close OUT;

  print "\r[$i] Writing test seqIDs to a file                                  ";
  my $testIDsFile = "$cvDir/test.seqIDs";
  open OUT, ">$testIDsFile" or die "Cannot open $testIDsFile for writing: $OS_ERROR";
  for my $id (@testIDs)
  {
    print OUT "$id\n";
  }
  close OUT;

  print "\r[$i] Creating training set fasta file                                 ";
  my $trFaFile = "$cvDir/train.fa";
  $cmd = "select_seqs.pl $quietStr -s $trIDsFile -i $faFile -o $trFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating test set fasta file                                      ";
  my $testFaFile = "$cvDir/test.fa";
  $cmd = "select_seqs.pl $quietStr -s $testIDsFile -i $faFile -o $testFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating training set taxon file                                  ";
  my $trTxFile = "$cvDir/train.tx";
  $cmd = "select_tx.pl $quietStr -s $trIDsFile -i $txFile -o $trTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating test taxon file                                           ";
  my $testTxFile = "$cvDir/test.tx";
  $cmd = "select_tx.pl $quietStr -s $testIDsFile -i $txFile -o $testTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building model tree and creating taxon's reference fasta files      ";
  my $mcDir = "$cvDir/mcDir";
  $cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $trSpLineageFile -i $trFaFile -t $trTxFile -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building MC models                                              ";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Estimating error thresholds                                      ";
  $cmd = "est_error_thlds --offset-coef $offsetCoef --tx-size-thld $txSizeThld -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Running classify on test fasta file                               ";
  $cmd = "classify $skipErrThldStr -d $mcDir -i $testFaFile -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Creating a table of true and classifier's taxonomy                 ";
  my %trainTx = readTbl($trTxFile);
  my %testTx = readTbl($testTxFile);
  my $clTxFile = "$cvDir/MC_order7_results.txt";
  my %clTx = readTbl($clTxFile);

  ## Checking if both tables have the same seqIDs as their keys
  my @testTxIDs = keys %testTx;

  if ( @testTxIDs != @testIDs )
  {
    warn "ERROR: size(testTxIDs) != size(testIDs)";
    print "size(testTxIDs): " . @testTxIDs . "\n";
    print "size(testIDs): " . @testIDs . "\n";
    exit 1;
  }

  my @clIDs   = keys %clTx;
  my @c = comm(\@testIDs, \@clIDs);
  if (@c != @testIDs || @c != @clIDs)
  {
    warn "WARNING: seq IDs of the testTx and clTx tables do not match";
    print "Number of elements in testTx: " . @testIDs . "\n";
    print "Number of elements in clTx: " . @clIDs . "\n";
    print "Number of elements common to both tables: " . @c . "\n\n";

    if ( @testIDs > @clIDs )
    {
      my @d = diff(\@testIDs, \@clIDs);
      print "\nElements in testTx that are not clTx:\n";
      for (@d)
      {
	print "\t$_\t" . $testTx{$_} . "\n";
      }
      print "\n\n";
    }

    if ( @clIDs > @testIDs )
    {
      my @d = diff(\@clIDs, \@testIDs);
      print "\nElements in clTx that are not in testTx:\n";
      for (@d)
      {
	print "\t$_\t" . $clTx{$_} . "\n";
      }
      print "\n\n";
    }
    exit 1;
  }

  my $combTxFile = "$cvReportsDir/comb_$nFolds" . "FCV_$offsetStr" . ".tx";
  if ($skipErrThld)
  {
    $combTxFile = "$cvReportsDir/comb_$nFolds" . "FCV_minSpSize$cvSpSizeThld" . "_skip_error_thld.tx";
  }

  if ($i==0)
  {
    open OUT, ">$combTxFile" or die "Cannot open $combTxFile for writing: $OS_ERROR";
  }
  else
  {
    open OUT, ">>$combTxFile" or die "Cannot open $combTxFile for appending: $OS_ERROR";
  }

  for (@c)
  {
    print OUT "$_\t" . $testTx{$_} . "\t" . $clTx{$_} . "\n";
  }
  close OUT;

  ## Identifying species present in the test and training sets
  my @trainSpp = values %trainTx;
  my @uqTrainSpp = unique(\@trainSpp);

  my @testSpp  = values %testTx;
  my @uqTestSpp = unique(\@testSpp);

  my @knownSpp = comm(\@uqTrainSpp, \@uqTestSpp);
  my $nKnownSpp = @knownSpp;
  my %knownSppTbl = map { $_ => 1 } @knownSpp;

  ## seq IDs in the test set from the known species
  my @testIDsKnownSpp = grep { exists $knownSppTbl{$testTx{$_}} } keys %testTx;
  my $nTestIDsKnownSpp = @testIDsKnownSpp;

  ## number of @testIDsKnownSpp with correct classification
  my $nCorrectClKnownSpp = 0;
  for (@testIDsKnownSpp)
  {
    if ( $testTx{$_} eq $clTx{$_} )
    {
      $nCorrectClKnownSpp++;
    }
  }

  ## percentage of correctly classified seq's from known species
  ##my $pCorrectClKnownSpp = 100.0 * $nCorrectClKnownSpp / $nTestIDsKnownSpp;
  $pCorrectClKnownSpp[$i] = 100.0 * $nCorrectClKnownSpp / $nTestIDs;

  ## Identifying species present in test but not training. Novel species from the point of view of the training dataset
  my @novelSpp = diff(\@uqTestSpp, \@uqTrainSpp);
  my $nNovelSpp = @novelSpp;
  my $pNovelSppClErrors = "NA";
  my $nNovelSppClErrors = "NA";
  my $nNovelSppTestIDs  = "NA";
  if (@novelSpp > 0)
  {
    my %novelSppTbl = map { $_ => 1 } @novelSpp;

    ## seq IDs in the test set from the novel species
    my @novelSppTestIDs = grep { exists $novelSppTbl{$testTx{$_}} } keys %testTx;
    $nNovelSppTestIDs = @novelSppTestIDs;

    ## number of misclassified novel species. A novel species is misclassified if
    ## it is classified to any species. We know that no species classification
    ## cannot be correct b/c these species are not present in the training dataset.
    $nNovelSppClErrors = 0;
    for (@novelSppTestIDs)
    {
      if ( $clTx{$_} !~ /^sg_/ && $clTx{$_} !~ /^\w_/ )
      {
	$nNovelSppClErrors++;
      }
    }

    ## percentage of correctly classified seq's of known species
    ##$pNovelSppClErrors = 100.0 * $nNovelSppClErrors / $nNovelSppTestIDs;
    $pNovelSppClErrors = 100.0 * $nNovelSppClErrors / $nTestIDs;
  }

  ## $cmd = "cmp_tx.pl -i $testTxFile -j $clTxFile -o $cvReportFile";

  print "\r[$i] No. known spp:                                       $nKnownSpp                                       \n";
  print "[$i] No. novel spp:                                       $nNovelSpp\n";
  print "[$i] No. test seq's:                                      $nTestIDs       \n";
  print "[$i] No. seq's from known species:                        $nTestIDsKnownSpp\n";
  print "[$i] No. correctly classified seq's from known species:   $nCorrectClKnownSpp\n";
  print "[$i] Perc. correctly classified seq's from known species: " . sprintf("%.2f%%", $pCorrectClKnownSpp[$i]) . "\n";

  print ROUT "$grPrefix\t$i\t$nTestIDs\t$nKnownSpp\t$nNovelSpp\t". sprintf("%.2f", $pCorrectClKnownSpp[$i]) . "\t$nCorrectClKnownSpp\t$nTestIDsKnownSpp\t";

  if (@novelSpp > 0)
  {
    print "[$i] No. seq's from novel species:                        $nNovelSppTestIDs\n";
    print "[$i] No. misclassified novel species:                     $nNovelSppClErrors\n";
    print "[$i] Perc. misclassified novel species:                   " . sprintf("%.2f%%", $pNovelSppClErrors) . "\n";
    print ROUT sprintf("%.2f", $pNovelSppClErrors) . "\t$nNovelSppClErrors\t$nNovelSppTestIDs\n";
  }
  else
  {
    print ROUT "NA\tNA\t0\n";
  }
  print "\n";

}
close ROUT;


my $meanPercCorrectClKnownSpp = sum(@pCorrectClKnownSpp)/@pCorrectClKnownSpp;
print "\nMean Perc. correctly classified seq's from known species: " . sprintf("%.2f%%", $meanPercCorrectClKnownSpp) . "\n";

$endRun = time();
$runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  $timeMin = int($runTime / 60);
  $timeSec = sprintf("%02d", $runTime % 60);
  $timeStr = "$timeMin:$timeSec";
}
else
{
  $runTime = sprintf("%02d", $runTime);
  $timeStr = "$timeMin:$runTime";
}

$startRun = time();
print "\n\tCompleted in $timeStr\n\n";
print "\tReport writtent to $grDir/$cvReportFile\n\n";

####################################################################
##                               SUBS
####################################################################

# read lineage table
sub readLineageTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readLineageTbl(): $file does not exist";
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
    ## test for '/' characters
    if ($t =~ /\//)
    {
      warn "\n\n\tERROR: Discovered '/' for id: $id\t$t";
      print "\n\n";
      exit 1;
    }
  }
  close IN;

  return %tbl;
}

sub get_seqIDs_from_fa
{
  my $file = shift;

  my $quiet = 1;
  my $startRun = time();
  my $endRun = time();

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <IN>;
  my $count = 1;
  my $timeStr = "";
  my @seqIDs;
  while (<IN>)
  {
    if ( !$quiet && ($count % 500 == 0) )
    {
      $endRun = time();
      my $runTime = $endRun - $startRun;
      if ( $runTime > 60 )
      {
	my $timeMin = int($runTime / 60);
	my $timeSec = sprintf("%02d", $runTime % 60);
	$timeStr = "$timeMin:$timeSec";
      }
      else
      {
	my $runTime = sprintf("%02d", $runTime);
	$timeStr = "0:$runTime";
      }
      print "\r$timeStr";
    }

    chomp;
    my ($id,@seqlines) = split /\n/, $_;
    push @seqIDs, $id;
    $count++;
  }
  close IN;
  $/ = "\n";

  return @seqIDs;
}

# common part of two arrays
sub comm{

  my ($a1, $a2) = @_;

  my @c; # common array
  my %count;

  foreach my $e (@{$a1}, @{$a2}){ $count{$e}++ }

  foreach my $e (keys %count)
  {
    push @c, $e if $count{$e} == 2;
  }

  return @c;
}

# read table with one column
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readArray(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( defined $hasHeader )
  {
    <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    push @rows, $_;
  }
  close IN;

  return @rows;
}

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

# read two column table; create a table that assigns
# elements of the first column to the second column
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

# difference of two arrays
sub diff{

  my ($a1, $a2) = @_;

  my (%aa1, %aa2);

  foreach my $e (@{$a1}){ $aa1{$e} = 1; }
  foreach my $e (@{$a2}){ $aa2{$e} = 1; }

  my @d; # dfference array

  foreach my $e (keys %aa1, keys %aa2)
  {
    push @d, $e if exists $aa1{$e} && !exists $aa2{$e};
  }

  return @d;
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTbl{

  my ($rTbl, $rSub) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  for (@args)
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "$_$pad" . $rTbl->{$_} . "\n";
  }
  print "\n";
}

exit 0;
