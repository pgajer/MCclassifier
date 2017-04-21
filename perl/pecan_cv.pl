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

  pecan_cv.pl --offset-coef 0.9 -n 10 -i Firmicutes_group_6_V3V4

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::MoreUtils qw( part );

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $offsetCoef = 0.99;
my $txSizeThld = 10;

GetOptions(
  "input-group|i=s"  => \my $grPrefix,
  "num-folds|n=i"    => \my $nFolds,
  "offset-coef|o=f"  => \$offsetCoef,
  "tx-size-thld|t=i" => \$txSizeThld,
  "verbose|v"       => \my $verbose,
  "debug"           => \my $debug,
  "dry-run"         => \my $dryRun,
  "help|h!"         => \my $help,
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
elsif ( !$nFolds )
{
  warn "\n\n\tERROR: Missing number of folds parameter";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

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

print "\r--- Extracting seq IDs from trimmed alignment fasta file              ";
my @seqIDs = get_seqIDs_from_fa($faFile);

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
my $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
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

##
## perform split of @seqIDs into nCV part
##

## permut @seqIDs
##my $nSeqs = @seqIDs;
##my @a = 1..$nSeqs;
fisher_yates_shuffle(\@seqIDs);

my $i = 0;
my @parts = part { $i++ % $nFolds } @seqIDs;

if ($debug)
{
  print "no parts: " . @parts . "\n";
  for my $p (@parts)
  {
    print "Part size: " . @{$p} . "\n";
  }
  print "\n";
  print "Seq's in part 1: @{$parts[0]}\n\n";
}

my $cvReportFile = "cv_report.txt";
open ROUT, ">$cvReportFile" or die "Cannot open $cvReportFile for writing: $OS_ERROR";
print ROUT "itr\tsize\tpTPcommSpp\tnTPcommSpp\tnIDsCommSpp\tpFPnovelSpp\tnFPnovelSpp\tnIDsNovelSpp\n";

print "\r                                                                                              ";

## cross-validation loop
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
  $cmd = "classify -d $mcDir -i $testFaFile -o $cvDir";
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

  my $combTxFile = "$cvDir/comb.tx";
  open OUT, ">$combTxFile" or die "Cannot open $combTxFile for writing: $OS_ERROR";
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

  my @commSpp = comm(\@uqTrainSpp, \@uqTestSpp);
  my %commSppTbl = map { $_ => 1 } @commSpp;

  ## seq IDs in the test set from the common species
  my @testIDsCommSpp = grep { exists $commSppTbl{$testTx{$_}} } keys %testTx;
  my $nTestIDsCommSpp = @testIDsCommSpp;

  ## number of @testIDsCommSpp with correct classification
  my $nCorrectClCommSpp = 0;
  for (@testIDsCommSpp)
  {
    if ( $testTx{$_} eq $clTx{$_} )
    {
      $nCorrectClCommSpp++;
    }
  }

  ## percentage of correctly classified seq's from common species
  my $pCorrectClCommSpp = 100.0 * $nCorrectClCommSpp / $nTestIDsCommSpp;

  print "\r[$i] Test size: $nTestIDs                                                                                                      \n";
  print "[$i] No. test seqIDs from spp present in test and train sets: $nTestIDsCommSpp\n";
  print "[$i] Percentage correctly classified test seq's from common species: " . sprintf("%.2f", $pCorrectClCommSpp) . "\n";

  print ROUT "$i\t$nTestIDs\t". sprintf("%.2f", $pCorrectClCommSpp) . "\t$nCorrectClCommSpp\t$nTestIDsCommSpp\t";

  ## Identifying species present in test but not training. Novel species from the point of view of the training dataset
  my @novelSpp = diff(\@uqTestSpp, \@uqTrainSpp);
  my %novelSppTbl = map { $_ => 1 } @novelSpp;

  ## seq IDs in the test set from the novel species
  my @novelSppTestIDs = grep { exists $novelSppTbl{$testTx{$_}} } keys %testTx;
  my $nNovelSppTestIDs = @novelSppTestIDs;

  ## number of misclassified novel species. A novel species is misclassified if
  ## it is classified to any species. We know that no species classification
  ## cannot be correct b/c these species are not present in the training dataset.
  my $nNovelSppClErrors = 0;
  for (@novelSppTestIDs)
  {
    if ( $clTx{$_} !~ /^sg_/ && $clTx{$_} !~ /^\w_/ )
    {
      $nNovelSppClErrors++;
    }
  }

  ## percentage of correctly classified seq's from common species
  my $pNovelSppClErrors = 100.0 * $nNovelSppClErrors / $nNovelSppTestIDs;

  print "[$i] No. misclassified novel species: $nNovelSppClErrors\n";
  print "[$i] Percentage misclassified novel species: " . sprintf("%.2f", $pNovelSppClErrors) . "\n";

  print ROUT sprintf("%.2f", $pNovelSppClErrors) . "\t$nNovelSppClErrors\t$nNovelSppTestIDs\n";

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
  print "[$i] Completed in $timeStr\n\n";
}
close ROUT;

print "\n\n\tReport writtent to $grDir/$cvReportFile\n\n";

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

exit 0;
