#!/usr/bin/env perl

=head1 NAME

  pecan_cv3.pl

=head1 DESCRIPTION

  Perform cross-validation on all bacteria models in given directory. This process considers 
  the entire lineage of the correct taxonomy as follows:
    - Each species has 7 levels of taxonomy, so correct classification to each of those levels
      equals a point of 1.00/8 = 0.125 x the level reached. 
    - Thus correct classification to each level receives a score as follows: 
        Domain (d_):    1 x 0 = 0
        Phylum (p_):    2 x 0.125 = 0.25
        Class (c_):     3 x 0.125 = 0.325
        Order (o_):     4 x 0.125 = 0.5
        Family (f_):    5 x 0.125 = 0.625
        Genus (g_):     6 x 0.125 = 0.75
        Subgenus (sg_): 7 x 0.125 = 0.875
        Species (s_):   8 x 0.125 = 1

        Where 1 is the best score, and 0.125 is the worst score. 

=head1 SYNOPSIS

  pecan_cv2.pl -o <output dir>

=head1 OPTIONS

=over

=item B<--input-dir, -i>
  Input directory containing all.fa, spp.tx, and spp.lineage files. 

=item B<--output-dir, -o>
  Output directory.

=item B<--num-folds, -n>
  Number of folds. That is the number of parts into which the data is split.

=item B<--offset-coef, -c>
  offset = log10(offsetCoef) is the amount of log10 posterior probability below the
  min(pp) that the error threshold is set.  Default value of offset coef: 0.9.

=item B<--tx-size-thld, -t>
  Size of taxon (number of its ref seq's) below which error threshold is computed
  differently than for the taxons of size at or above the threshold.

=item B<--skip-err-thld>
  Apply --skip-err-thld to the classifier

=item B<--cv-sp-size-thld, -s>
  Cross-validation species size threshold. Only species of that size or higher
  will have their ref seq's partitioned. Default value: 100.

=item B<--pp-embedding>
  Run classifier with the --pp-embedding flag.

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

  pecan_cv2.pl --cv-sp-size-thld 50 --offset-coef 0.9 -n 10 -o all_bacteria_CV_April_25_run1

  pecan_cv2.pl --skip-err-thld --cv-sp-size-thld 50 -o all_bacteria_CV_April_25_run2

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::MoreUtils qw( part );
use List::Util qw( sum );
use Data::Dumper qw(Dumper);

## use Math::Permute::Array;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $nFolds       = 10;
my $offsetCoef   = 0.9;
my $txSizeThld   = 10;
my $cvSpSizeThld = 50;

GetOptions(
  "input-dir|i=s"       => \my $mcDir,
  "output-dir|o=s"      => \my $outDir,
  "num-folds|n=i"       => \$nFolds,
  "offset-coef|c=f"     => \$offsetCoef,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "cv-sp-size-thld|s=i" => \$cvSpSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
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

if ( !$outDir )
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

my $ppEmbeddingStr = "";
if ($ppEmbedding)
{
  $ppEmbeddingStr = "--pp-embedding"
}



if ( ! -d $mcDir )
{
  warn "\n\n\tERROR: $mcDir does not exist";
  print "\n\n";
  my $mcDir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/all_bacteria_V3V4_MC_models_dir";
}

my $faFile    = $mcDir . "/select.fa";
my $txFile    = $mcDir . "/spp_new.tx";
my $spLgFile  = $mcDir . "/spp_new.lineage";


if ( ! -e $spLgFile )
{
  warn "\n\n\tERROR: $spLgFile does not exist";
  print "\n\n";
  exit 1;
}

elsif ( ! -e $faFile )
{
  warn "\n\n\tERROR: $faFile does not exist";
  print "\n\n";
  exit 1;
}

elsif ( ! -e $txFile )
{
  warn "\n\n\tERROR: $txFile does not exist";
  print "\n";
  exit 1;
}

# if ( -l $treeFile )
# {
#   $treeFile = readlink($treeFile);
# }

print "\n";

print "\r--- Creating cross-validation reports directory";
my $cvReportsDir = "$outDir/cv_reports_dir";
my $cmd = "rm -rf $outDir; mkdir -p $cvReportsDir";
print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


my $offsetStr = sprintf("o%d", int(100*$offsetCoef));
my $cvReportFile = "$cvReportsDir/accuracy_report_$nFolds" . "FCV_$offsetStr" . ".txt";
if ($skipErrThld)
{
  $cvReportFile = "$cvReportsDir/accuracy_report_$nFolds" . "FCV_minSpSize$cvSpSizeThld" . "_skip_error_thld.txt";
}



##
## Performing stratified (per species) partition of seq IDs
##

my %tx = readTbl($txFile);
## Create hash table of species => seqID
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

## For each species, consider the number of sequences representing it, 
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

if ($nBigSpp==0)
{
  my $nSpp = keys %spTbl;
  print "\n\nNo. all spp:                    $nSpp\n";
  print     "No. spp with at least $cvSpSizeThld seq's: $nBigSpp\n\n";
  print "\n\tNothing to be done. Exiting !!!\n\n";
  exit 0;
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
    ##print "Part size: " . (@{$p} + @base) . "\n";
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
my %clScore;
my @a;

foreach my $i (0..($nFolds-1))
{
  print "\r[$i] Creating cross-validation directory";
  ##my $cvDir = "$outDir/cvDir_$i";
  my $cvDir = "$outDir/cvDir";
  my $cmd = "rm -rf $cvDir; mkdir $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ## Create training set of seqIDs

  my @testIDs = @{$parts[$i]};  
  my @trainIDs;

  foreach my $j (0..($nFolds-1))
  {
    if ( $j != $i )
    {
      push @trainIDs, @{$parts[$j]};
    }
  }

  ##push @testIDs, @base; # this would obscure misclassification error on test sequences as @base would be correctly classified
  push @trainIDs, @base;

  my $nTestIDs = @testIDs;
  my $nTrainIDs = @trainIDs;

  ##print "i: $i\ttest size: " . @testIDs . "\ttrain size: " . @trainIDs . "\n";


  ## Print training seqIDs to seqID file 
  print "\r[$i] Writing training set seqIDs to a file                        ";
  my $trIDsFile = "$cvDir/train.seqIDs";
  print_seqID_to_file($trIDsFile, \@trainIDs);

  ## Print testing seqIDs to seqID file 
  print "\r[$i] Writing test seqIDs to a file                                  ";
  my $testIDsFile = "$cvDir/test.seqIDs";
  print_seqID_to_file($testIDsFile, \@testIDs);

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

  ## The training lineage file has to be made otherwise empty .fa files are made
  ## which will cause est_err_thld to fail
  print "\r[$i] Creating training set lineage file                                  ";
  my $trSpLineageFile = "$cvDir/train.lineage";
  $cmd = "select_fullTx.pl -t $trTxFile -f $spLgFile -o $trSpLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  
  my $mcDir = "$cvDir/mcDir";

  print "\r[$i] Creating seqID to lineage file                                  ";
  my $testSeqIDLineageFile = "$cvDir/sppSeqID.lineage";
  $cmd = "get_lineage.pl -a $testTxFile -l $spLgFile -o $testSeqIDLineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  ## the only difference is that models will be using fewer sequences for modeling big species

  #my $trSpLineageFile = $spLgFile; 

  print "\r[$i] Building model tree and creating taxon's reference fasta files      ";
  $cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $trSpLineageFile -i $trFaFile -t $trTxFile -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Building MC models                                              ";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ($ppEmbedding)
  {
    print "\r[$i] Generating log pp tables for internal nodes of the model tree                                              ";
    $cmd = "pp_embedding -d $mcDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    ##system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  print "\r[$i] Estimating error thresholds                                      ";
  $cmd = "est_error_thlds -d $mcDir --offset-coef $offsetCoef --tx-size-thld $txSizeThld";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\r[$i] Running classify on test fasta file                               ";
  $cmd = "time classify $ppEmbeddingStr $skipErrThldStr -d $mcDir -i $testFaFile -o $cvDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ($ppEmbedding)
  {
    $cmd = "cat $mcDir/classification_paths.txt > $outDir/classification_paths.txt";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  print "\r[$i] Creating a table of true and classifier's taxonomy                 ";
  my %trainTx = readTbl($trTxFile);
  my %testTx = readTbl($testTxFile);
  my $clTxFile = "$cvDir/MC_order7_results.txt";
  my %clTx = readTbl($clTxFile);
  my %testLineage = read2colTbl($testSeqIDLineageFile);

  if ($debug)
  {
    print Dumper \%testLineage;
  }

  ## Checking if both tables have the same seqIDs as their keys
  my @testTxIDs = keys %testTx;

  if ( @testTxIDs != @testIDs )
  {
    warn "ERROR: size(testTxIDs) != size(testIDs)";
    print "size(testTxIDs): " . @testTxIDs . "\n";
    print "size(testIDs): " . @testIDs . "\n";
    exit 1;
  }

  
  ## Checking if both the classified and test table have the same seqIDs as their keys
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

  ## Perform scoring system.
  my $clIDs;
  my %clScore = score_classification( \@clIDs, \%testLineage, \%clTx);
  my $clScore_sum = sum values %clScore;
  my $clScore_perc = 100.0 * $clScore_sum /  scalar keys %clScore;
  if ($skipErrThld) 
  {
    print "\n\nWithout error thresholds, the classification score is: " . $clScore_perc . "\n";
  }
  else 
  {
    print "\n\nUsing error thresholds, the classification score is: " . $clScore_perc . "\n";
  }
  my $clScoreFile = "$cvReportsDir/clScore_[$i].txt";
  my $clCompFile= "$cvReportsDir/cl_cmp_[$i].txt";
  open OUT, ">$clScoreFile" or die "Cannot open $clScoreFile for appending: $OS_ERROR";
  open OUT1, ">$clCompFile" or die "Cannot open $clCompFile for appending: $OS_ERROR";
  print OUT1 "seqID\tclassified_lineage\tcorrect_lineage_sg\tcorrect_lineage_g\tcorrect_lineage_f\tcorrect_lineage_o\tcorrect_lineage_c\tcorrect_lineage_p\tcorrect_lineage_d\tcorrect_lineage_sp\n";
  foreach my $key (keys %clScore)
  {
    print OUT1 $key . "\t". $clTx{$key} ."\t". $testLineage{$key} . "\n";
    print OUT $key . "\t" . $clScore{$key} . "\n";

  }
  close OUT;
  close OUT1;
}

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
print "\tReport written to $cvReportFile\n\n";

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

sub print_seqID_to_file
{
  my $IDsFile = shift;
  my $rIDs = shift;

  my @IDs = @{ $rIDs};

  open OUT, ">$IDsFile" or die "Cannot open $IDsFile for writing: $OS_ERROR";
  for my $id (@IDs)
  {
    print OUT "$id\n";
  }
  close OUT;
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
      if ($debug)
      {
        print "\nWARNING: $s not found in $rclTx\n";
      }
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
return %clScore;
}

exit 0;