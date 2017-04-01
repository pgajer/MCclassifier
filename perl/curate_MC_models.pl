#!/usr/bin/env perl

=head1 NAME

  curate_MC_models.pl

=head1 DESCRIPTION

  Given the prefix of a phylum group, the script computes random sequenes of all
  MC models, then selects only species sequences, runs classifier on them and
  then cmp_tx.pl.

=head1 SYNOPSIS

  curate_MC_models.pl -i <prefix> [Options]

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Groups prefix.

=item B<--num-rand-seqs, -n>
  Number of random sequences to be generated from each model.

=item B<--quiet>
  Do not print progress messages.

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

  curate_MC_models.pl -i Firmicutes_group_6_V3V4 -n 100 -l 429

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "group-prefix|i=s"  => \my $grPrefix,
  "num-rand-seqs|n=i" => \my $nRand,
  "rand-seqs-len|l=i" => \my $randSeqLen,
  "igs"               => \my $igs,
  "johanna"           => \my $johanna,
  "verbose|v"         => \my $verbose,
  "debug"             => \my $debug,
  "dry-run"           => \my $dryRun,
  "help|h!"           => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$grPrefix)
{
  warn "\n\n\tERROR: Missing group prefix";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

#my $quietStr = "--quiet";
my $quietStr = "";

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "\n\n\tERROR: $grDir does not exist";
  print "\n\n";
  exit 1;
}

chdir $grDir;
print "--- Changed dir to $grDir\n";

## from http://stackoverflow.com/questions/18532026/how-to-append-system-date-to-a-filename-in-perl
my @now = localtime();
my $timeStamp = sprintf("%04d-%02d-%02d_%02d_%02d_%02d",
			$now[5]+1900, $now[4]+1, $now[3],
			$now[2],      $now[1],   $now[0]);

my $mcDir  = $grPrefix . "_MC_models_dir";
my $outDir = $grPrefix . "_MC_rand_seqs_dir";


my $cmd = "rm -f $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

my $bbDir = $outDir; # bare bone out dir name (no time stamp)
$outDir =  $outDir . "_$timeStamp";

print "--- Generating random sequences of MC models\n";
$cmd = "buildMC --random-seq-length $randSeqLen --random-sample-size $nRand -t $mcDir/spp_paths.txt -k 8 -d $mcDir -o $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

$cmd = "ln -s $outDir $bbDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

if ($verbose)
{
  print "\n\tOutput written to $outDir\n\n";
}

## rsample.fa	rsample.tx are files in $outDir

print "--- Selecting only sequences from species models\n";

my $errorDir = $grPrefix . "_MC_models_clError_dir";
#$outDir = $errorDir;

my $txFile = $outDir . "/rsample.tx";
my %tx = readTbl($txFile);

# higher rank seq IDs
my %hrIDs;
my @sppIDs;
for (keys %tx)
{
  if ( $_ =~ /^\w_/ || $_ =~ /^sg_/ )
  {
    $hrIDs{$_} = 1;
  }
  else
  {
    push @sppIDs, $_;
  }
}
if ($debug)
{
  print "hrTbl:\n";
  printTbl(\%hrIDs);
  print "\n";
}

my $hrFile = $outDir . "/high_tx_ranks_rsample.seqIDs";
my @hrankIDs = keys %hrIDs;
writeArray(\@hrankIDs, $hrFile);

delete @tx{@hrankIDs};
#%tx = @tx{@sppIDs};
my $sppTxFile = $outDir . "/spp_rsample.tx";
writeTbl(\%tx, $sppTxFile);

print "--- Selecting only species sequences\n";
my $faFile    = $outDir . "/rsample.fa";
my $sppFaFile = $outDir . "/spp_rsample.fa";

$cmd = "select_seqs.pl $quietStr -e $hrFile -i $faFile -o $sppFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Running classify on $sppFaFile\n";
$cmd = "classify -d $mcDir -i $sppFaFile -o $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Comparing radn seq's taxonomy with the classification results\n";
$cmd = "cmp_tx.pl --verbose $quietStr -i $sppTxFile -j $outDir/MC_order7_results.txt -o $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


####################################################################
##                               SUBS
####################################################################

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
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

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

exit;
