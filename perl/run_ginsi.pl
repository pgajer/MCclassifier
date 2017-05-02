#!/usr/bin/env perl

=head1 NAME

  run_ginsi.pl

=head1 DESCRIPTION

  Run mafft's ginsi aligner on the degapped trimmed alignment fasta file and then
  do FastTree and rerooting.

=head1 SYNOPSIS

  run_ginsi.pl -i <input group name> -p <number of processors to use for ginsi (works only on linux)>

=head1 OPTIONS

=over

=item B<--input-group-name, -i>
  Prefix of input group. For example, if Firmicutes_group_6_dir is a directory
  Firmicutes_group_6 group, then the input group name is "Firmicutes_group_6".

=item B<--num-proc, -p>
  Number of processors to be used. Default value: 8.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--quiet>
  Do not print progress messages.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Firmicutes_dir

  run_ginsi.pl --debug -i Firmicutes_group_6 -p 8

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
  "input-group|i=s" 	=> \my $grPrefix,
  "num-proc|p=i"        => \my $nProc,
  "igs"                 => \my $igs,
  "johanna"             => \my $johanna,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "debug2"              => \my $debug2,## file name debug
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$grPrefix )
{
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $ginsi = "/usr/local/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)

my $igsStr = "";
if ( defined $igs )
{
  $igsStr = "--igs";
  $ginsi = "/usr/local/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)
}

my $johannaStr = "";
if ( defined $johanna )
{
  $johannaStr = "--johanna";
}

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

my @suffixes;

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

####################################################################
##                               MAIN
####################################################################


chdir $grDir;
print "--- Changed dir to $grDir\n";

my $trimmedAlgnFile = $grPrefix . "_algn_trimmed.fa";
my $outgroupFile    = $grPrefix . "_outgroup.seqIDs";

if ( ! -e $trimmedAlgnFile )
{
  warn "WARNING: $trimmedAlgnFile does not exist. Creating a symbolic link to $algnFile.\n";
  my $ap = abs_path( $algnFile );
  my $cmd = "rm -f $trimmedAlgnFile; ln -s $ap $trimmedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}
elsif ( ! -e $outgroupFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n\n";
  exit 1;
}


## Gathering outgroup data
print "--- Parsing $outgroupFile\n";
my @ogSeqIDs = readArray($outgroupFile);
my %ogInd = map{$_ =>1} @ogSeqIDs; # outgroup elements indicator table

print "--- Extracting seq IDs from trimmed alignment fasta file\n";
my @seqIDs = get_seqIDs_from_fa($trimmedAlgnFile);

print "--- Testing if outgroup sequences are part of seqIDs\n";
my @ogDiff = diff( \@ogSeqIDs, \@seqIDs );
@ogSeqIDs = comm( \@ogSeqIDs, \@seqIDs );

if ( scalar(@ogDiff) != 0 )
{
  warn "\n\tWARNING the following outgroup seq IDs are not in the trimmed alignment file:\n\n";
  printArray(\@ogDiff);
}

if ( scalar(@ogSeqIDs) == 0 )
{
  warn "\n\n\tERROR: All outgroup seq's were lost";
  print "\n\n";
  exit 1;
}

if ($debug)
{
  print "\n\tNumber of seq's in the trimmed alignment/lineage files: " . @seqIDs . "\n";
  print   "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n\n";
}

my $ginsiAlgnFile = $grPrefix . "_ginsi_algn.fa";
if ( ! -e $ginsiAlgnFile )
{
  print "--- Redoing alignment using ginsi\n";
  ## first remove gaps
  my $faFile2 = $grPrefix . "_from_algn_gapFree.fa";
  $cmd = "rmGaps -i $trimmedAlgnFile -o $faFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $ginsiAlgnFile; ginsi --inputorder $quietStr --thread $nProc $faFile2 > $ginsiAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

my $rrPrunedGinsiTreeFile = $grPrefix . "_ginsi_rr.tree";
if ( ! -e $rrPrunedGinsiTreeFile)
{
  print "--- Producing tree for the ginsi algn\n";
  my $ginsiTreeFile = $grPrefix . "_ginsi.tree";
  $cmd = "rm -f $ginsiTreeFile; FastTree -nt $trimmedAlgnFile > $ginsiTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;# && !$skipFastTree;

  # rerooting it
  print "--- Rerooting the ginsi tree\n";
  $cmd = "rm -f $rrPrunedGinsiTreeFile; nw_reroot $ginsiTreeFile @ogSeqIDs | nw_order -  > $rrPrunedGinsiTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

exit 0;
