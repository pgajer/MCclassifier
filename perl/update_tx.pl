#!/usr/bin/env perl

=head1 NAME

  update_tx.pl

=head1 DESCRIPTION

  The script updates genus (or any other higher than species taxonomic rank)
  taxonomy using a form of majority vote on vicut clusters.

=head1 SYNOPSIS

  update_tx.pl -a <tx file> -d <vicut directory> [Options]

=head1 OPTIONS

=over

=item B<--tx-file, -a>
  Tab or space delimited taxonomy file.

=item B<--vicut-dir, -d>
  vicut output directory, that will be used also as an output directory.

=item B<--use-long-spp-names>
  Use long species names for sequences from multi-species vicut clusters.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  update_tx.pl --debug -a genus.tx -d genus_vicut_dir

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Cwd 'abs_path';

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "tx-file|a=s"        => \my $txFile,
  "vicut-dir|d=s"      => \my $vicutDir,
  #"use-long-spp-names" => \my $useLongSppNames,
  "dry-run"            => \my $dryRun,
  "debug"              => \my $debug,
  "help|h!"            => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$txFile)
{
  print "\n\nERROR: Missing taxonomy file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$vicutDir)
{
  print "\n\nERROR: Missing vicut directory\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -f $txFile )
{
  print "\n\nERROR: $txFile does not exist\n\n\n";
  exit 1;
}

my $newTxFile = "$vicutDir/minNodeCut.cltrs";

if ( ! -f $newTxFile )
{
  print "\n\nERROR: $newTxFile does not exist\n\n\n";
  exit 1;
}

####################################################################
##                               MAIN
####################################################################


print "--- Parsing taxonomy table (in $txFile) from BEFORE vicut taxonomy modifications\n";
my %preVicutTxTbl = read2colTbl($txFile);

#print "\t- Updating min-node-cut taxonomy using majority vote\n";

my $cltrFile = "$vicutDir/minNodeCut.cltrs";
my ($rclTbl, $rclIDsTbl, $rtxTbl)  = readCltrTbl($cltrFile);

my %clIDsTbl  = %{$rclIDsTbl};    # cltrID => ref to array of IDs in the given cluster
my %clTbl     = %{$rclTbl};       # seqID => cltrID
my %txTbl     = %{$rtxTbl};       # seqID => taxonomy (NA for query seq's; otherwise should be consistent/the same as in preVicutTxTbl)

# Cluster frequency table
my %clTxFreqTbl;
for my $id ( keys %txTbl )
{
  my $tx = $txTbl{$id};
  my $cl = $clTbl{$id};
  $clTxFreqTbl{$cl}{$tx}++;
}

## Renaming of cluster members strategy

## If a cluster consists of at least one non-NA taxon, all seq's of that cluster
## are relabeled to the name of that taxon.

## There is only the NA taxon in the cluster (meaning this is a cluster
## consisting of query seq's). In this case, look at the frequency table of taxon
## names from before vicut run and pick as the new name of all seq's of that
## cluster, the name of the most abundant old taxon.

if ( $debug )
{
  print "\nFinding new taxonomy for all seq's of each vicut cluster\n";
  print "Going over clusters using clTxFreqTbl table\n";
}

my %newTxTbl; # new taxonomy table
for my $cl ( keys %clTxFreqTbl )
{
  # Sort taxons within a cluster by the number of seq's representing this taxon within the cluster
  # or if they have the same number of seq's, by the alphabetical order of their names
  my @txs = sort { $clTxFreqTbl{$cl}{$b} <=> $clTxFreqTbl{$cl}{$a} || $a cmp $b } keys %{$clTxFreqTbl{$cl}};

  if ($debug)
  {
    print "\nCluster $cl taxons:\n";
    map {print "\t$_\t" . $clTxFreqTbl{$cl}{$_} . "\n"} @txs;
    print "\n";
  }

  ## checking if the cluster has one or more members
  if ( @txs == 1 )
  {
    my $tx = $txs[0];
    print "tx: $tx\n" if $debug;

    if ( $tx eq "NA" )
    {
      print "The only member of the cluster is NA\n" if $debug;

      my @ids = @{$clIDsTbl{$cl}};
      my %preVicutTxFreqTbl;
      for ( @ids )
      {
	$preVicutTxFreqTbl{ $preVicutTxTbl{$_} }++;
      }

      my @preVicutTxs = sort { @{$preVicutTxFreqTbl{$b}} <=> @{$preVicutTxFreqTbl{$a}} } keys %preVicutTxFreqTbl;

      if ( $debug )
      {
	print "Pre-vicut taxons of the cluster\n";
	printFormatedTbl( \%preVicutTxFreqTbl, \@preVicutTxs );
	print "\n";
      }

      my $newTx = $preVicutTxs[0];
      print "New taxonomy: $newTx\n" if $debug;
      for my $id ( @{$clIDsTbl{$cl}} )
      {
	$newTxTbl{$id} = $newTx;
      }
    }
    else
    {
      ## The only taxon is not an NA.
      ## No change in taxonomy. All seq's of that cluster preserve their old
      ## taxonomy - the one of the only member of the cluster.
      for my $id ( @{$clIDsTbl{$cl}} )
      {
	$newTxTbl{$id} = $tx;
      }
    }
  }
  else
  {
    ## There is more than one taxon in the cluster.
    ## The non-NA cluster with the largest number of seq's is the winner.
    my $newTx = $txs[0];
    if ( $newTx eq "NA" )
    {
      $newTx = $txs[1];
    }

    print "New taxonomy: $newTx\n" if $debug;
    for my $id ( @{$clIDsTbl{$cl}} )
    {
      $newTxTbl{$id} = $newTx;
    }
  }

} # end of if ( @txs > 1 )

print "\n\n" if $debug;


# Frequencies of taxons in the new taxonomy table
my %txClFreqTbl;
my %txClIDsTbl;
for my $id ( keys %newTxTbl )
{
  my $tx = $newTxTbl{$id};
  my $cl = $clTbl{$id};
  $txClFreqTbl{$tx}{$cl}++;
  push @{ $txClIDsTbl{$tx}{$cl} }, $id;
}

# Defining multiplicity index for taxons appearing more than once

## For each taxon sort clusters by size and assign multiplicity index to each
## cluster (1= largest, 2=2nd largest, etc) Only for taxons with more than one
## cluster.

my %txClIdxTbl;
for my $tx ( keys %txClFreqTbl )
{
  my @cls = sort { $txClFreqTbl{$tx}{$b} <=> $txClFreqTbl{$tx}{$a} } keys %{ $txClFreqTbl{$tx} };
  if ( @cls > 1 )
  {
    my $idx = 1;
    for my $cl ( @cls )
    {
      $txClIdxTbl{$tx}{$cl} = $idx;
      $idx++;
    }
  }
  else
  {
    my $cl = shift @cls;
    $txClIdxTbl{$tx}{$cl} = 0;
  }
}


if ( $debug )
{
  print "\nTaxon multiplicity indices\n";
  print "Taxon\tMult\tcltrID\n";
  for my $tx ( keys %txClFreqTbl )
  {
    my @cls = sort { $txClFreqTbl{$tx}{$b} <=> $txClFreqTbl{$tx}{$a} } keys %{ $txClFreqTbl{$tx} };
    for my $cl ( @cls )
    {
      print "$tx\t$cl\t" . $txClIdxTbl{$tx}{$cl} . "\n";
    }
  }
  print "\n";
}


my $updatedTxFile = "$vicutDir/TDupdated.tx";
writeTbl2(\%newTxTbl, $updatedTxFile);

my $updatedTxFile2 = "$vicutDir/TDupdated2.tx";
open OUT, ">$updatedTxFile2" or die "Cannot open $updatedTxFile2 for writing: $OS_ERROR\n";
for my $tx ( keys %txClFreqTbl )
{
  my @cls = sort { $txClFreqTbl{$tx}{$b} <=> $txClFreqTbl{$tx}{$a} } keys %{ $txClFreqTbl{$tx} };
  for my $cl ( @cls )
  {
    my @ids = @{ $txClIDsTbl{$tx}{$cl} };
    for my $id ( @ids )
    {
      print OUT "$id\t$tx\t" . $txClIdxTbl{$tx}{$cl} . "\n";
    }
  }
}
close OUT;

print "\tUpdated taxonomy  written to $updatedTxFile2\n\n" if $debug;


####################################################################
##                               SUBS
####################################################################

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# read cltrs and new species assignment
# elements of the first column to the second column
# remove the header line
sub readCltrTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %clTbl;      # seqID => cltrID
  my %txTbl;      # seqID => taxonomy
  my %clIDsTbl;   # cltrID => ref to array of IDs in the given cluster
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  my $header = <IN>;
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    push @{$clIDsTbl{$cl}}, $id;
    $txTbl{$id} = $tx;
    $clTbl{$id} = $cl;
  }
  close IN;

  return (\%clTbl, \%clIDsTbl, \%txTbl);
}

# write hash table to a file
sub writeTbl{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} keys %tbl;
  close OUT;
}

# write hash table to a file
sub writeTbl2{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  for (keys %tbl)
  {
    if (exists $tbl{$_})
    {
      print OUT $_ . "\t" . $tbl{$_} . "\n"
    }
    else
    {
      warn "\n\n\tERROR: $_ does not exist in the table to be printed to $outFile";
      print "\n\n";
      exit 1;
    }
  }
  close OUT;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
  print "\n\n";
}


# print elements of a hash table
sub printTbl
{
  my ($rTbl, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
  print "\n\n";
}

# count L crispatus seq's
sub countLc
{
  my $rTbl = shift;
  my $i = 0;
  for (keys %$rTbl)
  {
    if( $rTbl->{$_} =~ /crispatus/ )
    {
      $i++;
    }
  }
  return $i;
}

# common part of two arrays
sub comm
{
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
}

exit 0;
