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

=item B<--tx-file, -i>
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
  exit;
}

if (!$txFile)
{
  print "\n\nERROR: Missing taxonomy file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$vicutDir)
{
  print "\n\nERROR: Missing vicut directory\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if ( ! -f $txFile )
{
  print "\n\nERROR: $txFile does not exist\n\n\n";
  exit;
}

my $newTxFile = "$vicutDir/minNodeCut.cltrs";

if ( ! -f $newTxFile )
{
  print "\n\nERROR: $newTxFile does not exist\n\n\n";
  exit;
}

####################################################################
##                               MAIN
####################################################################


print "--- Parsing tax table from before vicut taxonomy modifications\n";
my %tx = read2colTbl($txFile);

print "\t- Updating min-node-cut taxonomy using majority vote\n";

my $cltrFile = "$vicutDir/minNodeCut.cltrs";
if ( ! -f $cltrFile )
{
  print "\n\nERROR: $cltrFile does not exist\n\n\n";
  exit;
}

my ($rid2clTbl, $rclTbl, $rclFreqTbl, $rtxTbl)  = readCltrTbl($cltrFile);

##my %clFreqTbl = %{$rclFreqTbl}; # cltrID => table of taxon frequencies within the given cluster
my %clTbl     = %{$rclTbl};     # cltrID => ref to array of IDs in the given cluster
my %txTbl     = %{$rtxTbl};     # seqID => taxonomy
my %id2clTbl  = %{$rid2clTbl};  # seqID => cltrID

## Generating updated taxonomy file
## changing NAs to taxonomy in the original taxonomy table
for my $id ( keys %txTbl )
{
  if ( $txTbl{$id} eq "NA" )
  {
    if ( exists $tx{$id} )
    {
      $txTbl{$id} = $tx{$id};
    }
    else
    {
      warn "\n\tWARNING: $id not present in $txFile table";
      print "\n";
    }
  }
}

## Making sure there are no NA species
## Recall that NAs in $cltrFile come from query sequences
## the above loop should have replaced all NAs with their taxonomy as supplied to vicut
for my $id ( keys %txTbl )
{
  my $sp = $txTbl{$id};
  if ($sp eq "NA")
  {
    warn "\n\n\tERROR: Discovered NA species";
    print "\n\n";
    exit;
  }
}

if ($debug)
{
  my %txFreq;
  for my $id ( keys %txTbl )
  {
    my $sp = $txTbl{$id};
    $txFreq{$sp}++;
  }

  printTbl(\%txFreq, "txFreq");
}

my %txIDs;
for my $id ( keys %txTbl )
{
  my $sp = $txTbl{$id};
  push @{$txIDs{$sp}}, $id;
}


# Updating %clFreqTbl after the above modification
my %clFreqTbl;
for my $id ( keys %txTbl )
{
  my $sp = $txTbl{$id};
  my $cl = $id2clTbl{$id};
  $clFreqTbl{$cl}{$sp}++;
}

# Changing taxonomy of multi-species clusters based on majority vote

for my $cl ( keys %clFreqTbl )
{
  my @txs = sort { $clFreqTbl{$cl}{$b} <=> $clFreqTbl{$cl}{$a} || $a cmp $b } keys %{$clFreqTbl{$cl}};
  if ($debug)
  {
    print "Cluster $cl taxons:\n";
    map {print "\t$_\t" . $clFreqTbl{$cl}{$_} . "\n"} @txs;
    print "\n";
  }

  if ( @txs > 1 )
  {
    my $newTx = shift @txs;
    for (@txs)
    {
      $newTx .= "_" . $_;
    }

    if ($debug)
    {
      print "\tChanging taxonomy of all sequences of the cluster to $newTx\n\n";
    }
    for my $id ( @{$clTbl{$cl}} )
    {
      $txTbl{$id} = $newTx;
    }

  } # end of if ( @txs > 1 )
}

print "\n\n" if $debug;


my $updatedTxFile = "$vicutDir/updated.tx";
writeTbl2(\%txTbl, $updatedTxFile);

#print "\tUpdated taxonomy  written to $updatedTxFile\n\n";


####################################################################
##                               SUBS
####################################################################

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

  my $file = shift;

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

# read cltrs and new species assignment
# elements of the first column to the second column
# remove the header line
sub readCltrTbl
{
  my $file = shift;

  my %txTbl;     # seqID => taxonomy
  my %id2clTbl;  # seqID => cltrID
  my %clTbl;     # cltrID => ref to array of IDs in the given cluster
  my %clFreqTbl; # cltrID => table of taxon frequencies within the given cluster
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  my $header = <IN>;
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    push @{$clTbl{$cl}}, $id;
    $clFreqTbl{$cl}{$tx}++;
    $txTbl{$id} = $tx;
    $id2clTbl{$id} = $cl;
  }
  close IN;

  return (\%id2clTbl, \%clTbl, \%clFreqTbl, \%txTbl);
}

# write hash table to a file
sub writeTbl{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} keys %tbl;
  close OUT;
}

# write hash table to a file
sub writeTbl2{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  for (keys %tbl)
  {
    if (exists $tbl{$_})
    {
      print OUT $_ . "\t" . $tbl{$_} . "\n"
    }
    else
    {
      warn "\n\n\tERROR: $_ does not exist in the table to be printed to $outFile\n\n";
      exit;
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

exit;
