#!/usr/bin/env perl

=head1 NAME

  split_taxon.pl

=head1 DESCRIPTION

  Given a condenset tree at some taxonomic rank (say species) and a percentile
  threshold, the script performs a partition of a phylogenetic tree with the
  given percentile distance threshold using phyloPart and then does vicut
  clustering using all except 0 clusters of phyloPart as annotation and elements
  of cluster 0 as query leaves.

  If a lineage file is also give, the above clustering is added to the lineage
  data, creating sub-parent-rank (in the above example, sub-genus) taxonomy of
  the form, Lactobacillus_sub_3.

=head1 SYNOPSIS

  split_taxon.pl -i <group Prefix> [Options]

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Group prefix.

=item B<--taxon, -t>
  Taxon to be split. Expected is one of the following strings: "spp", "genus", "family", "order", "class"

=item B<--perc-thld, -p>
  Percentile threshold specified as a decimal between 0 and 1.

=item B<--tree-file, -j>
  A phylogenetic tree file in the Newick format.

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

  split_taxon.pl --debug -i Firmicutes_group_6_V3V4 -p 0.10 -t spp

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
  "group-prefix|i=s"     => \my $grPrefix,
  "taxon|t=s"            => \my $taxon,
  "perc-thld|p=f"        => \my $percThld,
  "tree-file|j=s"        => \my $treeFile,
  "igs"                  => \my $igs,
  "johanna"              => \my $johanna,
  "verbatim|v"           => \my $verbatim,
  "debug"                => \my $debug,
  "dry-run"              => \my $dryRun,
  "help|h!"              => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$grPrefix)
{
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$percThld)
{
  warn "\n\n\tERROR: Missing percentile threshold";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$taxon)
{
  warn "\n\n\tERROR: Missing taxon threshold";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my $phyloPart = "/Users/pgajer/devel/MCclassifier/PhyloPart_v2.1/PhyloPart_v2.1.jar";

if ($igs)
{
  $phyloPart = "???";
}

if ($johanna)
{
  $phyloPart = "???";
}

if ( ! -e $phyloPart )
{
  warn "\n\n\tERROR: $phyloPart does not exist";
  print "\n\n";
  exit;
}

my $debugStr = "";
$debugStr = "--debug" if $debug;

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "\n\n\tERROR: $grDir does not exist";
  print "\n\n";
  exit;
}

$grPrefix = "$grDir/$grPrefix";

if ( !defined $treeFile )
{
  $treeFile = $grPrefix . "_final_" . $taxon . "_condensed.tree";
  if ( ! -e $treeFile )
  {
    warn "\n\n\tERROR: Tree file $treeFile does not exist";
    print "\n\n";
    exit;
  }
}

my $lineageFile = $grPrefix . "_final.lineage";
if ( ! -f $lineageFile )
{
  warn "ERROR: $lineageFile does not exist\n";
  exit;
}


####################################################################
##                               MAIN
####################################################################

print "--- Running phylo partitioning on $treeFile at $percThld percentile thld\n";
my $partFile     = $grPrefix . "_phyloPart_$percThld" . ".txt";
my $phyloPartLog = $grPrefix . "_phyloPart.log";
my $cmd = "java -jar $phyloPart $treeFile $percThld -o$partFile > $phyloPartLog 2>&1";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Parsing phylo partitioning data\n";
my %part = read_part_tbl($partFile);

my $partCltrFile = $grPrefix . "_phyloPart_$percThld" . ".cltr";
print "--- Writing phylo partitioning to $partCltrFile\n";
writeTbl(\%part, $partCltrFile);

print "--- Generating annotation and query files\n";
my $annFile    = $grPrefix . "_phyloPart_$taxon" . "_ann.tx";
my $queryFile  = $grPrefix . "_phyloPart_$taxon" . "_query.seqIDs";
my $vicutDir   = $grPrefix . "_phyloPart_$taxon" . "_vicut_dir";
my $nQuerySeqs = 0;
my $nAnnSeqs   = 0;
my @queryTx;
open QOUT, ">$queryFile" or die "Cannot open $queryFile for writing: $OS_ERROR\n";
open ANNOUT, ">$annFile" or die "Cannot open $annFile for writing: $OS_ERROR\n";
for my $id ( keys %part )
{
  my $cl = $part{$id};

  if ( $cl == 0 )
  {
    print QOUT "$id\n";
    $nQuerySeqs++;
    push @queryTx, $id;
  }
  else
  {
    print ANNOUT "$id\t$cl\n";
    $nAnnSeqs++;
  }
}
close ANNOUT;
close QOUT;

print "--- Parsing tree leaves\n";
my $treeLeavesFile = "$grPrefix" . "_tree.leaves";
$cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my @treeLeaves = readArray($treeLeavesFile);

print "\n\n\tNumber of annotation seq's: $nAnnSeqs\n";
print     "\tNumber of query seq's:      $nQuerySeqs\n";
print     "\tSum:                        " . ($nAnnSeqs + $nQuerySeqs) . "\n";
print     "\tNumber of leaves: " . @treeLeaves . "\n\n";

printArray(\@queryTx, "Query taxons:\n");

print "--- Running vicut\n";
if ($nQuerySeqs)
{
  $cmd = "vicut -t $treeFile -a $annFile -q $queryFile -o $vicutDir";
}
else
{
  $cmd = "vicut -t $treeFile -a $annFile -o $vicutDir";
}
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running update_tx.pl to cluster singletons with larger clusters\n";
$cmd = "update_tx.pl $debugStr -a $partCltrFile -d $vicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $cltrFile = "$vicutDir/updated.tx";
my %part2 = readTbl($cltrFile);

my %cltr;
for my $tx (keys %part2)
{
  my $cl = $part2{$tx};
  push @{$cltr{$cl}}, $tx;
}

print "Phylo partition sizes\n";
my %partFreq; ## number of elements per phylo partition cluster
map { $partFreq{$_}++ } values %part;

printFormatedTbl(\%partFreq);

print "\nVicut updated phylo partition clusters:\n";
my @q = sort { @{$cltr{$b}} <=> @{$cltr{$a}} } keys %cltr;
for my $cl (@q)
{
  print "Cluster $cl:\n";
  for (@{$cltr{$cl}})
  {
    print "\t$_\n";
  }
}
print "\n\n";




####################################################################
##                               SUBS
####################################################################

# print elements of a hash table so that arguments are aligned
sub printFormatedTbl
{
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

# print elements of a hash table so that arguments are aligned
sub printFormatedTblToFile
{
  my ($rTbl, $rSub, $fh, $msg) = @_; # the second argument is a subarray of the keys of the table

  print $fh "\n$msg\n";

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
    print $fh "$_$pad" . $rTbl->{$_} . "\n";
  }
  print $fh "\n";
}

# parse a CSV partition table
sub read_part_tbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  my $headerStr = <IN>;
  foreach my $line (<IN>)
  {
    chomp $line;

    ##  $ clustername        : int  426 426 432 432 432 432 432 432 432 449 ...
    ##  $ bootstrap          : num  0.904 0.904 0.908 0.908 0.908 0.908 0.908 0.908 0.908 0.976 ...
    ##  $ leafname           : chr  "Lactobacillus_hordei" "Lactobacillus_mali_tRT_2" "Lactobacillus_nagelii" "Lactobacillus_vini" ...
    ##  $ branchPath         : num  0.0462 0.0525 0.0547 0.0546 0.0526 ...
    ##  $ medianOfDistances  : num  0.00651 0.00651 0.01502 0.01502 0.01502 ...
    ##  $ sequencesperCluster: int  2 2 7 7 7 7 7 7 7 2 ...

    my @f = split ",", $line;
    my $cltrId = shift @f;
    my $boot   = shift @f;
    my $leafId = shift @f;
    $tbl{ $leafId } = $cltrId;
  }
  close IN;

  return %tbl;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# read table with one column
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    print "\n\nERROR in readArray() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit;
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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR in readTbl() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit;
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


exit;
