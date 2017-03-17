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
  split_taxon.pl --debug -i Firmicutes_group_0_V3V4 -p 0.10 -t genus

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

my $readNewickFile = "/Users/pgajer/.Rlocal/read.newick.R";
my $phyloPart = "/Users/pgajer/devel/MCclassifier/PhyloPart_v2.1/PhyloPart_v2.1.jar";

if ($igs)
{
  $phyloPart = "???";
  $readNewickFile = "???";
}

if ($johanna)
{
  $phyloPart = "/Users/jholm/bin/PhyloPart_v2.1/PhyloPart_v2.1.jar";
  $readNewickFile = "";
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
  warn "\n\n\tERROR: $lineageFile does not exist";
  print "\n\n";
  exit;
}

my $newLineageFile = $grPrefix . "_final2.lineage";


####################################################################
##                               MAIN
####################################################################

print "--- Parsing lineage table\n";
my %lineageTbl = read2colTbl($lineageFile);

my %parent;
for my $id ( keys %lineageTbl )
{
  my $lineage = $lineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  my $or = pop @f;
  my $cl = pop @f;
  my $ph = pop @f;

  $parent{$sp} = $ge;
  $parent{$ge} = $fa;
  $parent{$fa} = $or;
  $parent{$or} = $cl;
  $parent{$cl} = $ph;
  $parent{$ph} = "d_Bacteria";
}

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

my $cltrFile = "$vicutDir/minNodeCut.cltrs";
print "--- Parsing vicut clustering file $cltrFile\n";
if ( ! -f $cltrFile )
{
  warn "\nERROR: $cltrFile does not exist";
  exit;
}

my %part2;  # seqID => cltrID
my %cltr;     # cltrID => ref to array of IDs in the given cluster
open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR\n";
my $header = <IN>;
foreach (<IN>)
{
  chomp;
  my ($id, $cl, $tx) = split /\s+/,$_;
  $part2{$id} = $cl;
  push @{$cltr{$cl}}, $id;
}
close IN;

# $cmd = "update_tx.pl $debugStr -a $partCltrFile -d $vicutDir";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# my $cltrFile = "$vicutDir/updated.tx";
# my %part2 = readTbl($cltrFile);

# my %cltr;
# for my $tx (keys %part2)
# {
#   my $cl = $part2{$tx};
#   push @{$cltr{$cl}}, $tx;
# }

print "Phylo partition sizes\n";
my %partFreq; ## number of elements per phylo partition cluster
map { $partFreq{$_}++ } values %part;

printFormatedTbl(\%partFreq);

print "After vicut phylo partition sizes\n";
my %part2Freq; ## number of elements per phylo partition cluster
map { $part2Freq{$_}++ } values %part2;

printFormatedTbl(\%part2Freq);


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


if ($taxon eq "spp")
{
  ## Each cluster will have the name sub_<genus>_<idx>
  ## if all species of the cluster are from the same genus <genus>
  ## and
  ## sub_<genus_1>_<genus_2>_...<genus_n>_<idx>
  ## if there is more than one genus within the cluster

  my %genera;
  my %clGenus;
  for my $cl (keys %cltr)
  {
    #print "Cluster $cl:\n";
    my %locGenera;
    for (@{$cltr{$cl}})
    {
      my ($g, $suffix) = split "_", $_;
      $locGenera{$g}++;
    }

    my @gen = sort { $locGenera{$b} <=> $locGenera{$a} } keys %locGenera;
    my $genStr = join "_", @gen;
    $genera{$genStr}++;
    $clGenus{$cl} = $genStr;
  }

  print "\nDiscovered Cluster Genera\n";
  my @a = sort { $genera{$b} <=> $genera{$a} } keys %genera;
  for (@a)
  {
    print "\t$_\n";
  }
  print "\n";

  print "--- Changing cluster names\n";
  my %genusIdx;
  my %cltr2;
  my %spSubGenus;
  my %spSubGenusIdx; # this is for species_idx tree only
  my @q = sort { @{$cltr{$b}} <=> @{$cltr{$a}} } keys %cltr;
  my $count = 1;
  for my $cl (@q)
  {
    my $tx = $clGenus{$cl};
    $genusIdx{$tx}++;
    if ( $genera{$tx}>1 )
    {
      $tx = "sub_$tx" . "_$genusIdx{$tx}";
    }
    else
    {
      $tx = "sub_$tx";
    }
    $cltr2{$tx} = $cltr{$cl};
    print "\tProcessing cluster $cl with new taxonomy $tx\n";
    for (@{$cltr{$cl}})
    {
      $spSubGenus{$_} = $tx;
      $spSubGenusIdx{$_} = $count;
    }
    $count++;
  }

  print "\nVicut updated phylo partition clusters with new names:\n";
  @q = sort { @{$cltr2{$b}} <=> @{$cltr2{$a}} } keys %cltr2;
  for my $cl (@q)
  {
    print "Cluster $cl:\n";
    for (@{$cltr2{$cl}})
    {
      print "\t$_\n";
    }
  }
  print "\n\n";

  print "--- Creating tree with taxon_cluster leaf names\n";
  my $spClFile = $grPrefix . ".sppCl";
  open OUT, ">$spClFile" or die "Cannot open $spClFile for writing: $OS_ERROR\n";
  for (keys %spSubGenusIdx)
  {
    print OUT "$_\t$_" . "_cl" . $spSubGenusIdx{$_} . "|\n";
  }
  close OUT;

  my $spClFileAbsPath = abs_path( $spClFile );
  writeTbl(\%spSubGenusIdx, $spClFileAbsPath);
  print "spSubGenusIdx written to $spClFileAbsPath\n" if $debug;

  my $treeFile2 = $grPrefix . "_final_" . $taxon . "_condensed_cltrs.tree";
  $cmd = "nw_rename $treeFile $spClFile | nw_order -c n  - > $treeFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "---Updating lineage table\n";
  my %lineageTbl2;
  for my $id ( keys %lineageTbl )
  {
    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $subGenus = $spSubGenus{$sp};
    my @t = (@f, $subGenus, $sp);
    $lineageTbl2{$id} = join ";", @t;
  }

  open OUT, ">$newLineageFile" or die "Cannot open $newLineageFile for writing: $OS_ERROR\n";
  for my $id (keys %lineageTbl2)
  {
    my $lineage = $lineageTbl2{$id};
    print OUT "$id\t$lineage\n";
  }
  close OUT;

  print "\n\n\tUpdated lineage written to $newLineageFile\n";
  print "\tSpecies_cluster leaves tree written to $treeFile2\n\n";
}
else
{
  ## Each cluster will have the name sub_<parent>_<idx>
  ## if all taxons of the cluster have the same parent
  ## and
  ## sub_<perent_1>_<parent_2>_...<parent_n>_<idx>
  ## otherwise

  my %pars;
  my %clParent;
  for my $cl (keys %cltr)
  {
    #print "Cluster $cl:\n";
    my %locPars;
    for (@{$cltr{$cl}})
    {
      my $p = $parent{$_};
      $p =~ s/_tRT_\d+//;
      my @f = split "_", $p;
      for (@f)
      {
	$locPars{$_}++;
      }
    }

    my @par = sort { $locPars{$b} <=> $locPars{$a} } keys %locPars;
    my $parStr = join "_", @par;
    $pars{$parStr}++;
    $clParent{$cl} = $parStr;
  }

  print "\nDiscovered Cluster Parents\n";
  my @a = sort { $pars{$b} <=> $pars{$a} } keys %pars;
  for (@a)
  {
    print "\t$_\n";
  }
  print "\n";

  print "--- Changing cluster names\n";
  my %parentIdx;
  my %cltr2;
  my %txSubParent;
  my %parentIdx;
  my %txSubParentIdx; # this is for species_idx tree only
  my @q = sort { @{$cltr{$b}} <=> @{$cltr{$a}} } keys %cltr;
  my $count = 1;
  for my $cl (@q)
  {
    my $tx = $clParent{$cl};
    $parentIdx{$tx}++;
    if ( $pars{$tx}>1 )
    {
      $tx = "sub_$tx" . "_$parentIdx{$tx}";
    }
    else
    {
      $tx = "sub_$tx";
    }
    $cltr2{$tx} = $cltr{$cl};
    $parentIdx{$tx} = $count;
    print "\tProcessing cluster $cl with new taxonomy $tx\n";
    for (@{$cltr{$cl}})
    {
      $txSubParent{$_} = $tx;
      $txSubParentIdx{$_} = $count;
    }
    $count++;
  }

  print "\nVicut updated phylo partition clusters with new names:\n";
  @q = sort { @{$cltr2{$b}} <=> @{$cltr2{$a}} } keys %cltr2;
  for my $cl (@q)
  {
    print "Cluster $cl => $parentIdx{$cl}:\n";
    for (@{$cltr2{$cl}})
    {
      print "\t$_\n";
    }
  }
  print "\n\n";

  print "--- Creating tree with taxon_cluster leaf names\n";
  my $spClFile = $grPrefix . ".sppCl";
  my $spClFile2 = abs_path( $grPrefix . ".sppCl2" );
  open OUT, ">$spClFile" or die "Cannot open $spClFile for writing: $OS_ERROR\n";
  open OUT2, ">$spClFile2" or die "Cannot open $spClFile2 for writing: $OS_ERROR\n";
  for (keys %txSubParentIdx)
  {
    print OUT "$_\t$_" . "_cl" . $txSubParentIdx{$_} . "|\n";
    print OUT2 "$_" . "_cl" . $txSubParentIdx{$_} . "|\t" . $txSubParentIdx{$_} . "\n";
  }
  close OUT;
  close OUT2;
  print "\n\n--> Cluster index tbl written to $spClFile2\n" if $debug;

  # my $spClFileAbsPath = abs_path( $spClFile );
  # writeTbl(\%txSubParentIdx, $spClFileAbsPath);
  # print "\n\n--> txSubParentIdx written to $spClFileAbsPath\n" if $debug;

  my $treeFile2 = $grPrefix . "_final_" . $taxon . "_condensed_cltrs.tree";
  $cmd = "nw_rename $treeFile $spClFile | nw_order -c n  - > $treeFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "---Updating lineage table\n";
  my %lineageTbl2;
  if ($taxon eq "genus")
  {
    for my $id ( keys %lineageTbl )
    {
      my $lineage = $lineageTbl{$id};
      my @f = split ";", $lineage;
      my $sp = pop @f;
      my $ge = pop @f;
      my $subPar = $txSubParent{$ge};
      my @t = (@f, $subPar, $ge, $sp);
      $lineageTbl2{$id} = join ";", @t;
    }
  }

  open OUT, ">$newLineageFile" or die "Cannot open $newLineageFile for writing: $OS_ERROR\n";
  for my $id (keys %lineageTbl2)
  {
    my $lineage = $lineageTbl2{$id};
    print OUT "$id\t$lineage\n";
  }
  close OUT;

  print "\n\n\tUpdated lineage written to $newLineageFile\n";
  print "\tTaxon_cluster leaves tree written to $treeFile2\n\n";

  my $pdfTreeFile = abs_path( $grPrefix . "_final_" . $taxon . "_condensed_cltrs_tree.pdf" );
  my $treeFile2AbsPath = abs_path( $treeFile2 );

  plotTree($treeFile2AbsPath, $spClFile2, $pdfTreeFile);

  if ( $OSNAME eq "darwin")
  {
    $cmd = "open $pdfTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }
}

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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
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

sub plotTree
{
  my ($treeFile, $clFile, $pdfFile) = @_;

  my $Rscript = qq~

clTbl <- read.table(\"$clFile\", header=F)
str(clTbl)

cltr <- clTbl[,2]
names(cltr) <- clTbl[,1]

source(\"$readNewickFile\")
require(phytools)

tr1 <- read.newick(file=\"$treeFile\")
tr1 <- collapse.singles(tr1)

tip.colors <- cltr[tr1\$tip.label]

pdf(\"$pdfFile\", width=6, height=12)
op <- par(mar=c(0,0,0,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr1,type=\"phylogram\", tip.color=tip.colors, no.margin=FALSE, show.node.label=F)
par(op)
dev.off()
~;

  runRscript( $Rscript );
}

  # execute an R-script
sub runRscript{

  my $Rscript = shift;

  my $outFile = "rTmp.R";
  open OUT, ">$outFile",  or die "cannot write to $outFile: $!\n";
  print OUT "$Rscript";
  close OUT;

  my $cmd = "R CMD BATCH $outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

  my $outR = $outFile . "out";
  open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
  my $exitStatus = 1;

  foreach my $line (<IN>)
  {
    if ( $line =~ /Error/ )
    {
      print "R script crashed at\n$line";
      print "check $outR for details\n";
      $exitStatus = 0;
      exit;
    }
  }
  close IN;
}

exit;
