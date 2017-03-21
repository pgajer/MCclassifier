#!/usr/bin/env perl

=head1 NAME

  cluster_taxons.pl

=head1 DESCRIPTION

  Given a condenset tree at some taxonomic rank (say species), a percentile
  threshold and a taxonomic parent table, the script performs a partition of a
  phylogenetic tree with the given percentile distance threshold using phyloPart
  and then does vicut clustering using all except 0 clusters of phyloPart as
  annotation and elements of cluster 0 as query leaves and uses the taxonomic
  parent table to assign taxonomically relevant names to the clusters.

=head1 SYNOPSIS

  cluster_taxons.pl -i <tree file> -p <perc thld> -f <tx parent file> -t <taxon> -o <output file> [Options]

=head1 OPTIONS

=over

=item B<--tree-file, -i>
  A phylogenetic tree file in the Newick format.

=item B<--perc-thld, -p>
  Percentile threshold specified as a decimal between 0 and 1.

=item B<--taxon, -t>
  Taxon to be split. Expected is one of the following strings: "species", "genus", "family", "order", "class"
  The only purpose of taxon parameter is to label phyloPart related input and output data/files.

=item B<--parent-file, -f>
  A table with <taxon> => <taxonomic parent> table for the given taxonomic rank (specified by --taxon flag).

=item B<--output-file, -o>
  Output file with <taxon> => <cluster name> table.

=item B<--show-tree>
  Open the pdf file with the tree used to do clustering.

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

  cluster_taxons.pl --debug -i Firmicutes_group_6_V3V4 -p 0.10 -t spp
  cluster_taxons.pl --debug -i Firmicutes_group_0_V3V4 -p 0.10 -t genus

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
  "tree-file|i=s"        => \my $treeFile,
  "perc-thld|p=f"        => \my $percThld,
  "taxon|t=s"            => \my $taxon,
  "parent-file|f=s"      => \my $parentFile,
  "output-file|o=s"      => \my $outFile,
  "show-boot-vals"       => \my $showBoostrapVals,
  "show-tree"            => \my $showTree,
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

if (!$treeFile)
{
  warn "\n\n\tERROR: Missing tree file";
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
elsif (!$parentFile)
{
  warn "\n\n\tERROR: Missing taxonomic parent table file";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outFile)
{
  warn "\n\n\tERROR: Missing output file name";
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
elsif ( ! -e $treeFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n\n";
  exit;
}
elsif ( ! -e $parentFile )
{
  warn "\n\n\tERROR: $parentFile does not exist";
  print "\n\n";
  exit;
}

my $debugStr = "";
$debugStr = "--debug" if $debug;

# if ( ! -d $outDir )
# {
#   my $cmd = "mkdir -p $outDir";
#   print "\tcmd=$cmd\n" if $dryRun || $debug;
#   system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
# }

# my $newLineageFile = $grPrefix . "_final2.lineage";

if (defined $showBoostrapVals)
{
  $showBoostrapVals = "T";
}
else
{
  $showBoostrapVals = "F";
}

my @inputTaxons = ("species", "genus", "family", "order", "class");
my %inputTaxonsTbl = map{$_ => 1} @inputTaxons;
if ( ! exists $inputTaxonsTbl{$taxon} )
{
  warn "\n\n\tERROR: $taxon has to be one of the following strings: @inputTaxons";
  print "\n\n";
  exit;
}

####################################################################
##                               MAIN
####################################################################

print "--- Parsing parent table\n";
my %parent = read2colTbl($parentFile);

my %parentFreq; ## number of elements of parent clusters
map { $parentFreq{$_}++ } values %parent;

my %parentTbl; ## parent => children
map { push @{$parentTbl{$parent{$_}}}, $_ } keys %parent;

my $nParents = keys %parentFreq;

if ($debug)
{
  print "\n\nNumber of parents: $nParents\n";

  print "\nParent table cluster sizes:\n";
  printFormatedTbl(\%parentFreq);

  print "\nParent table clusters:\n";
  for my $cl (keys %parentTbl)
  {
    print "$cl:\n";
    for ( @{$parentTbl{$cl}} )
    {
      print "\t$_\n";
    }
  }
  print "\n";
}

print "--- Rooting the tree\n";
my $cmd = "root_tree.pl -i $treeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running phylo partitioning on $treeFile at $percThld percentile thld\n";
my $partFile     = "phyloPart_$taxon" . "_$percThld" . ".txt";
my $phyloPartLog = "phyloPart_$taxon.log";
$cmd = "rm -f $phyloPartLog; java -jar $phyloPart $treeFile $percThld -o$partFile > $phyloPartLog ";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Parsing phylo partitioning data\n";
my %part = read_part_tbl($partFile);

my %partTbl; ## cluster => elements of the cluster
map { push @{$partTbl{$part{$_}}}, $_ } keys %part;

my %partFreq; ## number of elements per phylo partition cluster
map { $partFreq{$_}++ } values %part;

my $nPhyloParts = keys %partFreq;

if ($debug)
{
  print "\n\nNumber of phyloPart clusters: $nPhyloParts\n";

  print "\nPhylo partition sizes\n";
  printFormatedTbl(\%partFreq);
  print "\n";

  print "\nphyloPart part table\n";
  printFormatedTbl(\%part);

  print "\nphyloPart clusters:\n";
  for my $cl (keys %partTbl)
  {
    print "$cl:\n";
    for ( @{$partTbl{$cl}} )
    {
      print "\t$_\n";
    }
  }
  print "\n";
}

my $partCltrFile = "phyloPart_$taxon" . "_$percThld" . ".cltr";
print "--- Writing phylo partitioning to $partCltrFile\n";
writeTbl(\%part, $partCltrFile);

if ( $nPhyloParts > 1 )
{
  print "--- Generating annotation and query files\n";
  my $annFile    = "phyloPart_$taxon" . "_ann.tx";
  my $queryFile  = "phyloPart_$taxon" . "_query.seqIDs";
  my $vicutDir   = "phyloPart_$taxon" . "_vicut_dir";
  my $nQuerySeqs = 0;
  my $nAnnSeqs   = 0;
  my @queryTx;
  my @annTx;
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
      push @annTx, $id;
    }
  }
  close ANNOUT;
  close QOUT;

  print "--- Parsing tree leaves\n";
  my $treeLeavesFile = "phyloPart_$taxon" . "_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my @treeLeaves = readArray($treeLeavesFile);

  print "\n\n\tNumber of annotation seq's: $nAnnSeqs\n";
  print     "\tNumber of query seq's:      $nQuerySeqs\n";
  print     "\tSum:                        " . ($nAnnSeqs + $nQuerySeqs) . "\n";
  print     "\tNumber of leaves: " . @treeLeaves . "\n";

  my %part2;  # leaf label => cltrID
  my %cltr;   # cltrID => ref to array of IDs in the given cluster
  if ($nAnnSeqs > 0)
  {
    my $nParent = keys %parent;
    print     "\tNumber of elemets of the parent table: $nParent\n\n";

    if ( $nParent != @treeLeaves )
    {
      warn "\n\n\tERROR: $treeFile and $parentFile should have the same number of elements";
      print "\n\n";
      exit;
    }

    # Checking if keys of parent and leaves of the tree are the same sets
    my @parentElts = keys %parent;
    my @commElts = comm(\@parentElts, \@treeLeaves);
    if (@commElts != @parentElts || @commElts != @treeLeaves)
    {
      warn "\n\n\tERROR: $treeFile and $parentFile should have the same elements";
      print "Number of elements in $treeFile: " . @treeLeaves . "\n";
      print "Number of elements in $parentFile: $nParent\n";
      print "Number of common elements: " . @commElts . "\n\n";
      exit;
    }

    printArray(\@annTx, "Annotation taxons:");
    print "\n";

    printArray(\@queryTx, "Query taxons:");
    print "\n";

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
  }
  else
  {
    warn "\n\n\tERROR: annotation sequences missing";
    print "\n\n";
    exit;
    # # putting each leaf to its own cluster
    # my $count = 1;
    # %part2 = map { $_ => $count++ } keys %part;

    # for my $id (keys %part2)
    # {
    #   my $cl = $part2{$id};
    #   push @{$cltr{$cl}}, $id;
    # }
  }

  print "\nPhylo partition sizes\n";
  printFormatedTbl(\%partFreq);

  print "After vicut phylo partition sizes\n";
  my %part2Freq; ## number of elements per phylo partition cluster
  map { $part2Freq{$_}++ } values %part2;

  printFormatedTbl(\%part2Freq);

  my $nVicutCltrs = keys %part2Freq;

  if ( $nVicutCltrs <= $nPhyloParts )
  {
    print "--- Number of vicut clusters, $nVicutCltrs, is not greater than the number of phylo parts, $nPhyloParts\n";
    print "    Using phyloParts partition\n";

    my $clCounter = 1;
    undef %part2;
    undef %cltr;
    for my $id (keys %part)
    {
      if ($part{$id} != 0)
      {
	$part2{$id} = $part{$id};
      }
      else
      {
	$part2{$id} = $part{$id} . "_$clCounter";
	$clCounter++;
      }
      my $cl = $part2{$id};
      push @{$cltr{$cl}}, $id;
    }
  }

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
      $p =~ s/_\d+//;
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
  my %cltr2;
  my %txSubParent;
  my %parentIdx;
  my %txSubParentIdx; # this is for species_idx tree only
  @q = sort { @{$cltr{$b}} <=> @{$cltr{$a}} } keys %cltr;
  my $count = 1;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  for my $cl (@q)
  {
    my $tx = $clParent{$cl};
    $parentIdx{$tx}++;
    if ( $pars{$tx}>1 )
    {
      $tx = "$tx" . "_$parentIdx{$tx}";
    }
    else
    {
      $tx = "$tx";
    }
    $cltr2{$tx} = $cltr{$cl};
    $parentIdx{$tx} = $count;
    print "\tProcessing cluster $cl with new taxonomy $tx\n";
    for (@{$cltr{$cl}})
    {
      $txSubParent{$_} = $tx;
      print OUT "$_\t$tx\n";
      $txSubParentIdx{$_} = $count;
    }
    $count++;
  }
  close OUT;

  print "\nVicut updated phylo partition clusters with new names:\n";
  @q = sort { @{$cltr2{$b}} <=> @{$cltr2{$a}} } keys %cltr2;
  for my $cl (@q)
  {
    print "Cluster $cl => $parentIdx{$cl} (" . @{$cltr2{$cl}} . "):\n";
    for (@{$cltr2{$cl}})
    {
      print "\t$_\n";
    }
  }
  print "\n\n";

  print "\nVicut updated phylo partition cluster sizes:\n";
  @q = sort { @{$cltr2{$b}} <=> @{$cltr2{$a}} } keys %cltr2;
  for my $cl (@q)
  {
    print "Cluster $cl " . @{$cltr2{$cl}} . "\n";
  }
  print "\n\n";

  print "--- Creating tree with taxon_cluster leaf names\n";
  my $spClFile = "phyloPart_$taxon" . ".sppCl";
  my $spClFile2 = abs_path( "phyloPart_$taxon" . ".sppCl2" );
  open OUT, ">$spClFile" or die "Cannot open $spClFile for writing: $OS_ERROR\n";
  open OUT2, ">$spClFile2" or die "Cannot open $spClFile2 for writing: $OS_ERROR\n";
  for (keys %txSubParentIdx)
  {
    print OUT "$_\t$_" . "_cl_" . $txSubParentIdx{$_} . "\n";
    print OUT2 "$_" . "_cl_" . $txSubParentIdx{$_} . "\t" . $txSubParentIdx{$_} . "\n";
  }
  close OUT;
  close OUT2;
  print "\n\n--> Cluster index tbl written to $spClFile2\n" if $debug;

  my $treeFile2 = "phyloPart_$taxon" . "_final_condensed_cltrs.tree";
  $cmd = "nw_rename $treeFile $spClFile | nw_order -c n  - > $treeFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $nLeaves = scalar( keys %txSubParentIdx );

  ## calculating pdf file height
  ## the following formula is a linear model for the height
  ## so that when the number of leaves is 50 the height is 12in
  ## and when its 100, the height is 18in.

  ## from http://stackoverflow.com/questions/18532026/how-to-append-system-date-to-a-filename-in-perl
  my @now = localtime();
  my $timeStamp = sprintf("%04d-%02d-%02d_%02d_%02d_%02d",
			  $now[5]+1900, $now[4]+1, $now[3],
			  $now[2],      $now[1],   $now[0]);

  my $pdfTreeFile = abs_path( "phyloPart_$taxon" . "_cltrs_condensed_tree_$timeStamp.pdf" );
  my $treeFile2AbsPath = abs_path( $treeFile2 );

  plot_tree($treeFile2AbsPath, $spClFile2, $pdfTreeFile);

  if ( $showTree && $OSNAME eq "darwin")
  {
    $cmd = "open $pdfTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }
}
else
{
  ## $nPhyloParts < $nParents
  ## using parent table

  print "--- Phylo partition generated only one cluster\n";
  print "    Putting all taxons in one cluster\n";

  my @pars = sort { $parentFreq{$b} <=> $parentFreq{$a} } keys %parentFreq;
  my $tx = join "_", @pars;

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  for ( keys %parent )
  {
    print OUT "$_\t$tx\n";
  }
  close OUT;
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

sub plot_tree
{
  my ($treeFile, $clFile, $pdfFile, $title) = @_;

  my $readNewickFile = "/Users/pgajer/.Rlocal/read.newick.R";
  my $showBoostrapVals = "F";

  if (!defined $title)
  {
    $title = "";
  }

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

(nLeaves <- length(tr1\$tip.label))

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr1,type=\"phylogram\", no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8, tip.color=tip.colors, main=\"$title\")
par(op)
dev.off()
~;

  runRscript( $Rscript );
}

sub plotTree
{
  my ($treeFile, $clFile, $pdfFile, $figH) = @_;

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

pdf(\"$pdfFile\", width=6, height=$figH)
op <- par(mar=c(0,0,0,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr1,type=\"phylogram\", tip.color=tip.colors, no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8)
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

exit;
