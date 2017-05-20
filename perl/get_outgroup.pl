#!/usr/bin/env perl

=head1 NAME

  get_outgroup.pl

=head1 DESCRIPTION

  Given a full taxonomic lineage file, a taxonomic rank (say phylum), and the maximum
  desired number of sequences for each group at that rank, generate a list of outgroup
  sequence IDs for each group containing fewer than the desired maximum number of sequences.
  An outgroup sequence ID is randomly chosen from each sister group at that rank and placed
  into the generated *_outgroup.seqID file.

=head1 SYNOPSIS

  get_outgroup.pl -i <full tx file> -t <taxonomic rank> -s <group size> [Options]

=head1 OPTIONS

=over

=item B<--in-tx-file, -i>
  Input lineage file.

=item B<--tx-rank, -t>
  Taxonomic rank. Possible values: domain, phylum, class, order, family, genus, species.

=item B<--group-size, -s>
  Maximum size of group.

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

  cd /local/scratch/jbh_pecan_new_tx/
  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir

  get_outgroup.pl -i ../rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.lineage -t phylum -s 10000 -o rdp_Bacteria_phylum_OGs_dir

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper qw(Dumper);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $minLim = 1000;

GetOptions(
  "input-file|i=s"  => \my $inFile,
  "output-dir|o=s"  => \my $outDir,
  "group-size|s=i"  => \my $grpSz,
  "min-limit|m=i"   => \$minLim,
  "tx-rank|t=s"     => \my $txRank,
  "verbatim|v"      => \my $verbatim,
  "debug"           => \my $debug,
  "dry-run"         => \my $dryRun,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$inFile)
{
  print "\n\nERROR: Missing input file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$txRank)
{
  print "\n\nERROR: Missing taxonomic rank\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outDir)
{
  print "\n\nERROR: Missing output directory\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my %txRankTbl;
my @txRanks = ("domain", "phylum", "class", "order", "family", "genus", "species");
@txRankTbl{@txRanks} = 1..7;

my @txPrefix = ("root_", "d_", "p_", "c_", "o_", "f_", "g_", "s_");

if ( !exists $txRankTbl{$txRank} )
{
  print "\n\nERROR: Unrecognized taxonomic rank: $txRank\n\n\n";
  exit;
}

if ( ! -f $inFile )
{
  print "\n\nERROR: input file does not exist\n\n\n";
  exit;
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

my $masterLineageFile = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.lineage";
my $masterFaFile      = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.fa";
my $masterSeqLenFile  = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.seqLen";

print "--- Parsing master lineage table\n";
my %mLineageTbl = read2colTbl($masterLineageFile);

print "--- Parsing master seqLen table\n";
my %mSeqLen = read2colTbl($masterSeqLenFile);

print "--- Generating master children/parent tables\n";
my ($rmParent, $rmChildren, $rid2sp) = get_parent_children_tbls(\%mLineageTbl);

my %id2sp = %{$rid2sp};


my $perc = 0;
my $fileSize = -s $inFile; # file size in bytes

my $idx = $txRankTbl{$txRank}; # family
my $pidx = $idx - 1;           # order
my $ppidx = $idx - 2;          # class
my $p3idx = $idx - 3;          # phylum
my $p4idx = $idx - 4;          # domain

print "\nidx: $idx\n\n$ppidx\n" if $debug;
#exit;

print "--- Deciphering parent/child relationships from the full taxonomic lineage file\n";
open (IN, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
my $count = 0;
my %children;
my %parent;
my %freq;
my %ids;
while (<IN>)
{
  if ($count % 10000 == 0)
  {
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

    my $perc = sprintf("%.1f%%", 100 * (tell IN) / $fileSize);
    print "\r$timeStr [$perc]";
  }
  $count++;

  chomp;
  my ($id, $t) = split /\s+/,$_; ## split the input taxonomy at the space into seqID and lineage
  my @f = split ";", $t; ## place the lineage in an array split the lineage into rank

  if ($debug)
  {
    print "\nline: $_\n";
    print "id: $id\tt: $t\n";
    print "f[$idx]: $f[$idx]\n";
    print "f[$pidx]: $f[$pidx]\n";
    ##exit;
  }

  push @{$ids{$f[$idx]}}, $id;
  push @{$ids{$f[$pidx]}}, $id;

  $children{$f[$pidx]}{$f[$idx]}++;
  $parent{$f[$idx]} = $f[$pidx];

  if ($ppidx > -1)
  {
    $parent{$f[$pidx]} = $f[$ppidx];
    $children{$f[$ppidx]}{$f[$pidx]}++;
    push @{$ids{$f[$ppidx]}}, $id;

    if ($p3idx > -1)
    {
      $parent{$f[$ppidx]} = $f[$p3idx];
      $children{$f[$p3idx]}{$f[$ppidx]}++;
      push @{$ids{$f[$p3idx]}}, $id;

      if ($p4idx > -1)
      {
	$parent{$f[$p3idx]} = $f[$p4idx];
	$children{$f[$p4idx]}{$f[$p3idx]}++;
      }
    }
  }
  $freq{$f[$idx]}++;
}
close IN;

if ($debug)
{
	print Dumper \%freq;
	print Dumper \%parent;
	print Dumper \%children;
	#exit;
}
print "\r--- Finding outgroup sequences for groups with 1,000-$grpSz sequences.\n";

my $cmd = "mkdir -p $outDir";
print "\tcmd=$cmd\n" if $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

chdir $outDir;

my @largerPhyla;
my @smallerPhyla;
my @sPhyla;
my $outFile;
my $outFile2;
my $nMidRks = 0;
for my $rk (keys %freq ) ## For each group at the specified rank - use phylum,
{
  print "\n rk: $rk\tfreq: $freq{$rk}\n";
  if ( $freq{$rk} <= $grpSz && $freq{$rk} > $minLim )
  {
    $nMidRks++;
    $outFile = $txPrefix[$idx].$rk."_outgroup.seqID";
    $outFile2 = $txPrefix[$idx].$rk."_outgroup.tx";
    open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
    open OUT2, ">$outFile2" or die "Cannot open $outFile2 for writing: $OS_ERROR\n";
    my $p;
    if ( exists $parent{$rk}) ## if a parent exists for the group (say group = phylum)
    {
      $p = $parent{$rk}; ## set p to the parent
    }
    else
    {
      print "\n\nERROR: parent for $rk not defined\n\n\n";
      exit;
    }

    my @ch;
    if ( exists $children{$p} ) ## if children exist for p
    {
      @ch = keys %{$children{$p}}; ## place all of the children in array ch
    }
    else
    {
      print "\n\nERROR: children for $p not defined\n\n\n";
      exit;
    }

    if ($debug)
    {
      print "\nrk: $rk\tp: $p\tch: " .  @ch . "\n"; ## print out the group, its parent and the number of children the parent has
      printArray(\@ch, "ch");
      print "\n";
    }

    if ( @ch > 1 ) ## If there is greater than one child of the parent rank,
    {
      my @sib;
      for my $c (@ch)  ## For each child
      {
	if ($c ne $rk)	## If the child is not equal to rk
	{
	  push @sib, $c;	## Push that child to the sibling array
	}
      }

      foreach my $sib (@sib)
      {
	my @a = @{$ids{$sib}}; ## get all ids equal to the sibling
	my $id = $a[rand @a];  ## pick a random id from the sibling id array NOTE: this biases species with lots of sequences
	print OUT "$id\n";
	print OUT2 "$id\t" . $id2sp{$id} . "\n";
      }
      print "--- Outgroup sequence IDs for $txRank $rk written to $outFile\n";
    }
    close OUT;
    close OUT2;
  }
  elsif ($freq{$rk} <= $minLim)
  {
    push @smallerPhyla, "$rk\t".$freq{$rk}." sequences\n";
    push @sPhyla, $rk;
  }
  else
  {
    push @largerPhyla,  "$rk\t".$freq{$rk}." sequences\n";
  }
}

## getting outgroup seq's for @smallerPhyla

my $p;
my $rk = $sPhyla[0];
if ( exists $parent{$rk}) ## if a parent exists for the group (say group = phylum)
{
  $p = $parent{$rk}; ## set p to the parent
}
else
{
  print "\n\nERROR: parent for $rk not defined\n\n\n";
  exit;
}

my @ch;
if ( exists $children{$p} ) ## if children exist for p
{
  @ch = keys %{$children{$p}}; ## place all of the children in array ch
}
else
{
  print "\n\nERROR: children for $p not defined\n\n\n";
  exit;
}

my @sibs = diff(\@ch, \@sPhyla);
$outFile = "phyla_lessthen_1k_outgroup.seqIDs";
$outFile2 = "phyla_lessthen_1k_outgroup.tx";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
open OUT2, ">$outFile2" or die "Cannot open $outFile2 for writing: $OS_ERROR\n";
foreach my $sib (@sibs)
{
  my @a = @{$ids{$sib}}; ## get all ids equal to the sibling
  my $id = $a[rand @a]; ## pick a random id from the sibling id array
  print OUT "$id\n"; ## Push the seqID.
  print OUT2 "$id\t". $id2sp{$id} . "\n";
}
close OUT;
close OUT2;
print "--- Outgroup sequence IDs for phyla_lessthen_1k written to $outFile\n";


my $nRk = keys %freq;
my $larger = "phyla_greater_than_$grpSz.txt";
open LARGER, ">$larger" or die "Cannot open $larger for writing: $OS_ERROR\n";
print LARGER @largerPhyla;
close LARGER;

my $smaller = "phyla_less_than_1k.txt";
open SMALLER, ">$smaller" or die "Cannot open $smaller for writing: $OS_ERROR\n";
print SMALLER @smallerPhyla;
close SMALLER;

print "\n\n\tNumber of $txRank members:                            $nRk\n";
print     "\tNumber of $txRank members of size > 1k and <= $grpSz:  $nMidRks\n";
print     "\tNumber of $txRank members of size > $grpSz:            " . @largerPhyla ." \n";
print     "\tNumber of $txRank members of size <= 1k:              " . @smallerPhyla ." \n\n";

#print "The following $txRank contained </= 1000 sequences:\n $complete\n\n";
print "Groups containing > $grpSz sequences and are listed in $larger.\n";
print "Groups containing < 1,000 sequences and are listed in $smaller.\n\n";

## report timing
$endRun = time();
$runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  $timeMin = int($runTime / 60);
  $timeSec = sprintf("%02d", $runTime % 60);
  print "\rCompleted in $timeMin:$timeSec\n"
}
else
{
  print "\rCompleted in $runTime seconds\n"
}

print "\n";

####################################################################
##                               SUBS
####################################################################

## given a taxon name of the form f_Lactobacillaceae
## extract the taxon prefix, in this example 'f'
sub get_tx_prefix
{
  my $tx = shift;
  my @f = split "_", $tx;
  my $prefix = shift @f;
  return $prefix;
}

sub get_parent_children_tbls
{
  my $rLineage = shift;

  my %lineageTbl = %{$rLineage};

  my %id2sp;

  # my %spTbl;
  # my %geTbl;
  # my %faTbl;
  # my %orTbl;
  # my %clTbl;
  # my %phTbl;

  my %children;
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

    $id2sp{$id} = $sp;

    ##$sp = "s_$sp";
    ##$sp .= "_OG" if ( exists $ogInd{$id} );
    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $parent{$sp} = $ge;
    $parent{$ge} = $fa;
    $parent{$fa} = $or;
    $parent{$or} = $cl;
    $parent{$cl} = $ph;
    $parent{$ph} = "d_Bacteria";

    $children{"d_Bacteria"}{$ph}++;
    $children{$ph}{$cl}++;
    $children{$cl}{$or}++;
    $children{$or}{$fa}++;
    $children{$fa}{$ge}++;
    $children{$ge}{$sp}++;

    # $spTbl{$sp}++;
    # $geTbl{$ge}++;
    # $faTbl{$fa}++;
    # $orTbl{$or}++;
    # $clTbl{$cl}++;
    # $phTbl{$ph}++;
  }

  return (\%parent, \%children, \%id2sp);
}

## get random sequence from ONE randomly selected species of the given taxon
sub get_one_rand_sp_seq
{
  my ($tx, $rSpTbl, $rSeqLen, $perc) = @_;

  my $prop = $perc / 100.0;
  my %seqLen = %{$rSeqLen};
  my %spTbl = %{$rSpTbl};

  my $c = get_tx_prefix($tx);
  my $id;
  if ($c eq 's')
  {
    $id = get_rand_seqID_from_sp($tx, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'g')
  {
    my @sps = keys %{$children{$tx}};
    # pick randomly a species and then select from it
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'f')
  {
    my @ges = keys %{$children{$tx}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'o')
  {
    my @fas = keys %{$children{$tx}};
    my $fa = $fas[rand @fas];
    my @ges = keys %{$children{$fa}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'c')
  {
    my @ors = keys %{$children{$tx}};
    my $or = $ors[rand @ors];
    my @fas = keys %{$children{$or}};
    my $fa = $fas[rand @fas];
    my @ges = keys %{$children{$fa}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'p')
  {
    my @cls = keys %{$children{$tx}};
    my $cl = $cls[rand @cls];
    my @ors = keys %{$children{$cl}};
    my $or = $ors[rand @ors];
    my @fas = keys %{$children{$or}};
    my $fa = $fas[rand @fas];
    my @ges = keys %{$children{$fa}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    ##$id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
    $id = get_rand_seqID_from_sp($sp, $perc);
  }

  return $id;
}

## get random sequence ID from seq IDs of a given species here the selection is using master lineage of all bacteria
## the selected seq should have length greater than in the perc of the max length
sub get_rand_seqID_from_sp
{
  my ($sp, $rSpTbl, $perc) = @_;

  my %spTbl = %{$rSpTbl};

  my $c = get_tx_prefix($sp);
  if ($c ne 's')
  {
    print "ERROR in get_rand_seqID_from_sp(): input taxon $sp is not a species\n";
    exit;
  }

  if (!defined $spTbl{$sp})
  {
    print "ERROR in get_rand_seqID_from_sp(): Could not find $sp in masterSpTbl\n";
    exit;
  }

  my @s = @{$spTbl{$sp}};
  my @len = @mSeqLen{@s};
  my $imax = argmax(\@len);
  my $maxLen = $len[$imax];

  my $prop = $perc / 100.0;
  my $goodLen = $prop * $maxLen;
  my @ss = grep{ $mSeqLen{$_} > $goodLen } @s;
  my $id = $ss[rand @ss];

  return $id;
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

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "$header\n" if $header;
  map {print "\t$_\n"} @{$a};
}

# difference of two arrays
sub diff
{
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

exit;
