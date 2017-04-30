#!/usr/bin/env perl

=head1 NAME

  mm_report_uq_spp.pl

=head1 DESCRIPTION

 Report unique species and their phylo-groups of PECAN taxonomic assignment of
 the M&M sequences

=head1 SYNOPSIS

  mm_report_uq_spp.pl

=head1 OPTIONS

=over

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  mm_report_uq_spp.pl

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "help|h!"                  => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
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

my $mmDir = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/";

## out file
##my $uqSppReport = $mmDir . "mm_uq_spp_report.txt";
my $uqSppReport = $mmDir . "mm_uq_spp_report_no_controls.txt";
open my $ROUT, ">$uqSppReport" or die "Cannot open $uqSppReport for writing: $OS_ERROR";

print "\r--- Parsing file with PECAN generated taxonomy on M&M's sequences\n";
##my $qTxFile = "/Users/pgajer/devel/MCextras/data/mm_pecan_tx2.txt";
my $qTxFile = "/Users/pgajer/devel/MCextras/data/mm_pecan_tx2_no_controls.txt"; # using only non-control seq's

## $qTxFile file content: <seqID> <tx classification of seqID> <pp> <its phylo-group>

# $ head /Users/pgajer/devel/MCextras/data/mm_pecan_tx2.txt

#                                                tx   pp                        phGr
# 1642.V1_0                     Lactobacillus_iners 0.93     Firmicutes_group_6_V3V4
# 0980.V2_1                     Lactobacillus_iners 0.97     Firmicutes_group_6_V3V4
# 1670.V2_2 Lactobacillus_crispatus_kefiranofaciens 0.98     Firmicutes_group_6_V3V4
# 0711.V3_3                       Atopobium_vaginae 0.56 Actinobacteria_group_0_V3V4
# 1149.V1_4                     Lactobacillus_iners 0.94     Firmicutes_group_6_V3V4
# 1386.V1_5                     Prevotella_buccalis 0.85  Bacteroidetes_group_2_V3V4

my ($rspIdsTbl, $rphGrSppTbl, $rppTbl) = readQtxTbl($qTxFile);

my %spIdsTbl   = %{$rspIdsTbl};   # sp   => ref of array with seqIDs of seq's classified to sp
my %phGrSppTbl = %{$rphGrSppTbl}; # phGr => species of phGr
my %ppTbl      = %{$rppTbl};      # seqID => posterior probability of the best model

for my $phGr ( keys %phGrSppTbl )
{
  my @uqSpp = unique($phGrSppTbl{$phGr});
  my @spp = sort { @{$spIdsTbl{$b}} <=> @{$spIdsTbl{$a}} } @uqSpp;
  for my $spIdx ( 0..$#spp )
  {
    my $sp = $spp[$spIdx];
    my @ids = @{$spIdsTbl{$sp}}; # seq IDs of $sp
    my $nSp = scalar(@ids); # number of sequences of the given species

    print "Processing $sp (n=" . commify($nSp) . ")\n";

    print $ROUT "$sp\t$nSp\t$phGr\n";

  } ## end of    for my $spIdx (0..
} ## end of   for my $phGr ( keys %phGrSppTbl )

close $ROUT;

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

print "\n\n\tReport written to $uqSppReport\n\n";

####################################################################
##                               SUBS
####################################################################

# parse a clstr2 file
# output table: refId -> number of elements in the corresponding cluster
sub parseClstr2
{
  my $inFile = shift;

  my %tbl;
  open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
  foreach my $rec (<IN>)
  {
    chomp $rec;
    my @ids = split ",", $rec;
    my $refId = shift @ids;
    ##$tbl{$refId} = \@ids;
    $tbl{$refId} = @ids; # we are only interested in the size of the cluseter
  }
  close IN;

  return %tbl;
}

##
## read 3 column taxon table
##

# file format

#                                                tx   pp                        phGr
# 1642.V1_0                     Lactobacillus_iners 0.93     Firmicutes_group_6_V3V4
# 0980.V2_1                     Lactobacillus_iners 0.97     Firmicutes_group_6_V3V4
# 1670.V2_2 Lactobacillus_crispatus_kefiranofaciens 0.98     Firmicutes_group_6_V3V4
# 0711.V3_3                       Atopobium_vaginae 0.56 Actinobacteria_group_0_V3V4
# 1149.V1_4                     Lactobacillus_iners 0.94     Firmicutes_group_6_V3V4
# 1386.V1_5                     Prevotella_buccalis 0.85  Bacteroidetes_group_2_V3V4
# ...
sub readQtxTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readQtxTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %spIdsTbl;
  my %phGrSppTbl;
  my %ppTbl;
  ##my %sp2phGrSppTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    next if /^$/;
    chomp;
    my ($id, $sp, $pp, $gr) = split /\s+/,$_;
    push @{$spIdsTbl{$sp}}, $id;
    push @{$phGrSppTbl{$gr}}, $sp;
    $ppTbl{$id} = $pp;
    ##$sp2phGrSppTbl{$sp} = $gr;
  }
  close IN;

  return (\%spIdsTbl, \%phGrSppTbl, \%ppTbl);
}

# read 3 column clstrs table
sub readCltrsTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %vCltrTbl;
  my %txTbl;
  my %txTbl2;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    if (defined $id)
    {
      $vCltrTbl{$id} = "c" . $cl;
      $txTbl{$id} = $tx;
      if ($tx ne "NA")
      {
	$txTbl2{$id} = $tx;
      }
      else
      {
	$txTbl2{$id} = "c" . $cl;
      }
    }
  }
  close IN;

  return (\%vCltrTbl, \%txTbl, \%txTbl2);
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
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t, $r) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

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
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
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

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR";
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

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
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
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  #print "\n";
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTblToFile{

  my ($rTbl, $rSub, $fh) = @_; # the second argument is a subarray of the keys of the table

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

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify {
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}

# parse a CSV partition table
sub read_part_tbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit 1;
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

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

## Testing if two arrays are identical in a set-theoretic sense. That is that
## they have exactly the same set of elements.
sub setEqual
{
  my ($rA, $rB) = @_;

  my @a = @{$rA};
  my @b = @{$rB};
  my @c = comm(\@a, \@b);

  my $ret = 1;

  if (@c != @a || @c != @b)
  {
    warn "\n\n\tERROR: Elements of the two arrays do not match";
    print "\n\tNumber of elements in the first array: " . @a . "\n";
    print "\tNumber of elements in the second array: " . @b . "\n";
    print "\tNumber of common elements: " . @c . "\n";

    # writeArray(\@a, "a.txt");
    # writeArray(\@b, "b.txt");
    #print "\n\tNew taxon keys and fasta IDs written to a.txt and b.txt, respectively\n\n";

    if (@a > @b)
    {
      my @d = diff(\@a, \@b);
      print "\nElements a but not b:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    if (@b > @a)
    {
      my @d = diff(\@b, \@a);
      print "\nElements in b that are not a:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    $ret = 0;
  }

  return $ret;
}

exit 0;
