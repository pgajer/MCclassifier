#!/usr/bin/env perl

=head1 NAME

  rm_seqs_from_phGr.pl

=head1 DESCRIPTION

  This script removes a set of ref sequences from a given phylogroup updating the
  corresponding taxon, lineage, outgropu and alignment. The original taxon,
  lineage, outgropu and alignment files are moved to a backup dir.

=head1 SYNOPSIS

  rm_seqs_from_phGr.pl -e <seq file> -i <group prefix> [Options]

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Phylo group prefix.

=item B<--seq-file, -e>
  File with sequences to be removed.

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

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Proteobacteria_dir/tmp1

  rm_seqs_from_phGr.pl --debug -e Proteobacteria_group_10_V3V4_dir/bad.seqIDs -i Proteobacteria_group_10_V3V4

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
  "group-prefix|i=s" => \my $grPrefix,
  "seq-file|e=s"     => \my $seqFile,
  "keep-tmp-files"   => \my $keepTmpFiles,
  "quiet"            => \my $quiet,
  "igs"              => \my $igs,
  "verbose|v"        => \my $verbose,
  "debug"            => \my $debug,
  "dry-run"          => \my $dryRun,
  "help|h!"          => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$grPrefix)
{
  print "\n\nERROR: Missing group prefix\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$seqFile)
{
  print "\n\nERROR: Missing seq file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $debugStr = "";
my $quietStr = "--quiet";
if ($debug)
{
  $debugStr  = "--debug";
  $quietStr  = "";
}

my $verboseStr = "";
if ($verbose)
{
  $verboseStr  = "--verbose";
}

my $FastTree    = "FastTree";
my $select_seqs = "select_seqs.pl";
my $select_tx   = "select_tx.pl";

if ( defined $igs )
{
  $FastTree     = "/home/pgajer/bin/FastTree_no_openMP";
  $select_seqs  = "/home/pgajer/devel/MCclassifier/perl/select_seqs.pl";
  $select_tx    = "/home/pgajer/devel/MCclassifier/perl/select_tx.pl";
}

####################################################################
##                               MAIN
####################################################################

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  exit 1;
}

print "--- Parsing file with seq's to be removed\n" if !$quiet;
my @selSeqs = read_array( $seqFile );
#my %seqTbl = map { $_ => 1 } @selSeqs;


print "--- Changing dir to $grDir\n";

chdir $grDir;
##$grPrefix = "$grDir/$grPrefix";

$seqFile =~ s/$grDir\///;

my $liFile    = $grPrefix . ".lineage";
my $txFile    = $grPrefix . ".tx";
my $treeFile  = $grPrefix . ".tree";
my $algnFile  = $grPrefix . "_ginsi_algn.fa";
my $ogFile    = $grPrefix . "_outgroup.seqIDs";


print "--- Creating directory for storing the original files\n" if !$quiet;
my @now = localtime();
my $timeStamp = sprintf("%04d%02d%02d_%02d%02d%02d",
			$now[5]+1900, $now[4]+1, $now[3],
			$now[2],      $now[1],   $now[0]);

my $bkpDir = "orig_data_dir_$timeStamp";
my $cmd = "mkdir $bkpDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Copying the original files to the backup dir\n" if !$quiet;

$cmd = "cp $liFile $bkpDir/$liFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "cp $txFile $bkpDir/$txFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "cp $algnFile $bkpDir/$algnFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "cp $ogFile $bkpDir/$ogFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "cp $treeFile $bkpDir/$treeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Removing given seq's from the lineage file\n" if !$quiet;

my $liCountBefore = line_count( "$bkpDir/$liFile" );
$cmd = "$select_tx $quietStr -e $seqFile -i $bkpDir/$liFile -o $liFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
my $liCountAfter = line_count( $liFile );

if ( !$quiet )
{
  print "\nLineage file line count BEFORE: $liCountBefore\n";
  print   "Lineage file line count AFTER:  $liCountAfter\n\n";
}

print "--- Removing given seq's from the taxonomy file\n" if !$quiet;
my $txCountBefore = line_count( "$bkpDir/$txFile" );
$cmd = "$select_tx $quietStr -e $seqFile -i $bkpDir/$txFile -o $txFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
my $txCountAfter = line_count( $txFile );

if ( !$quiet )
{
  print "\nTaxonomy file line count BEFORE: $txCountBefore\n";
  print   "Taxonomy file line count AFTER:  $txCountAfter\n\n";
}


print "--- Removing given seq's from the alignment file\n" if !$quiet;
my $seqCountBefore = seq_count( "$bkpDir/$algnFile" );
$cmd = "$select_seqs $quietStr -e $seqFile -i $bkpDir/$algnFile -o $algnFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
my $seqCountAfter = seq_count( $algnFile );

if ( !$quiet )
{
  print "\nAlignment file line count BEFORE: $seqCountBefore\n";
  print   "Alignment file line count AFTER:  $seqCountAfter\n\n";
}


print "--- Removing given seq's from the OG file\n" if !$quiet;
my @og = read_array( "$bkpDir/$ogFile" );
my @d  = diff( \@og, \@selSeqs );

if ( @d )
{
  write_array( \@d, $ogFile );
}
else
{
  warn "\n\n\tWARNING: All outgroup seq's were removed\n";
  print "\n\n";
}

if ( $seqCountBefore != $seqCountAfter )
{
  print "--- Rebuilding tree using pruned alignment\n" if !$quiet;
  $cmd = "rm -f $treeFile; $FastTree -nt $algnFile > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
else
{
  print "\n\n\tWARNING: Do not rebuilding tree as the alignment file was not modified\n\n";
}


####################################################################
##                               SUBS
####################################################################

# extract line count of a file
sub line_count
{
  my $file = shift;

  my $wcline = qx/ wc -l $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
}

# extract sequence count of a fasta file
sub seq_count
{
  my $file = shift;

  my $wcline = qx/ grep -c '>' $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
}

# read table with one column
sub read_array
{
  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -e $file )
  {
    warn "\n\n\nERROR: $file does not exist";
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

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;

  local $" = ', ';
  ##local $, = ',';
  print "$header: " if $header;
  print "@{$a}\n";
}

# print array to stdout
sub print_arrayByRow
{
  my ($a, $header) = @_;

  local $" = '\n ';
  ##local $, = ',';
  print "$header:\n" if $header;
  map { print "$_\n" } @{$a};
  print "\n\n";
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify
{
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}

## print two strings (key and scalar value of a hash table) formated using $maxKeyLen, $maxValLen generated using previeous routine
## (4, $fa, scalar(@{$faTbl{$fa}}), $faKVL);
sub printF
{
  my ($nTabs, $key, $val, $rKV ) = @_;

  my $tabs = "";
  my $i = 0;
  while ( $i < $nTabs)
  {
    $tabs .= "\t";
    $i++;
  }

  my $maxKeyLen = $rKV->[0];
  my $maxValLen = $rKV->[1];
  my $n = ($maxKeyLen - length($key)) + ($maxValLen - length($val));
  my $pad = ": ";
  for (my $i=0; $i<$n; $i++)
  {
    $pad .= " ";
  }
  print "$tabs$key$pad$val\n";
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
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

sub get_taxons
{
  my $rLineage = shift;

  my %lineageTbl = %{$rLineage};

  my @taxons;

  my %spTbl;
  my %geTbl;
  my %faTbl;
  my %orTbl;
  my %clTbl;
  my %phTbl;

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

    $spTbl{$sp}++;
    $geTbl{$ge}++;
    $faTbl{$fa}++;
    $orTbl{$or}++;
    $clTbl{$cl}++;
    $phTbl{$ph}++;
  }


  my @phs = keys %phTbl;
  for my $ph (@phs)
  {
    push @taxons, $ph;
    my @cls = keys %{$children{$ph}};
    for my $cl (@cls)
    {
      push @taxons, $cl;
      my @ors = keys %{$children{$cl}};
      for my $or (@ors)
      {
	push @taxons, $or;
	my @fas = keys %{$children{$or}};
	push @taxons, @fas;
      }
    }
  }

  return (\@taxons, \%parent, \%children);
}

# use a two-ranks-up rule for OG selection: find the highest taxonomic rank (say
# class), go to grandparent (the parent of the parent) (Bacteria in the case of
# the class), take for OG children of the grandparent and subtract from them our
# phylum. Those are taxonomic ranks from which the OGs will be selected. In the
# case of family, we start from genus and go two levels up, so we are in the
# order and so all other families of that order would be the OGs of the given
# family

# For Group 0 of Firmicutes - Find the highest taxonomic rank within Group 0 - it
# is class. Go to the grandparent - Bacteria. Take for OGs the children of the
# grandparent (phyla) that are not the phylum of Group 0.

# For Bacilli, which is separated at lower ranks, say family, outgroups would be
# chosen from sister orders excluding the order containing the family of interest.

sub get_OG_families
{
  my ($rTaxons, $rParent, $rChildren) = @_;

  my @taxons = @{$rTaxons};
  my %parent = %{$rParent};
  my %children = %{$rChildren};

  my %txInt;
  $txInt{'p'} = 0;
  $txInt{'c'} = 1;
  $txInt{'o'} = 2;
  $txInt{'f'} = 3;
  $txInt{'g'} = 4;
  $txInt{'s'} = 5;

  my %grTx; # group taxons
  my %pTbl;
  my @d;
  my @ch;
  my $imin = 5;
  my $cmin;

  my $debug = 1;

  for my $tx ( @taxons )
  {
    my $c = get_tx_prefix($tx);
    my $i = $txInt{$c};
    push @{$grTx{$c}}, $tx;
  }

  print_arrayTbl(\%grTx, "grTx") if $debug;

  ## pick the lowest taxon with the single element
  my @cc = ('s','g','f','o','c','p');
  for my $c (@cc)
  {
    if ( exists $grTx{$c} && @{$grTx{$c}}==1 )
    {
      $cmin = $c;
      last;
    }
  }
  my @t;
  my $gp; # parent
  my $p;  # parent
  if ( exists $grTx{$cmin} )
  {
    @t = @{$grTx{$cmin}}; # extracting highest level taxons
    $gp = $parent{$t[0]};
    $p = $t[0];
  }
  else
  {
    print "\nERROR: grTx{$cmin} does not exist\n";
    print "Group Taxons\n";
    for my $tx (@taxons)
    {
      my $c = get_tx_prefix($tx);
      my $i = $txInt{$c};
      print "tx: $tx\tc: $c\ti: $i\n";
      if ( $i > $imin )
      {
	$imin = $i;
	$cmin = $c;
      }
    }
    print "imin: $imin\n";
    print "cmin: $cmin\n";
    print_arrayTbl(\%grTx, "grTbl");
    print "\n";
    exit 1;
  }

  print_array(\@t, "t=grTx{cmin}") if $debug;

  # for (@t)
  # {
  #   $pTbl{$parent{$_}}++ if exists $parent{$_};
  # }
  # ##$p = $parent{$t[0]};# its enough to take the parent of the first taxonomic rank as all taxonomic ranks at the lowest level have to have the same parent

  # printTbl(\%pTbl, "pTbl") if $debug;

  # for my $p (keys %pTbl)
  # {
  #   push @ch, keys %{$children{$p}};
  # }

  @ch = keys %{$children{$gp}};

  print_array(\@ch, "ch") if $debug;

  @t = ($p);
  @d = diff(\@ch, \@t);

  print_array(\@d, "d") if $debug;

  if (@d == 0)
  {
    print "\nERROR: Difference between ch and t is zero\n";
    print "Taxons:\n";
    for my $tx (@taxons)
    {
      print "$tx\n";
    }
    print "\n";

    if (@taxons > 1)
    {
      printTbl(\%pTbl, "pTbl");
      print "\n";
    }
    else
    {
      my $p = $parent{$rTaxons->[0]};
      print "p: $p\n";
    }

    print "ch: @ch\n";
    print "d: @d\n";
    exit 1;
  }

  return @d;
}

# get families of the given taxonomic rank
sub get_families
{
  my ($tx, $rParent, $rChildren) = @_;

  my %parent = %{$rParent};
  my %children = %{$rChildren};

  my $c = get_tx_prefix($tx);

  my %txInt;
  $txInt{'p'} = 0;
  $txInt{'c'} = 1;
  $txInt{'o'} = 2;
  $txInt{'f'} = 3;
  $txInt{'g'} = 4;
  $txInt{'s'} = 5;

  my $i = $txInt{$c};
  if ($i>3)
  {
    print "ERROR: in get_families() at __LINE__ the input taxon $tx is genus or species\n";
    return "";
  }

  my @fas;

  if ($c eq 'f')
  {
    ##print "WARNING: in get_families() the input taxon is a family\n";
    push @fas, $tx;
  }
  elsif ($c eq 'o')
  {
    @fas = keys %{$children{$tx}};
  }
  elsif ($c eq 'c')
  {
    my @ors = keys %{$children{$tx}};
    for (@ors)
    {
      push @fas, keys %{$children{$_}};
    }
  }
  elsif ($c eq 'p')
  {
    my @cls = keys %{$children{$tx}};
    for my $cl (@cls)
    {
      my @ors = keys %{$children{$cl}};
      for my $or (@ors)
      {
	push @fas, keys %{$children{$_}};
      }
    }
  }

  return @fas;
}


## given a taxon name of the form f_Lactobacillaceae
## extract the taxon prefix, in this example 'f'
sub get_tx_prefix
{
  my $tx = shift;
  my @f = split "_", $tx;
  my $prefix = shift @f;
  return $prefix;
}

# print elements of a hash table
sub printTbl
{
  my ($rTbl, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
  print "\n";
}

# write array to a file (one column format)
sub write_array
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

exit 0;
