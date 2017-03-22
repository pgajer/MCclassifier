#! /usr/bin/perl

=head1 NAME

  cmp_tx.pl

=head1 DESCRIPTION

  Compare taxonomic assignments

=head1 SYNOPSIS

  cmp_tx.pl -i <file 1 with tx> -j <file 2 with tx> -o <output directory> [Options]

=head1 OPTIONS


=over

=item B<--input-file1, -i>
  Input file 1 with two columns seqID, taxonomy

=item B<--input-file2, -j>
  Input file 2 with two or three columns seqID, taxonomy, error_probability

=item B<--output-dir, -o>
  Output dir with the following files.

  comb.tx
     seqID tx1 tx2 matchIndicator
     where matchIndicator = 1 if assignment was the same in both files and 0 otherwise.

  spp.summary
  <file1 species name> <number of sequences for given species>
   file2 classification list for sequences of the species

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /Users/pgajer/devel/packages/vaginal-0.2.1.refpkg

  cmp_tx.pl -i vaginal_27F_534R_nr.tx -j vaginal_27F_534R_nr.MC.order3.otu -o .

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
  "input-file1|i=s"  => \my $inFile1,
  "input-file2|j=s"  => \my $inFile2,
  "output-dir|o=s"   => \my $outDir,
  "dry-run"          => \my $dryRun,
  "help|h!"          => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$inFile1)
{
  print "ERROR: Missing input file 1\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$inFile2)
{
  print "ERROR: Missing input file 2\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outDir)
{
  print "ERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################

# creating output directory
my $cmd = "mkdir -p $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

if ( -l $inFile1 )
{
  $inFile1 = readlink($inFile1);
}
if ( -l $inFile2 )
{
  $inFile2 = readlink($inFile1);
}

# reading file 1 and file 2 tables
my %tbl1 = read2colTbl($inFile1);
my %tbl2 = read2colTbl($inFile2);
#my ($rT, $rE) = read3colTbl($inFile2);
#my %tbl2 = %{$rT};

my %sppTbl; # for each species from file 1 assign frequencies of taxonomic
	    # assignments to the sequences of that species as given in file2

my $outFile1 = "$outDir/comb.tx";
open OUT, ">$outFile1" or die "Cannot open $outFile1 for writing: $OS_ERROR\n";
my $count = 0;
my $match = 0;
my @mismatchedSpp;
for my $id (keys %tbl1)
{
  next if !exists $tbl1{$id} || !exists $tbl2{$id};

  $count++;
  # print "$id";
  # print "\t" . $tbl1{$id};
  # print "\t" . $tbl2{$id} . "\n";
  my $b = int($tbl1{$id} eq $tbl2{$id});
  $match += $b;
  if ($b==0)
  {
    push @mismatchedSpp, $tbl1{$id};
  }
  print OUT "$id\t" . $tbl1{$id} . "\t" . $tbl2{$id} . "\n";
  $sppTbl{$tbl1{$id}}{$tbl2{$id}}++;
}
close OUT;


@mismatchedSpp = unique(\@mismatchedSpp);

print "\nNumber of taxons:\t" . commify($count) . "\n";

my $pMatch = sprintf("%.2f",100 * $match / $count);
##print "Percentage of matches: $pMatch%\n";
print "Matches:\t\t " . commify($match) . " ($pMatch\%)\n";

my $pMismatch = sprintf("%.2f",100 * ($count - $match) / $count);
##print "Percentage of mismatches: $pMismatch\n";
print "Mismatches:\t\t " . commify($count - $match) . " ($pMismatch\%)\n";

my $outFile2 = "$outDir/spp.summary";
open OUT, ">$outFile2" or die "Cannot open $outFile2 for writing: $OS_ERROR\n";
for my $sp (sort keys %sppTbl)
{
  my $size = 0;
  for (keys %{$sppTbl{$sp}})
  {
    $size += $sppTbl{$sp}{$_};
  }

  print OUT "$sp\t$size\n";
  for (keys %{$sppTbl{$sp}})
  {
    print OUT "\t$_\t" . $sppTbl{$sp}{$_} . "\n";
  }
}
close OUT;

my $outFile3 = "$outDir/mismatch_spp.summary";
open OUT, ">$outFile3" or die "Cannot open $outFile3 for writing: $OS_ERROR\n";
for my $sp (sort @mismatchedSpp)
{
  my $size = 0;
  for (keys %{$sppTbl{$sp}})
  {
    $size += $sppTbl{$sp}{$_};
  }
  print OUT "$sp\t$size\n";
  for ( sort { $sppTbl{$sp}{$b} <=> $sppTbl{$sp}{$a}} keys %{$sppTbl{$sp}})
  {
    print OUT "\t" . sprintf("%-30s",$_) . "\t" . sprintf("%s",$sppTbl{$sp}{$_}) . "\t" . sprintf("%.1f%%", 100 * $sppTbl{$sp}{$_} / $size ) . "\n";
  }
}
close OUT;

print "\nOutput written to\n$outFile1\n$outFile2\n$outFile3\n\n";

####################################################################
##                               SUBS
####################################################################

# read two column table
sub read2colTbl{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# read three column table
sub read3colTbl{

  my $file = shift;

  my %tx;
  my %er;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t, $e) = split /\s+/,$_;
    $tx{$id} = $t;
    $er{$id} = $e;
  }
  close IN;

  return (\%tx, \%er);
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
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

exit;
