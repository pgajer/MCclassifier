#!/usr/bin/env perl

=head1 NAME

  multi_select_tx_bigMem_v2.pl

=head1  DESCRIPTION

  A taxon counterpart of multi_select_tx_bigMem_v2.pl putting each lineage file
  in its own directory phGr_dir/phGr.fa

=head1 SYNOPSIS

    multi_select_tx_bigMem_v2.pl -s <file with 2 column table of seqIDs and taxonomies> -i <lineage file>

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Linege/taxon file

=item B<--tx-file, -t>
  A file with 2 or more columns (only the first two columns will be used) of
  sequence IDs and their taxonomy. All sequences of the same taxon will be put in
  the same fasta file in the ouput directory.  file with a list of sequence IDs
  to be selected from the input fasta file

=item B<-h|--help>

    Print help message and exit successfully.

=item B<--quiet>
  Do not print progress messages.

=back

=head1  EXAMPLES

  cd ~/projects/PECAN/data/phylo_groups/v0.3/cx_hb_rdp_FL_5500_phGr_dir

  multi_select_tx_bigMem_v2.pl --debug -t thld_5500_phGr_tbl_v3.txt -i cx_hb_rdp_thld_5500.lineage

=cut

## Copyright (C) 2017 Pawel Gajer pgajer@gmail.com
##
## Permission to use, copy, modify, and distribute this software and its
## documentation with or without modifications and for any purpose and
## without fee is hereby granted, provided that any copyright notices
## appear in all copies and that both those copyright notices and this
## permission notice appear in supporting documentation, and that the
## names of the contributors or copyright holders not be used in
## advertising or publicity pertaining to distribution of the software
## without specific prior permission.
##
## THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
## WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
## CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
## OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
## OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
## OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
## OR PERFORMANCE OF THIS SOFTWARE.
##

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
  "input-file|i=s"  => \my $inFile,
  "tx-file|t=s"     => \my $txFile,
  "dry-run"         => \my $dryRun,
  "debug"           => \my $debug,
  "quiet"           => \my $quiet,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$inFile )
{
  print "\n\nERROR: Missing input file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$txFile )
{
   print "\n\nERROR: Missing selection file\n\n\n";
   pod2usage(verbose => 2,exitstatus => 0);
   exit 1;
}

if ( ! -e $inFile )
{
  warn "\n\n\tERROR: $inFile does not exist";
  print "\n\n";
  exit 1;
}

if ( -l $inFile )
{
  $inFile = readlink($inFile);
}

if ( -l $txFile )
{
   $txFile = readlink($txFile);
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

print "--- Parsing taxon file\n";
my $rspIDsTbl = parse_tx_tbl( $txFile ); # sp => ref to array of seqIDs of the given species
my %spIDsTbl = %{ $rspIDsTbl } ;

my $nSpp = keys %spIDsTbl;

if ( $debug )
{
  print "\nNumber of species in the taxon table: $nSpp\n\n";
}

print "\t--- Parsing lineage file\n";
my %liTbl = read_tbl( $inFile );

print "\r--- Writing species fasta files to outDir\n";
my $spCount = 1;
my $perc;
for my $sp ( keys %spIDsTbl )
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
    $timeStr = "$timeMin:$runTime";
  }
  $perc = sprintf("%.1f%%", 100 * $spCount / $nSpp );
  print "\r" . commify($spCount) . " $timeStr [$perc]";
  $spCount++;

  my $outDir = $sp . "_dir";
  my $outFile = $outDir . "/" . $sp . ".lineage";
  print "\rWriting $sp to $outFile";
  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  for ( @{ $spIDsTbl{$sp} } )
  {
    if ( !exists $liTbl{$_} )
    {
      warn "\n\n\tERROR: $_ does not exist in liTbl";
      print "\n\n";
      exit 1;
    }
    print OUT ">$_\n" . $liTbl{$_} . "\n";
  }
  close OUT;
}

print "\r                                     ";

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


####################################################################
##                               SUBS
####################################################################

##
## parse 3 column tbl
##

## 1642.V1_0	Lactobacillus_iners	0.93
## 0980.V2_1	Lactobacillus_iners	0.97
## 1670.V2_2	Lactobacillus_helveticus_acidophilus	0.98
## 0711.V3_3	Atopobium_vaginae	0.56
## 1149.V1_4	Lactobacillus_iners	0.94
## 1386.V1_5	Prevotella_buccalis	0.85
## 1119.V2_6	BVAB1	0.79
## 1449.V1_7	BVAB1	0.97
## 1600.V1_8	BVAB3	0.93

sub parse_tx_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in parse_tx_tbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %spIDsTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  my $count = 1;
  foreach (<IN>)
  {
    next if /^$/;
    if ( $count % 500 == 0 )
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
	$timeStr = "$timeMin:$runTime";
      }
      print "\r" . commify($count) . " $timeStr";
    }
    $count++;
    chomp;
    my ($id, $sp, $pp) = split /\s+/,$_;
    push @{$spIDsTbl{$sp}}, $id;
  }
  close IN;

  return \%spIDsTbl;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\n"} @{$a};
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my ($file, $skipHeader) = @_;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( $skipHeader )
  {
    my $header = <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# read table with one column
sub readArray
{
  my $file = shift;
  my %rows;

  #local $/ = '\n';
  open (INPUT, "<", $file) or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<INPUT>)
  {
    chomp;
    my ($id) = split /\s+/, $_;
    #print "id: $id\n"; exit 1;
    $rows{$id}=1;
  }
  close INPUT;

  return %rows;
}

# print elements of a hash table
sub printTbl
{
  my $rTbl = shift;
  map {print "|$_|\t|" . $rTbl->{$_} . "|\n"} keys %$rTbl;
}

# print n keys of a hash table
sub printKeys
{
  my ($rTbl, $n) = @_;
  my $i = 0;
  map {$i++; print "$_\n" if $i < $n} keys %$rTbl;
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

exit 0;
