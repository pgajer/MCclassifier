#!/usr/bin/env perl

=head1 NAME

  multi_select_seqs_bigMem.pl

=head1  DESCRIPTION

  A version of multi_select_seqs.pl that reads the whole fasta file into memory
  and then prints out one species at a time to its fasta file.

=head1 SYNOPSIS

    multi_select_seqs_bigMem.pl -s <file with 2 column table of seqIDs and taxonomies> -i <fastaFile> -o <outDir>

=head1 OPTIONS

=over

=item B<--fasta-file, -i>
  Fasta file

=item B<--tx-file, -t>
  A file with 2 or more columns (only the first two columns will be used) of
  sequence IDs and their taxonomy. All sequences of the same taxon will be put in
  the same fasta file in the ouput directory.  file with a list of sequence IDs
  to be selected from the input fasta file

=item B<--output-dir, -o>
  Output diretory holding all taxon's fasta files.

=item B<-h|--help>

    Print help message and exit successfully.

=item B<--quiet>
  Do not print progress messages.

=back

=head1  EXAMPLES

  multi_select_seqs_bigMem.pl --debug -t /local/scratch/MM/MM_MC_order7_results.txt -i /local/scratch/MM/MM_all_updated_seqID_no_ctrls.fa -o mm_spp_dir

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
  "output-dir|o=s"  => \my $outDir,
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
elsif ( !$outDir )
{
   print "\n\nERROR: Missing output directory\n\n\n";
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

if ( ! -e $outDir )
{
  my $cmd = "mkdir -p $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?";# if !$dryRun;
}

print "--- Parsing taxon file\n";
my $rspIDsTbl = parse_tx_tbl( $txFile ); # sp => ref to array of seqIDs of the given species
my %spIDsTbl = %{ $rspIDsTbl } ;

my $nSpp = keys %spIDsTbl;

if ( $debug )
{
  print "\nNumber of species in the taxon table: $nSpp\n\n";
}

print "\t--- Parsing big fat fasta file\n";
my %faTbl; # seqID => its sequence
my $count = 0;
my $selCount = 0;
open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
$/ = ">";
my $junkFirstOne = <FASTA>;
while (<FASTA>)
{
  chomp;
  if ( $count % 500 == 0 )
  {
    print "\r$count" if $debug;
  }
  $count++;

  my ($def,@seqlines) = split /\n/, $_;
  my $seq = join '', @seqlines;
  my ($id) = split /\s+/, $def;
  $faTbl{$id} = $seq;
}
close FASTA;


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

  my $outFile = $outDir . "/" . $sp . ".fa";
  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  for ( @{ $spIDsTbl{$sp} } )
  {
    print OUT ">$_\n" . $faTbl{$_} . "\n";
  }
  close OUT;
}

if (!$quiet)
{
  print "\r\n\n\tFasta files of selected sequences written to $outDir\n\n";
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
      print "\r" . commify($spCount) . " $timeStr";
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
