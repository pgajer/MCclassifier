#!/usr/bin/env perl

=head1 NAME

  trim_align.pl

=head1 DESCRIPTION

  trim Multiple Sequence Alignment from 'start' to 'end' position in the alignment

=head1 SYNOPSIS

  trim_align.pl -i <input file> -o <output file> -s <start_end> [Options]

=head1 OPTIONS


=over

=item B<--sel-seqs-range, -s>
  1-based range of selected sequences separated by '_', as in 399_405

=item B<--input-file, -i>
  Input Multiple Sequence Alignment fasta file.

=item B<--output-file, -o>
  Output Multiple Sequence Alignment fasta file.

=item B<--summary-file, -j>
  mothur summary file generated by summary.seqs()

  $ head Lactobacillaceae_V3V4_good95_wStrepOG_SATe.summary

  seqname	      start	end	nbases	ambigs	polymer	numSeqs
  JQ517277.1.1474	2	640	429	0	4	1
  KC561106.1.1500	3	640	427	0	5	1
  EF120375.1.1550	2	640	428	0	5	1
  HM217944.1.1489	2	640	428	0	6	1


=item B<--criteria, -c>
  An integer between 0 and 100. The alignment will be trimmed at the positions
  where percentage of sequences that start and end is at least the criteria value.

=item B<--quiet>
  Do not print progress messages.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/devel/MCclassifier/data/SILVA_123/Lactobacillaceae_ss100_dir/Lactobacillaceae_V3V4_good95_wStrepOG_dir
  trim_align.pl --verbose -c 95 -j Lactobacillaceae_V3V4_good95_wStrepOG_SATe.summary -i Lactobacillaceae_V3V4_good95_wStrepOG_SATe.good.aln -o Lactobacillaceae_V3V4_good95_wStrepOG_SATe_good_aln_trimmed.fa

  cd ~/devel/MCclassifier/data/SILVA_123/Lactobacillaceae_ss100_dir
  trim_align.pl -c 95 -j satejob.marker001.Lactobacillaceae_good_full.fa.good95.summary -i satejob.marker001.Lactobacillaceae_good_full.fa.good95.aln -o satejob.marker001.Lactobacillaceae_full_good95_aln_trimmed.fa

  Alignment length: 2330
  Alignment trimmed at the positions: 115, 2173
  Output written to satejob.marker001.Lactobacillaceae_full_good95_aln_trimmed.fa

  cd ~/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013
  trim_align.pl -i /Users/pgajer/projects/16S_rRNA_pipeline/vaginal-0.2.1.refpkg/vaginal_aln.fasta -o vaginal_V1V3.align.fa -s 245_763

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::Util qw( sum );

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-file|i=s"     => \my $inFile,
  "output-file|o=s"    => \my $outFile,
  "sel-seqs-range|s=s" => \my $selRange,
  "summary-file|j=s"   => \my $summaryFile,
  "criteria|c=i"       => \my $criteria,
  "quiet"              => \my $quiet,
  "verbose|v"          => \my $verbose,
  "help|h!"            => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$inFile)
{
  print "ERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$outFile)
{
  print "ERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
# elsif (!$selRange)
# {
#   print "ERROR: Missing selection range\n\n";
#   pod2usage(verbose => 2,exitstatus => 0);
#   exit 1;
# }

if ( ! -f $inFile )
{
  print "\n\nERROR: $inFile does not exist\n\n\n";
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

my $startPos;
my $endPos;

if ( $summaryFile )
{
  if ( ! $criteria )
  {
    print "\n\nERROR: $summaryFile present but criteria not defined.\n\n\n";
    exit 1;
  }

  if ( $criteria < 0 || $criteria > 100 )
  {
    print "\n\nERROR: $criteria has to be an integer between 0 and 100.\n\n\n";
    exit 1;
  }

  ($startPos, $endPos) = readSummaryFile($summaryFile);
}
else
{
  ($startPos, $endPos) = split '_', $selRange;
}

print "\nstartPos: $startPos\nendPos: $endPos\n" if $verbose;

## print "--- Generating fasta files of selected clusters ... ";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
open (IN, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
$/ = ">";
my $junkFirstOne = <IN>;
my $count = 1;
my $step = 100;
my $seq;
while (<IN>)
{
  if ( !$quiet && ($count % $step == 0) )
  {
    print "\r$count";
  }
  $count++;

  chomp;

  my ($def,@seqlines) = split /\n/, $_;
  $seq = join '', @seqlines;
  my ($id) = split /\s+/, $def;
  my $s = substr( $seq, $startPos-1, $endPos - $startPos + 1);
  print OUT ">$id\n$s\n";
}
$/ = "\n";
close OUT;
close IN;

if ( !$quiet )
{
  print "\r                                          \n";
  print "\n\tAlignment length: " . length($seq) . "\n";
  print "\tAlignment trimmed at the positions: $startPos, $endPos\n";
  print "\tOutput written to $outFile\n\n";
}


####################################################################
##                               SUBS
####################################################################

# parse summary file, extract only start/end position data
sub readSummaryFile
{
  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR: $file does not exist\n\n\n";
    exit 1;
  }

  my %startTbl;
  my %endTbl;
  open IN, "<$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  my $headerStr = <IN>;
  my $nSeqs = 0;
  foreach  (<IN>)
  {
    chomp;
    $nSeqs++;
    my ($id, $start, $end) = split /\s+/, $_;
    ##print "line: $_\nid: $id  start: $start  end: $end\n"; exit 1;
    $startTbl{$start}++;
    $endTbl{$end}++;
  }
  close IN;

  # printTbl(\%startTbl, "startTbl");
  # printTbl(\%endTbl, "endTbl");
  # exit 1;

  #### ---------- start position -------------
  my @startPerc; # array of cummulative start position percentages, it ends when
                 # the cummulative percentage reaches 100%

  my %cummStartTbl;
  my $cummStart = 0;
  my @startPositions = sort { $a <=> $b } keys %startTbl;
  my @startFreqs = @startTbl{@startPositions};

  if ($verbose)
  {
    printArray(\@startPositions, "startPositions");
    printArray(\@startFreqs, "startFreqs");
  }

  my @cummStartFreqs;
  for my $i ( 0..$#startPositions )
  {
    $cummStartFreqs[$i] = sum( @startTbl{@startPositions[0..$i]} );
  }

  if ($verbose)
  {
    printArray(\@cummStartFreqs, "cummStartFreqs");
  }

  my $criteriaNum = (100 - $criteria) * $nSeqs/100.0;
  # print "criteriaNum: $criteriaNum\n";

  my $startPos;
  for my $i ( 0..$#startPositions )
  {
    if ( $cummStartFreqs[$i] >= $criteriaNum )
    {
      $startPos = $startPositions[$i];
      last;
    }
  }
  # print "startPos: $startPos\n";

  #### ---------- end position -------------
  my @endPerc; # array of cummulative end position percentages, it ends when
               # the cummulative percentage reaches 100%

  my %cummEndTbl;
  my $cummEnd = 0;
  my @endPositions = sort { $b <=> $a } keys %endTbl; ## note sorting in the decreasing order
  my @endFreqs = @endTbl{@endPositions};

  if ($verbose)
  {
     printArray(\@endPositions, "endPositions");
     printArray(\@endFreqs, "endFreqs");
  }

  my @cummEndFreqs;
  for my $i ( 0..$#endPositions )
  {
    $cummEndFreqs[$i] = sum( @endTbl{@endPositions[0..$i]} );
  }

  if ($verbose)
  {
    printArray(\@cummEndFreqs, "cummEndFreqs");
  }

  my $endPos;
  for my $i ( 0..$#endPositions )
  {
    if ( $cummEndFreqs[$i] >= $criteriaNum )
    {
      $endPos = $endPositions[$i];
      last;
    }
  }

  return ($startPos, $endPos);
}

# print elements of a hash table
sub printTbl{

  my ($rTbl, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# print array to stdout
sub printArray{

  my ($a, $header) = @_;
  print "$header\n" if $header;
  map {print "$_,"} @{$a};
  print "\n";
}

exit 0;
