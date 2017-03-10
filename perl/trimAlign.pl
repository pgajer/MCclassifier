#! /usr/bin/perl

=head1 NAME

  trimAlign.pl

=head1 DESCRIPTION

  Trim multiple sequence alignment reads to the clear range
  and remove columns with more that <nGaps> gaps in them,
  where <nGaps> = int( nSeqs * <percGaps> / 100 ),
  if <percGaps> is defined.

=head1 SYNOPSIS

  trimAlign.pl -i <multiple sequence alignment file> -o <output file> [-g <percGaps>] [-h]

=head1 OPTIONS


=over

=item B<--clr-rage-type, -c>
  type of clear range. Supported types: median, min. Default: min.

=item B<--msa-file, -i>
  multiple sequence alignment file

=item B<--max-perc-gaps, -g>
  maximum percentage of gaps allowed in the multiple sequence alignment (MSA)

=item B<--rm-gap-cols, -G>
  remove columns consisting of only gaps

=item B<--output-file, -o>
  Output file.

=item B<--print-removed-pos, -p>
  print 1-based positions in the alignment that were removed

=item B<--seq-ranges-file, -r>
  print 1-based start and end positions of each sequence in the alignment to a file

=item B<--sel-seqs-range, -s>
  1-based range of selected sequences separated by '_', as in 399_405

=item B<--trim-range, -t>
  0-based trim range of selected sequences separated by '_', as in 398_404

=item B<--outlier-end, -e>
  0-based position 'pos' such that if a sequence's end position is less than
  'pos', than the sequence is removed

=item B<--quiet>
  Do not print progress messages.

=item B<-h|--help>
    Print help message and exit successfully.

=back


=head1 EXAMPLE

  pushd /Users/pgajer/projects/douchingStudy/speciationData_july10_trim220/Lactobacillus/Lactobacillus.ribotypeTbls500
  trimAlign.pl -i L.iners_top500.filter.fasta -o L.iners_top500.filter.fasta.trimmed --seq-ranges-file L.iners_top500.filter.fasta.seqRange

  pushd /Users/pgajer/projects/douchingStudy/speciationData_april23/Lactobacillus
  trimAlign.pl -i Lactobacillus.cloneTbls/L.gasseri_top20.align.fa -o Lactobacillus.cloneTbls/L.gasseri_top20_trimmed.fa --max-perc-gaps 60

  pushd /Users/pgajer/devel/speciesClassifier/data/U01_400s/speciationData_394s_Oct7_2009/Ignavigranum
  trimAlign.pl -i Ignavigranum_nr.align.fa --seq-ranges-file Ignavigranum_nr.seqRange

  (see seqRange.R)

  Start       End
  Min.   :0   Min.   : 813
  1st Qu.:0   1st Qu.:1101
  Median :0   Median :1121
  Mean   :0   Mean   :1108
  3rd Qu.:0   3rd Qu.:1135
  Max.   :0   Max.   :1150


  trimAlign.pl -i Ignavigranum_nr.align.fa -t 0_1100 -e 1000 -o Ignavigranum_nr.align.trimmed.fa



  pushd /home/pgajer/speciesClassifier/data/U01_400s/394s_comb (on IGS network)
  trimAlign.pl -i comb_nr.align.filtered --seq-ranges-file comb_nr.seqRange


  pushd ~/devel/speciesClassifier/data/RAKAI/speciationData_Dec16.2009/Lactobacillus
  trimAlign.pl -i Lactobacillus_plusTrgSeqs_nr.algn.fa --seq-ranges-file Lactobacillus_plusTrgSeqs_nr.algn.range -s 3989_4315

  trimAlign.pl -i Lactobacillus_plusTrgSeqs_nr.algn.fa -o Lactobacillus_plusTrgSeqs_nr.algn.trimmed.fa -t 2_598


  pushd /Users/pgajer/devel/speciesClassifier/data/douchingStudy/speciationData_Sep1.2009/Lactobacillus/Lactobacillus.richness
  trimAlign.pl -i L_crispatus_top500.fa -o L_crispatus_top500_trimmed2.fa --print-removed-pos

=cut

use strict;
use warnings;
use Pod::Usage;
use List::Util qw( max min );
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Perl;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::LocatableSeq;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

# my $percGaps = 20;
my $selStart; # start index of $selSeqs
my $selEnd;   # end index of $selSeqs
my $trimStart;
my $trimEnd;
my $clrRangeType = "min";

GetOptions(
  "clr-rage-type|c=s"    => \$clrRangeType,
  "msa-file|i=s"         => \my $msaFile,
  "max-perc-gaps|g=s"    => \my $percGaps,
  "rm-gap-cols|G=s"      => \my $rmGapCols,
  "print-removed-pos"    => \my $printRemovedPos,
  "output-file|o=s"      => \my $outFile,
  "seq-ranges-file|r=s"  => \my $seqRangesFile,
  "sel-seqs-range|s=s"   => \my $selRange,
  "trim-range|t=s"       => \my $trimRange,
  "outlier-end|e=s"      => \my $outlierEnd,
  "min-seq-len|l=i"      => \my $minLen,
  "quiet"                => \my $quiet,
  "help|h!"              => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help || !$msaFile)
{
    pod2usage(verbose => 2,exitstatus => 0);
    exit;
}
elsif (!$minLen)
{
  print "ERROR: Missing sequence min length threshold\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################


# get alignment format
my $format = guess_format($msaFile);

# read MSA
print "Parsing alignment ... " if !$quiet;
my $seqIO = Bio::AlignIO->new(-file => $msaFile,
			      -format=>$format);

my $aln = $seqIO->next_aln; # a Bio::Align::AlignI object
                            # ther is only one alignment

my $lastEl = $aln->length - 1;
my $nSeqs  = $aln->num_sequences();
print "done\n" if !$quiet;

# if ($rmGapCols)
# {
#   exit;
# }

if ($selRange)
{
  ($selStart, $selEnd) = split '_', $selRange;
  # print "selStart=$selStart\tselEnd=$selEnd\n";  exit;
}

if ($trimRange)
{
  ($trimStart, $trimEnd) = split '_', $trimRange;
}


my @alnTbl;
my @seqIds;
my @starts;
my @ends;
my %selSeqs;
my $nSel = 0;

print "Determining alignemnt clear range ... " if !$quiet;
foreach my $i (1..$nSeqs)
{
  if ($i % 10 == 0)
  {
    print "\rDetermining clear alignemnt range ... $i/$nSeqs" if !$quiet;
  }

  my $locseq = $aln->get_seq_by_pos($i); # a Bio::LocatableSeq object
  # Gets a sequence based on its position in the alignment.
  # Numbering starts from 1.  Sequence positions larger than
  # no_sequences() will throw an error.

  push @seqIds, $locseq->id();

  my @seq = split "", $locseq->seq();

  # find first non-gapped position
  my $startPos = 0;
  while ( $seq[$startPos] !~ /\w/ )
  {
    $startPos++;
  }

  # find last non-gapped position
  my $endPos = $lastEl;
  while ( $seq[$endPos] !~ /\w/ )
  {
    $endPos--;
  }

  if ( ($selRange && ($i < $selStart || $i > $selEnd)) ||
       ($outlierEnd && $endPos < $outlierEnd) )
  {
    $selSeqs{$locseq->id()}=0;
  }
  else
  {
    $selSeqs{$locseq->id()}=1;
  }

  if ( $selSeqs{$locseq->id()} )
  {
    $alnTbl[$nSel] = \@seq;
    $nSel++;

    push @starts, $startPos;
    push @ends, $endPos;
  }
}
print "\rDetermining clear alignemnt range ... done\n" if !$quiet;

if ($selRange && $nSel != ($selEnd - $selStart+1))
{
  warn "Warning!: nSel!=(selEnd-selStart+1); nSel=$nSel\tselEnd=$selEnd\tselStart=$selStart";
}

if ($seqRangesFile)
{
  open OUT, ">$seqRangesFile" or die "Cannot open $seqRangesFile for reading: $OS_ERROR\n";
  foreach my $i (0..$#starts)
  {
    print OUT $starts[$i] . "\t" . $ends[$i] . "\n";
  }
  close OUT;
  print "Sequence ranges written to $seqRangesFile\n" if !$quiet;
}

# remove columns of the MSA that contain more than or equal to
# int( nSeqs * <percGaps> / 100 ) gaps
# where nSeqs is the number of sequences in the MSA

my @goodPos;

if ($percGaps)
{
  my $maxGaps = int( $nSeqs * $percGaps / 100 );
  my @gapPos;    # 0-based positions with more than $maxGaps
  my @selGapPos; # 0-based positions with more than $maxGaps and all gaps in the selected sequences

  print "Finding gapped columns ... " if !$quiet;
#  foreach my $pos ($maxStart..$minEnd)
  foreach my $pos (1..$nSeqs)
  {
    if ($pos % 10 == 0)
    {
      print "\rFinding gapped columns ... $pos" if !$quiet;
    }

    my %count;
    my %selCount;

    foreach my $seq ($aln->each_seq)
    {
      my $res = $seq->subseq(($pos+1), ($pos+1));
      $count{$res}++;

      if ($selSeqs{$seq->id()})
      {
	$selCount{$res}++;
      }
    }

    my @gapChars = grep {!/\w/} keys %count;

    if (@gapChars > 1)
    {
      print "ERROR in at line " . __LINE__ . ": gapChars cannot have more than one element!\n";
      print "gapChars: ";
      foreach (@gapChars) {

	print "$_ ";
      }
      print "\n";
      exit;
    }

    if ( @gapChars < 1 || $count{shift @gapChars} < $maxGaps )
    {
      push @goodPos, $pos;
    }
    else
    {
      push @gapPos, $pos;
    }

    if ($selRange)
    {
      my @selGapChars = grep {!/\w/} keys %selCount;
      # my $n = keys %selCount; my $s = @selGapChars; my $c = shift @selGapChars; print "n=$n\ts=$s\n"; print "c=$c\n"; exit;
      my $c;

      if ( defined ($c = shift @selGapChars) && exists $count{$c} && $count{$c} > $maxGaps && $nSel == @selGapChars)
      {
	push @selGapPos, $pos;
      }
    }
  }
  print "\rFinding gapped columns ... done\n" if !$quiet;

  if ($selRange)
  {
    my $selGapFile = "selGaps.txt";
    open OUT, ">$selGapFile" or die "Cannot open $selGapFile for reading: $OS_ERROR\n";
    foreach (@selGapPos)
    {
      print OUT "$_\n";
    }
    close OUT;
    print "A list of gapped positions (with gaps in all selected sequences) was written to $selGapFile\n" if !$quiet;
  }
}

if ($trimRange)
{
  @goodPos = $trimStart..$trimEnd;
  print "clear range: trimStart=$trimStart\ttrimEnd=$trimEnd\n" if !$quiet;
}
elsif ( $clrRangeType eq "min") # clear range selection
{
  my $maxStart = max( @starts );
  my $minEnd   = min( @ends );
  @goodPos = $maxStart..$minEnd;
  print "clear range: maxStart=$maxStart\tminEnd=$minEnd\n" if !$quiet;
}
else
{
  my $medStart = median(\@starts);
  my $medEnd   = median(\@ends);
  @goodPos = $medStart..$medEnd;  # 0-based positions with no more than $maxGaps
  print "clear range: medStart=$medStart\tmedEnd=$medEnd\n" if !$quiet;
}

if ($outFile)
{
  my $base = basename($outFile,(".fa",".fasta"));
  my $dir  = dirname($outFile);
  my $noGapFile = $dir . "/" . $base . "_noGap.fa";

  print "Writting trimmed alignment to $outFile ... " if !$quiet;
  # create an output Bio::SeqIO object
  my $out = Bio::SeqIO->new(-file => ">$outFile" ,
			    -format => 'Fasta');

  my $outNoGap = Bio::SeqIO->new(-file => ">$noGapFile" ,
				 -format => 'Fasta');

  # generate trimmed aligment with no gapped columns
  foreach my $i (1..$nSel)
  {
    my $seqStr = join "", @{$alnTbl[$i-1]}[@goodPos];

    #print "\nbefore tr seqStr: $seqStr\n";
    my $seqStrNoGaps = ($seqStr =~ s/-//g);
    $seqStrNoGaps =~ s/\.//g;

    print "seqStr: $seqStr\nseqStrNoGaps: $seqStrNoGaps";
    exit if $i > 100;

    if (length($seqStrNoGaps) > $minLen)
    {
      my $seq = Bio::Seq->new( -display_id => $seqIds[$i-1],
			       -seq        => $seqStr);
      $out->write_seq($seq);

      $seq->seq($seqStrNoGaps);
      $outNoGap->write_seq($seq);
    }
  }
  print "done\n" if !$quiet;

  if ( defined $printRemovedPos ) {

    my %goodPoss = map{$_ => 1} @goodPos;
    my @badPos = grep !exists $goodPoss{$_}, 0..$lastEl;

    print "removed positions:"Y
    foreach (@badPos) {

      print " $_" if !$quiet;
    }
    print "\n" if !$quiet;
  }

  print "Output written to:\n\t$outFile\nand\n\t$noGapFile\n" if !$quiet;
}



####################################################################
##                               SUBS
####################################################################

# median of an array
sub median
{
    my ($aRef) = @_;
    my $aSize  = @{ $aRef };

    if ( $aSize == 1 )
    {
        return $aRef->[0];
    }
    elsif ( $aSize == 2 )
    {
        return ( $aRef->[0] + $aRef->[1] ) / 2;
    }

    my @b = sort { $a <=> $b } @{ $aRef };

    if( $aSize % 2 == 0 )
    {
	return ( $b[$aSize/2-1] + $b[$aSize/2] ) / 2;
    }
    else
    {
	return $b[$aSize/2];
    }
}

# minor modification of _guess_format() of Bio::AlignIO
sub guess_format {

    return unless $_ = shift;

    return 'fasta'   if /\.(fasta|fast|seq|fa|fsa|nt|aa)$/i;
    return 'msf'     if /\.(msf|pileup)$/i;
    return 'pfam'    if /\.(pfam|pfm)$/i;
    return 'selex'   if /\.(selex|slx|selx|slex|sx)$/i;
    return 'phylip'  if /\.(phylip|phlp|phyl|phy|phy|ph)$/i;
    return 'nexus'   if /\.(nexus|nex)$/i;
    return 'mega'    if /\.(meg|mega)$/i;
}

exit;
