#! /usr/bin/perl

=head1 NAME

  select_seqs.pl

=head1  DESCRIPTION

  It is a faster version of select_seqs.pl

  select sequences from a single fasta file given a list of sequence IDs
  selected sequences are place in the new fasta file in the order they are listed

=head1  INPUT

    file with readIds of sequences to be selected
    fasta file to select from

=head1  OUTPUT

    fasta file of selected sequences

=head1 OPTIONS

=over

=item B<-h|--help>

    Print help message and exit successfully.

=item B<--fasta-file, -i>

    fasta file

=item B<--sels-file, -s>

  a file with a list of sequence IDs to be selected from the input fasta file

=item B<--excl-file, -e>

  a file with a list of sequence IDs to be excluded from the input fasta file

=item B<--output-file, -o>

    Output file.

=item B<--quiet>
  Do not print progress messages.

=back

=head1 SYNOPSIS

    select_seqs.pl -s <list of readIDs file> -i <fastaFile> -o <outFile>


=head1  EXAMPLES

    select_seqs.pl -s Unclassified_nr.1k.sample.align.cmd.3d.cl2 -i Unclassified_nr.1k.sample.align.fa  -o Unclassified_nr.1k.sample.align.cmd.3d.cl2.fa

=cut

## Copyright (C) 2015 Pawel Gajer pgajer@gmail.com
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
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

$OUTPUT_AUTOFLUSH = 1;

GetOptions(
  "h|help"          => \my $help,
  "fasta-file|i=s"  => \my $inFile,
  "sels-file|s=s"   => \my $selsFile,
  "excl-file|e=s"   => \my $exclFile,
  "output-file|o=s" => \my $outFile,
  "quiet"           => \my $quiet,
  )
  || pod2usage();


if ( $help || !$inFile || !$outFile )
{
  pod2usage( { -exitval => 0, -verbose => 2 } );
  exit;
}

if ( !$selsFile && !$exclFile )
{
  print "Either -s or -e option has to be used.\n";
  pod2usage( { -exitval => 0, -verbose => 2 } );
  exit;
}

if ( -l $inFile )
{
  $inFile = readlink($inFile);
}

if ( $selsFile && -l $selsFile )
{
  $selsFile = readlink($selsFile);
}

if ( $exclFile && -l $exclFile )
{
  $exclFile = readlink($exclFile);
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

my $fileSize = -s $inFile; # file size in bytes

open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
my $count = 0;
my $selCount = 0;

if ($selsFile)
{
  print "\r--- Reading seq IDs from $selsFile" if !$quiet;

  my %selRead = readArray($selsFile);

  my $nSelSeqs = keys %selRead;

  # printKeys(\%selRead,100);
  # exit;
  # printTbl(\%selRead); exit;
  # print "size(selRead): " . keys(%selRead) . "\n"; exit;

  if (!$quiet)
  {
    print "\r                                  ";
  }

  $/ = ">";
  my $junkFirstOne = <FASTA>;
  while (<FASTA>)
  {
    if ( !$quiet && ($count % 500 == 0) )
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

      my $perc = sprintf("%.1f%%", 100 * (tell FASTA) / $fileSize);
      print "\r$timeStr [$perc]";

    }
    $count++;

    chomp;
    my ($def,@seqlines) = split /\n/, $_;
    my $seq = join '', @seqlines;
    my ($id, $size) = split /;|\s+/, $def;
    ##print "id=$id\tdef=$def\n";exit;
    ##$size =~ s/size=//;
    ##print "size=$size\n"; exit;

    if ( exists $selRead{$id} && $selRead{$id}==1 )
    {
      if ($size)
      {
	print OUT ">$def\n$seq\n";
      }
      else
      {
	print OUT ">$id\n$seq\n";
      }
      $selCount++;
      $selRead{$id}++;

      last if $selCount==$nSelSeqs;
    }
  } ## end of while()

  if ($selCount!=$nSelSeqs)
  {
    print "\r                                                           ";
    warn "\n\tWARNING: Number of sequences selected:       $selCount";
    warn "  \tWARNING: Number of sequences to be selected: $nSelSeqs";
  }
}
else
{
  my %exclRead = readArray($exclFile);

  $/ = ">";
  my $junkFirstOne = <FASTA>;

  while (<FASTA>)
  {
    if (!$quiet && ($count % 500 == 0))
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

      my $perc = sprintf("%.1f%%", 100 * (tell FASTA) / $fileSize);
      print "\r$timeStr [$perc]";

    }
    $count++;

    chomp;
    my ($def,@seqlines) = split /\n/, $_;
    my $seq = join '', @seqlines;
    my ($id) = split /\s+/, $def;
    if ( !exists $exclRead{$id} )
    {
      print OUT ">$id\n$seq\n";
      $selCount++;
    }
  }
}
close OUT;
$/ = "\n";

if (!$quiet)
{
  print "\r                                                           ";
  print "\n\tNumer of sequences in the input fasta file:  $count\n";
  print "  \tNumer of sequences in the output fasta file: $selCount\n\n";

  print "\rSelected sequences written to $outFile\n";

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
}


####################################################################
##                               SUBS
####################################################################

# read table with one column
sub readArray{

  my $file = shift;
  my %rows;

  $/ = "\n";
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id) = split /\s+/, $_;
    #print "id: $id\n"; exit;
    $rows{$id}=1;
  }
  close IN;
#  $/ = ">";

  return %rows;
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "|$_|\t|" . $rTbl->{$_} . "|\n"} keys %$rTbl;
}

# print n keys of a hash table
sub printKeys{

  my ($rTbl, $n) = @_;
  my $i = 0;
  map {$i++; print "$_\n" if $i < $n} keys %$rTbl;
}

exit;
