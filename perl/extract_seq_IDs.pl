#! /usr/bin/perl

=head1 NAME

  extract_seq_IDs.pl

=head1 DESCRIPTION

  extract sequence IDs from a fasta file

=head1 SYNOPSIS

  extract_seq_IDs.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS


=over

=item B<--input-file, -i>
  Input fasta file.

=item B<--output-file, -o>
  Output file.

=item B<--quiet>
  Do not print progress messages.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  extract_seq_IDs.pl --quiet -i test.fa -o test.seqIDs

=cut

use strict;
use warnings;
use Pod::Usage;
use List::Util qw( first max min sum );
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-file|i=s"  => \my $inFile,
  "output-file|o=s" => \my $outFile,
  "quiet"           => \my $quiet,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help || !$inFile || !$outFile)
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

my $fileSize = -s $inFile; # file size in bytes

open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
$/ = ">";
my $junkFirstOne = <IN>;
my $count = 0;
while (<IN>)
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

    my $perc = sprintf("%.1f%%", 100 * (tell IN) / $fileSize);
    print "\r$timeStr [$perc]";

  }
  $count++;

  chomp;
  my ($def) = split /\n/, $_;
  my ($id) = split /\s+/, $def;
  print OUT "$id\n";

  if ($id =~ /oxidizing/)
  {
    print "count: $count\tid: $id\tdef: $def\n\n";
  }
}
close OUT;
close IN;

print "\rOutput written to $outFile\n" if !$quiet;

exit 0;
