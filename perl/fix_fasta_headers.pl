#! /usr/bin/perl

=head1 NAME

  fix_fasta_headers.pl

=head1 DESCRIPTION

  Change seq header from >seqID';size=\d+; to >seqID' size=\d+
  More precisely, the script changes
  >seqID;str1;str2;...;strN
  to
  >seqID str1 str2 ... strN

=head1 SYNOPSIS

  fix_fasta_headers.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS


=over

=item B<--input-file, -i>
  Input fasta file.

=item B<--output-file, -o>
  Output fasta file.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  fix_fasta_headers.pl -i test.fa -o test_fixed.fa

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
  "input-file|i=s"  => \my $inFile,
  "output-file|o=s" => \my $outFile,
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
  print "ERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outFile)
{
  print "ERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################
my $startRun = time();

open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
$/ = ">";
my $junkFirstOne = <FASTA>;

open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
my $count = 1;

while (<FASTA>)
{
  if ($count % 100 == 0)
  {
    print "\r$count";
  }
  $count++;

  chomp;
  my ($def,@seqlines) = split /\n/, $_;
  my $seq = join '', @seqlines;
  my @f = split /\;/, $def;

  print OUT ">@f\n$seq\n";
}
close OUT;
$/ = "\n";
close FASTA;

print "\rOutput written to $outFile\n";
exit;
