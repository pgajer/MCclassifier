#!/usr/bin/env perl

=head1 NAME

  tx_freq.pl

=head1 DESCRIPTION

  This script generates a table of taxon frequencies. Input file is a two column
  table with seqIDs in the first column and taxons/or whatever in the second
  column. Frequency table of the labels present in the second column is generated.

=head1 SYNOPSIS

  tx_freq.pl -i <taxon file> -o <freq file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

=item B<--print>
  Print the frequency table to stdout.

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

  cd ~/projects/PECAN/data/Banfield_contax/FL

  tx_freq.pl --print -i cxhb_FL_QCed_trAlgn_phylum.tx -o cxhb_FL_QCed_trAlgn_phylum.freq

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
    "input-file|i=s"  => \my $inFile,
    "output-file|o=s" => \my $outFile,
    "print"           => \my $print,
    "verbose|v"       => \my $verbose,
    "debug"           => \my $debug,
    "dry-run"         => \my $dryRun,
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
  print "\n\n\tERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outFile )
{
  print "\n\n\tERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $inFile )
{
  warn "\n\n\tERROR: $inFile does not exist";
  print "\n\n";
  exit;
}

####################################################################
##                               MAIN
####################################################################

my %tx = read_tbl( $inFile );

my %freq; ## table of number of sequences per species
map { $freq{$_}++ } values %tx;

if ( $print )
{
    print_formated_freq_tbl( \%freq );
}

write_formated_freq_tbl( \%freq, $outFile );

####################################################################
##                               SUBS
####################################################################

# this is a version of print_formated_tbl() where sorting w/r values are
# performed within this routine
sub write_formated_freq_tbl
{
  my ($rTbl, $outFile ) = @_;

  my @args = sort { $rTbl->{$b} <=> $rTbl->{$a} } keys %{$rTbl};
  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  for (@args)
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "WARNING: tbl value not defined for $_\n" if !exists $rTbl->{$_};
    print OUT "$_$pad" . $rTbl->{$_} . "\n";
  }
  close OUT;

  return @args;
}

# this is a version of print_formated_tbl() where sorting w/r values are
# performed within this routine
sub print_formated_freq_tbl
{
  my $rTbl = shift;

  my @args = sort { $rTbl->{$b} <=> $rTbl->{$a} } keys %{$rTbl};
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
    print "WARNING: tbl value not defined for $_\n" if !exists $rTbl->{$_};
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  print "\n";

  return @args;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my $file = shift;

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

exit 0;
