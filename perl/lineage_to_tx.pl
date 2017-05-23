#!/usr/bin/env perl

=head1 NAME

  lineage_to_tx.pl

=head1 DESCRIPTION

  Generates seqID => species taxonomy table from a lineage table

=head1 SYNOPSIS

  lineage_to_tx.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS

=over

=item B<--li-file, -i>
  Lineage file.

=item B<--tx-file, -o>
  Taxon file file.

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

  lineage_to_tx.pl -i master_V3V4_no_outliers.lineage -o master_V3V4_no_outliers.tx

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
  "li-file|i=s"     => \my $liFile,
  "tx-file|o=s"     => \my $txFile,
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

if ( !$liFile )
{
  print "\n\n\tERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$txFile )
{
  print "\n\n\tERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $liFile )
{
  warn "\n\n\tERROR: $liFile does not exist";
  print "\n\n";
  exit;
}

####################################################################
##                               MAIN
####################################################################

open IN, "$liFile" or die "Cannot open $liFile for reading: $OS_ERROR";
open OUT, ">$txFile" or die "Cannot open $txFile for writing: $OS_ERROR\n";
for my $lineage (<IN>)
{
  chomp $lineage;
  my ($id, $li) = split /\s+/, $lineage;
  my @f = split ";", $li;
  my $sp = pop @f;
  print OUT "$id\t$sp\n";
}
close IN;
close OUT;

exit 0;
