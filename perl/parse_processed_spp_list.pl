#!/usr/bin/env perl

=head1 NAME

  parse_processed_spp_list.pl

=head1 DESCRIPTION

  Parse

  /Users/pgajer/projects/M_and_M/new_16S_classification_data/mm_validate_reports_dir_2017-04-29_01_02_13/processed_spp_list.txt

  File format

  --
  Enterococcus_termitis

  n:     1
  n(nr): 1
  --
  Enterococcus_asini

  n:     1
  n(nr): 1


  Skipping blank lines, "--", n: and n(nr): lines. Recording only species.

=head1 SYNOPSIS

  parse_processed_spp_list.pl

=head1 OPTIONS

=over

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  parse_processed_spp_list.pl

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
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}


####################################################################
##                               MAIN
####################################################################

my $inFile = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mm_validate_reports_dir_2017-04-29_01_02_13/processed_spp_list.txt";

open (IN, "<", "$inFile") or die "Cannot open $inFile for reading: $OS_ERROR";
while (<IN>)
{
  next if /^$/;
  next if /^n/;
  next if /^--/;
  ##chomp;
  print "$_";
}
close IN;

####################################################################
##                               SUBS
####################################################################

exit 0;
