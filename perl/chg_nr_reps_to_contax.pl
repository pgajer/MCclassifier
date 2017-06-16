#!/usr/bin/env perl

=head1 NAME

  chg_nr_reps_to_contax.pl

=head1 DESCRIPTION

  Change representative sequences of 100% identity clusters to contax seq's.

=head1 SYNOPSIS

  chg_nr_reps_to_contax.pl -i <cltr2 file> -o <seqIds file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

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

  chg_nr_reps_to_contax.pl -i Banfield_medoids_FL_QCed_tr_v2.clstr2 -o Banfield_medoids_FL_QCed_tr_v2.seqIDs

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd qw(abs_path);
use File::Temp qw/ tempfile /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-file|i=s"  => \my $inFile,
  "output-file|o=s" => \my $outFile,
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

my %tbl = parse_clstr2_file( $inFile );
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
my $nChanges = 0;
for my $refId ( keys %tbl )
{
    my @ids = split ",", $tbl{$refId};
    if ( $refId !~ /^Con/ && @ids > 1 )
    {
        my @cIds = grep{ /^Con/ } @ids;
        if ( @cIds > 0 )
        {
            $refId = shift @cIds;
            $nChanges++;
        }
    }
    print OUT "$refId\n";
}
close OUT;

print "\nNumber of changes made: $nChanges\n\n";

####################################################################
##                               SUBS
####################################################################

# parse a clstr2 file
sub parse_clstr2_file
{
  my $inFile = shift;

  my %tbl;
  open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
  foreach my $rec (<IN>)
  {
    chomp $rec;
    my @ids = split ",", $rec;
    my $refId = shift @ids;
    $tbl{$refId} = \@ids;
  }
  close IN;

  return %tbl;
}

exit 0;
