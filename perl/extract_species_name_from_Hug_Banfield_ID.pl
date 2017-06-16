#!/usr/bin/env perl

=head1 NAME

  extract_species_name_from_Hug_Banfield_ID.pl

=head1 DESCRIPTION

  Extract species name from Hug Banfield seqID

=head1 SYNOPSIS

  extract_species_name_from_Hug_Banfield_ID.pl [Options]

=head1 OPTIONS

=over

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

  cd /Users/pgajer/projects/PECAN/data/Banfield_contax

  extract_species_name_from_Hug_Banfield_ID.pl

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


####################################################################
##                               MAIN
####################################################################

my $hbTxFile = "/Users/pgajer/projects/PECAN/data/Banfield_contax/Hug_Banfield.tx";
my %hbTx = read_tbl( $hbTxFile );

my $outFile = "/Users/pgajer/projects/PECAN/data/Banfield_contax/Hug_Banfield_spp.csv";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
#my $count = 0;
for my $id ( keys %hbTx )
{
  my $li = $hbTx{$id};
  ## find the the first pair of words (going from the end of $li) such that the
  ## first word starts from the capital letter and then second from the lower
  ## case letter
  my @f = split "_", $li;
  my $i = $#f;
  do
  {
    $i--;
  } while ( $i > 0 && ! ( $f[$i] =~ /^\p{Uppercase}/ && $f[$i+1] =~ /^\p{Lowercase}/ ) );
  $f[$i+1] =~ s/\.//;
  my $sp = $f[$i] . "_" . $f[$i+1];

  print OUT "$id,$sp,$li\n";
  #exit if $count == 100;
  #$count++;
}
close OUT;


print "\nOutput written to $outFile\n\n";

####################################################################
##                               SUBS
####################################################################

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my ($file, $skipHeader) = @_;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( $skipHeader )
  {
    my $header = <IN>;
  }
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
