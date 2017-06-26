#!/usr/bin/env perl

=head1 NAME

  setdiff.pl

=head1 DESCRIPTION



=head1 SYNOPSIS

  setdiff.pl <file1> <file2>

=head1 OPTIONS

=over

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /Users/pgajer/projects/PECAN/data/phylo_groups/v0.3/cx_hb_rdp_FL_5500_phGr_dir
  for d in *V3V4_dir; do echo $d; done | cut -f1 -d'_' | sed 's/phGr//g' > V3V4_dir.list
  ls *V3V4_dir/*V3V4_final.fa | cut -f1 -d'_' | sed 's/phGr//' > good_V3V4.ids

  setdiff.pl V3V4_dir.list good_V3V4.ids

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

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

if ( ! -e $file1 )
{
  warn "\n\n\tERROR: $file1 does not exist";
  print "\n\n";
  exit;
}

if ( ! -e $file2 )
{
  warn "\n\n\tERROR: $file2 does not exist";
  print "\n\n";
  exit;
}

my @a1 = read_array( $file1 );
my @a2 = read_array( $file2 );

my @d = diff( \@a1, \@a2 );
@d = sort @d;

print_array( \@d );


####################################################################
##                               SUBS
####################################################################

# print array to stdout
sub print_array
{
  my $a = shift;
  map {print "$_\n"} @{$a};
  #print "\n";
}


# difference of two arrays
sub diff
{
  my ($a1, $a2) = @_;

  my (%aa1, %aa2);

  foreach my $e (@{$a1}){ $aa1{$e} = 1; }
  foreach my $e (@{$a2}){ $aa2{$e} = 1; }

  my @d; # dfference array

  foreach my $e (keys %aa1, keys %aa2)
  {
    push @d, $e if exists $aa1{$e} && !exists $aa2{$e};
  }

  return @d;
}

# read table with one column
sub read_array
{
  my $file = shift;

  my @rows;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    push @rows, $_;
  }
  close IN;

  return @rows;
}


exit 0;
