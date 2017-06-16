#!/usr/bin/env perl

=head1 NAME

  taxon_from_lineage.pl

=head1 DESCRIPTION

  This script generates a taxon file from a lineage file.

=head1 SYNOPSIS

  taxon_from_lineage.pl -l <lineage file> -t <taxon> -o <taxon file> [Options]

=head1 OPTIONS

=over

=item B<--lineage-file, -l>
  Lineage file.

=item B<--taxon, -t>
  Taxonomic rank. Possible values: phylum, class, order, family, genus, species.

=item B<--output-file, -o>
  Output taxon file.

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

  taxon_from_lineage.pl -l cxhb_FL_QCed_trAlgn.lineage -t phylum -o cxhb_FL_QCed_trAlgn_phylum.tx

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
    "lineage-file|l=s" => \my $liFile,
    "taxon|t=s"        => \my $taxon,
    "output-file|o=s"  => \my $outFile,
    "verbose|v"        => \my $verbose,
    "debug"            => \my $debug,
    "dry-run"          => \my $dryRun,
    "help|h!"          => \my $help,
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
elsif ( !$outFile )
{
  print "\n\n\tERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$taxon )
{
  print "\n\n\tERROR: Missing taxon\n\n";
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

open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
open IN, "<$liFile" or die "Cannot open $liFile for reading: $OS_ERROR";

if ( $taxon eq "phylum" )
{
    while (<IN>)
    {
        chomp;
        my ($id, $li) = split /\s+/;
        my @f = split ";", $li;
        shift @f; # Root
        shift @f; # domain
        my $ph = shift @f;
        print OUT "$id\t$ph\n";
    }
}
elsif ( $taxon eq "class" )
{
    while (<IN>)
    {
        chomp;
        my ($id, $li) = split /\s+/;
        my @f = split ";", $li;
        shift @f; # Root
        shift @f; # domain
        shift @f; # phylum
        my $cl = shift @f;
        print OUT "$id\t$cl\n";
    }
}
elsif ( $taxon eq "order" )
{
    while (<IN>)
    {
        chomp;
        my ($id, $li) = split /\s+/;
        my @f = split ";", $li;
        shift @f; # Root
        shift @f; # domain
        shift @f; # phylum
        shift @f; # class
        my $or = shift @f;
        print OUT "$id\t$or\n";
    }
}
elsif ( $taxon eq "family" )
{
    while (<IN>)
    {
        chomp;
        my ($id, $li) = split /\s+/;
        my @f = split ";", $li;
        pop @f; # species
        pop @f; # genua
        my $fa = pop @f;
        print OUT "$id\t$fa\n";
    }
}
elsif ( $taxon eq "genus" )
{
    while (<IN>)
    {
        chomp;
        my ($id, $li) = split /\s+/;
        my @f = split ";", $li;
        my $sp = pop @f;
        my $ge = pop @f;
        print OUT "$id\t$ge\n";
    }
}
elsif ( $taxon eq "species" )
{
    while (<IN>)
    {
        chomp;
        my ($id, $li) = split /\s+/;
        my @f = split ";", $li;
        my $sp = pop @f;
        print OUT "$id\t$sp\n";
    }
}
else
{
    warn "\n\n\tERROR: unrecognized taxon string $taxon";
    print "\n\n";
    exit 1;
}
close IN;
close OUT;


####################################################################
##                               SUBS
####################################################################


exit 0;
