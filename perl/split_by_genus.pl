#!/usr/bin/env perl

=head1 NAME

  split_by_genus.pl

=head1 DESCRIPTION

  Given CXHG and subRDP lineage and fasta files + CXHG genus tree, generate
  - outgroup sequence for each genus
  - lineage file of CXHG+subRDP sequences of that genus + OG sequence
  - fasta file of CXHG+subRDP sequences of that genus + OG sequence

=head1 SYNOPSIS

  split_by_genus.pl -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--output-dir, -o>
  Output directory.

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

  split_by_genus.pl -o cxhb_subRDP_FL_dir

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
  "output-dir|o=s"  => \my $outDir,
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

if ( !$outDir )
{
  print "\n\n\tERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

# creating output directory
my $cmd = "mkdir -p $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $tmpDir = $outDir . "/temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my $debugStr = "";
my $quietStr = "--quiet";
if ( $debug )
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ( $verbose )
{
  $verboseStr = "--verbose";
}


####################################################################
##                               MAIN
####################################################################


my $cxhgFaFile = "";


my $ourFaFile = "";

my $treeFile = "";

my %geOGtbl = get_genus_og_tbl( "" );

print "--- Parsing CXHB lineage table\n";
my $cxhgLiFile = "";
my %cxhbGeTbl = get_genus_tbl( $cxhgLiFile );

print "--- Parsing subRDP lineage table\n";
my $ourLiFile = "";
my %ourGeTbl  = get_genus_tbl( $ourLiFile );

# If we didn't have to add OG seq's to each genus ref seq's, we could simply
# concatenate CXHB and subRDP genus taxonomy tables (that would have to be
# extracted from the lineage tables) and then use multi_select_seq/tx to create
# fasta/lineage files for each genus.

# The followed loop creates a genus table with the only modlification of OG
# sequence being added to each genus.

print "--- Generating genus table\n";
my $geTblFile = $outDir . "/genus_with_OG.tx";
open OUT, ">$geTblFile" or die "Cannot open $geTblFile for writing: $OS_ERROR\n";
for my $ge ( keys %cxhbGeTbl  )
{
    my @ids = @{$cxhbGeTbl{$ge}};
    if ( exists $ourGeTbl{$ge} )
    {
        push @ids, @{$ourGeTbl{$ge}};
    }
    # adding OG
    if ( exists $geOGtbl{$ge} )
    {
        push @ids, $geOGtbl{$ge};
    }
    else
    {
        warn("\n\n\tERROR: $ge does not have have OG sequence");
        print("\n\n");
        exit 1;
    }

    for my $id ( @ids )
    {
        print OUT "$id\t$ge\n";
    }
}
close OUT;


####################################################################
##                               SUBS
####################################################################


exit 0;
