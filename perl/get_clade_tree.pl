#!/usr/bin/env perl

=head1 NAME

  get_clade_tree.pl

=head1 DESCRIPTION

  Given a taxonomy file, a taxonomic rank and a tree, the script generates a
  tree of the clade of the taxonomic rank.

  This is useful when the taxonomic rank does not form a monophyletic clade.

=head1 SYNOPSIS

  get_clade_tree.pl -a <taxon file> -b <taxon> -t <tree file> -o <clade tree file> [Options]

=head1 OPTIONS

=over

=item B<--tx-file, -a>
  Taxon file.

=item B<--taxon, -b>
  Taxon.

=item B<--tree-file, -t>
  Tree file.

=item B<--output-file, -o>
  Output clade tree file.

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

  get_clade_tree.pl -a phylum.tx -b Tenericutes -t hb_FL.tree -o Tenericutes_clade.tree

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
  "tx-file|a=s"     => \my $txFile,
  "taxon|b=s"       => \my $taxon,
  "tree-file|t=s"   => \my $treeFile,
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

if ( !$txFile )
{
  print "\n\n\tERROR: Missing taxon file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outFile )
{
  print "\n\n\tERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$treeFile )
{
  print "\n\n\tERROR: Missing tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$taxon )
{
  print "\n\n\tERROR: Missing taxon\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $txFile )
{
  warn "\n\n\tERROR: $txFile does not exist";
  print "\n\n";
  exit;
}

if ( ! -e $treeFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n\n";
  exit;
}

####################################################################
##                               MAIN
####################################################################

get_clade_tree( $treeFile, $txFile, $taxon, $outFile );

####################################################################
##                               SUBS
####################################################################

sub get_clade_tree
{
  my ($treeFile, $txFile, $taxon, $outFile) = @_;

exit 0;
