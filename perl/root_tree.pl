#!/usr/bin/env perl

=head1 NAME

  root_tree.pl

=head1 DESCRIPTION

  Given a phylogenetic tree in the Newick format, root it by changing )0.990:0.00612; to ); at the end of the tree string.
  The change is done in place !

=head1 SYNOPSIS

  root_tree.pl -i <input file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Firmicutes_dir/Firmicutes_group_6_dir

  root_tree.pl -i Firmicutes_group_6_pruned_Lactobacillus_Pediococcus.tree

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
  "verbatim|v"      => \my $verbatim,
  "debug"           => \my $debug,
  "dry-run"         => \my $dryRun,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$inFile)
{
  warn "\n\n\tERROR: Missing input file";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $inFile )
{
  warn "\n\n\tERROR: $inFile does not exist";
  print "\n\n";
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
my $tree = <IN>;
close IN;

$tree =~ s/\d\.\d+:\d\.\d+;$/;/;

my $outFile = $inFile;
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
print OUT $tree;
close OUT;

exit 0;
