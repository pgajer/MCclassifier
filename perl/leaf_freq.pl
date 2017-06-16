#!/usr/bin/env perl

=head1 NAME

  leaf_freq.pl

=head1 DESCRIPTION

  This script generates a table of leaf label frequencies from a phylo tree in the newick format.

=head1 SYNOPSIS

  leaf_freq.pl -t <tree file> -o <freq file> [Options]

=head1 OPTIONS

=over

=item B<--tree-file, -t>
  Tree file.

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

  leaf_freq.pl --print -i cxhb_FL_QCed_trAlgn_condensed_phylum.tree -o cxhb_FL_QCed_trAlgn_condensed_phylum.freq

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
    "tree-file|t=s"   => \my $treeFile,
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

if ( !$treeFile )
{
  print "\n\n\tERROR: Missing tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outFile )
{
  print "\n\n\tERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $treeFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n\n";
  exit;
}

my $tmpDir = "leaf_freq_temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

####################################################################
##                               MAIN
####################################################################

my @leaves = get_leaves( $treeFile );

my %freq;
map { $freq{$_}++ } @leaves;

if ( $print )
{
  print "\nLeaf Frequencies\n";
  print_formated_freq_tbl( \%freq );
  print "\n";
}

write_formated_freq_tbl( \%freq, $outFile );

remove_tmp_dir();

####################################################################
##                               SUBS
####################################################################

sub remove_tmp_dir
{
  my $cmd = "rm -rf $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

sub get_leaves
{
  my $treeFile = shift;

  my ($fh, $leavesFile) = tempfile("leaves.XXXX", SUFFIX => '', OPEN => 1, DIR => $tmpDir);
  close $fh;

  my $cmd = "nw_labels -I $treeFile > $leavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @leaves = read_array( $leavesFile );

  return @leaves;
}

# read table with one column
sub read_array
{
  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -e $file )
  {
    warn "\n\n\nERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( defined $hasHeader )
  {
    <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    push @rows, $_;
  }
  close IN;

  return @rows;
}

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
    print "$_$pad" . $rTbl->{$_} . "\n";
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
