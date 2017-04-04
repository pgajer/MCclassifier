#!/usr/bin/env perl

=head1 NAME

  rm_incertae_sedis.pl

=head1 DESCRIPTION

  Remove "_incertae_sedis" from a given lineage file. move the input file to
  '_with_incertae_sedis.lineage' and move the updated one to the input file name.

=head1 SYNOPSIS

  rm_incertae_sedis.pl -i <input file>

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


  rm_incertae_sedis.pl -i Firmicutes_group_0_dir/Firmicutes_group_0.lineage

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

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
  print "\n\nERROR: Missing input file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -f $inFile )
{
  print "\n\nERROR: $inFile does not exist\n\n\n";
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

my @suffixes = (".lineage");
my ($filename, $dir, $suffix) = fileparse($inFile, @suffixes);
#my $inBasename = basename($inFile, @suffixes);

my $newFile = "$dir$filename" . "_with_incertae_sedis.lineage";

my $cmd = "mv $inFile $newFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "sed 's/_incertae_sedis//g' $newFile > $inFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "\n\n\tInput file moved to $newFile\n";
print "\tUpdated (no '_incertae_sedis') lineage file is in $inFile\n\n";

exit 0;
