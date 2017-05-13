#! /usr/bin/perl

=head1 NAME

  uc2clstr2.pl

=head1 DESCRIPTION

  Produce a clstr2 file from a uc file

=head1 SYNOPSIS

  uc2clstr2.pl -i <uc file> -o <clstr2 output file> [Options]

=head1 OPTIONS


=over

=item B<--input-file, -i>
  Input uc file.

=item B<--output-file, -o>
  Output clstr2 file.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /usr/local/projects/pgajer/projects/HMP/RAV151_100.0

  uc2clstr2.pl -i seqs.uc -o seqs.clstr2

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Temp qw/ :POSIX /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-file|i=s"  => \my $inFile,
  "output-file|o=s" => \my $outFile,
  "igs"             => \my $igs,
  "dry-run"         => \my $dryRun,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$inFile)
{
  print "ERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outFile)
{
  print "ERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my $usearch4        = "usearch4";
my $clstr_to_clstr2 = "clstr_to_clstr2.pl";

if ( defined $igs )
{
  $usearch4        = "/usr/local/projects/pgajer/bin/usearch4";
  $clstr_to_clstr2 = "/home/pgajer/devel/MCclassifier/perl/clstr_to_clstr2.pl";
}

####################################################################
##                               MAIN
####################################################################

if ( $inFile !~ /uc$/)
{
  print "ERROR: Input file should have uc suffix!\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my $startRun = time();

print "--- Generating clstr file\n";
my $file = tmpnam();
my $cmd = "$usearch4 --uc2clstr $inFile --output $file";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating clstr2 file\n";
$cmd = "$clstr_to_clstr2 -i $file -o $outFile; rm -f $file";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $endRun = time();
my $runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  my $timeMin = int($runTime / 60);
  my $timeSec = $runTime % 60;
  print "Completed in $timeMin:$timeSec\n"
}
else
{
  print "Completed in $runTime seconds\n"
}

exit;
