#!/usr/bin/env perl

=head1 NAME

  phylip2fasta.pl

=head1 DESCRIPTION

  Translation of PHYLIP to FASTA format.

=head1 SYNOPSIS

  phylip2fasta.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input interleaved PHYLIP alignment file.

=item B<--output-file, -o>
  Output FASTA file.

=item B<--interleaved, -j>
  Set to 0 if the input file is in the sequencial PHYLIP format.

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

  phylip2fasta.pl -i Banfield_medoids_FL_algn_i.phy.reduced -o Banfield_medoids_FL_algn_reduced.fa

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Bio::AlignIO;
use Bio::SimpleAlign;

# http://doc.bioperl.org/bioperl-live/Bio/AlignIO.html
# http://doc.bioperl.org/bioperl-live/Bio/AlignIO/phylip.html
# http://www.bioperl.org/wiki/PHYLIP_multiple_alignment_format

# Taken from
# https://stackoverflow.com/questions/15397672/how-to-convert-phylip-format-to-fasta

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $interleaved = 1;

GetOptions(
  "input-file|i=s"  => \my $inFile,
  "output-file|o=s" => \my $outFile,
  "interleaved|j=i" => \$interleaved,
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

my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

my $in = Bio::AlignIO->new(
    -file => $inFile,
    -format => 'phylip',
    -interleaved => $interleaved,
    -longid => 1
    );

## “-idlength” is for writing, “-longid” is for reading.

my $out = Bio::AlignIO->new(-file   => ">$outFile",
			    -format => 'fasta');

while ( my $aln = $in->next_aln() )
{
  $out->write_aln($aln);
}

## report timing
$endRun = time();
$runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  $timeMin = int($runTime / 60);
  $timeSec = sprintf("%02d", $runTime % 60);
  print "\rCompleted in $timeMin:$timeSec\n"
}
else
{
  print "\rCompleted in $runTime seconds\n"
}


print "Output written to $outFile\n";


exit 0;
