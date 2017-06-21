#!/usr/bin/env perl

=head1 NAME

  clean_master_lineage_v2.pl

=head1 DESCRIPTION

  This script modifies species names so that they are of the form Genus_speciesname
  where Genus is capitilized and speciesname is a lowercase

  '/', '-' and '_' characters (except the '_' between genus and speciesname) are changed to '|'

  _gp1/2/3 changes to |gp1/2/3

  Escherichia-Shigella_sp => Escherichia|Shigella_sp

  sp.suffix => sp

  Genus name and prefix of species name have to be the same

=head1 SYNOPSIS

  clean_master_lineage_v2.pl -i <lineage file> -o <output lineage file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input lineage file.

=item B<--output-file, -o>
  Output lineage file.

=item B<--verbose, -v>
  Prints content of some output files. Default value: 5000.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/RDP/r2.2

  clean_master_lineage_v2.pl -i rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis_nr.lineage -o rdp_Bacteria_curated.lineage

=cut

use strict;
use warnings;
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
  "input-file|i=s"      => \my $inFile,
  "output-file|o=s"      => \my $outFile,
  "quiet"               => \my $quiet,
  "igs"                 => \my $igs,
  "verbose|v"           => \my $verbose,
  "dry-run"             => \my $dryRun,
  "debug"               => \my $debug,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$inFile )
{
  print "ERROR: Missing input lineage file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outFile )
{
  print "ERROR: Missing output file name\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
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

my $suffix;

open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
for ( <IN> )
{
  chomp;
  my ($id, $li) = split /\s+/;

  my @f = split ";", $li;
  my $sp = pop @f;
  my $ge = pop @f;

  my $origSp = $sp;
  #$sp = "AAA-BBB_gp1_sp.sldur-23X";
  #$sp = "Curtobacterium_MG-2011-84-GV";
  #$sp = "Mycobacterium_n";
  #print "sp BEFORE: $sp\n";
  $sp =~ s/\//\|/g;
  $sp =~ s/\-/\|/g;
  $sp =~ s/_gp/|gp/;
  $sp =~ s/\..+//;
  $sp =~ s/_[[:upper:]].+/_sp/;
  $sp =~ s/_\w$/_sp/;
  #print "sp AFTER: $sp\n";  exit 1;

  if ( $sp !~ /[[:upper:]][\||\w]+_[[:lower:]]\w+/ && $sp !~ /BVAB/ )
  {
    warn "\n\n\tERROR: incorrect format sp: $sp";
    print "origSp: $origSp\n\n";
    exit 1;
  }

  ($ge, $suffix) = split "_", $sp;
  $li = join ";", (@f, $ge, $sp);
  print OUT "$id\t$li\n";
  #print "$id\t$li\n"; exit 1;
}
close IN;
close OUT;

exit 0;
