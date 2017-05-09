#!/usr/bin/env perl

=head1 NAME

  rename_pecan_vaginal_spp.pl

=head1 DESCRIPTION

  Change the names of some species for vaginal samples

=head1 SYNOPSIS

  rename_pecan_vaginal_spp.pl -i <input taxon file> -o <output taxon file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

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

  cd /Users/pgajer/projects/M_and_M/new_16S_classification_data

  rename_pecan_vaginal_spp.pl -i mm_phylo_clfy_dir/query.tx -o mm_phylo_clfy_dir/query_chgd.tx

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
  "input-file|i=s"  => \my $inFile,
  "output-file|o=s" => \my $outFile,
  "igs"             => \my $igs,
  "verbose|v"       => \my $verbose,
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
  print "\n\n\tERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$outFile)
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

my $trTbl = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/vag_pecan_tr_tbl.txt";

if ($igs)
{
  $trTbl = "/usr/local/projects/pgajer/devel/MCextras/data/vag_pecan_tr_tbl.txt";
}

####################################################################
##                               MAIN
####################################################################

print "--- Parsing translation table\n";
my %trTbl = readTbl( $trTbl );

my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

print "--- Checking the number of lines in the taxonomy file\n";
my $wcline = qx/ wc -l $inFile /;
$wcline =~ s/^\s+//;
my ($nLines, $str) = split /\s+/, $wcline;
print "nLines: $nLines\n";

print "--- Translating the taxonomy table\n";
open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
my $counter = 1;
for (<IN>)
{
  if ($counter % 500 == 0)
  {
    $endRun = time();
    $runTime = $endRun - $startRun;
    if ( $runTime > 60 )
    {
      $timeMin = int($runTime / 60);
      $timeSec = sprintf("%02d", $runTime % 60);
      $timeStr = "$timeMin:$timeSec";
    }
    else
    {
      $runTime = sprintf("%02d", $runTime);
      $timeStr = "$timeMin:$runTime";
    }
    amy $perc = sprintf("%.1f%%", 100 * $counter / $nLines);
    print "\r$timeStr [$perc]";
  }
  $counter++;

  chomp;
  my ($id, $tx) = split /\s+/, $_;
  if ( exists $trTbl{$tx} )
  {
    print OUT "$id\t" . $trTbl{$tx} . "\n";
  }
  else
  {
    print OUT "$_\n";
  }
}
close OUT;
close IN;


####################################################################
##                               SUBS
####################################################################

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
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
