#!/usr/bin/env perl

=head1 NAME

  MED_node_rep_rename.pl

=head1 DESCRIPTION

  Given PECAN output for MED node representatives and a count table file, 
  add PECAN taxon names to MED node representative labels in a count table.
  
=head1 SYNOPSIS

  MED_node_rep_rename.pl -a <taxonomic annotation file> -c <count table file> []

=head1 OPTIONS

=over

=item B<--taxonomic-annotation file, -a>
  List of selected outgroup sequence IDs. 

=item B<--count-table-file, -c>
  The starting, source file containing lineage information for all sequences.
 
=item B<--output file, -o>
 Out file name.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE


  
=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper qw(Dumper);
use File::Basename;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "taxonomy-file|a=s"   => \my $txFile,
  "count-tbl|c=s"    => \my $countTbl,
  "out-file|o=s"          => \my $outFile,
  "igs"                   => \my $igs,
  "verbose|v"             => \my $verbose,
  "debug"                 => \my $debug,
  "debug2"                => \my $debug2,## file name debug
  "dry-run"               => \my $dryRun,
  "help|h!"               => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if ( !$txFile )
{
  print "\n\nERROR: Missing taxonomic annotation file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$countTbl)
{
  print "\n\nERROR: Missing source lineage file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
if (!$outFile)
{
    print "\n\nOutput file name not provided";
    print "\n Writing to new.lineage";
    $outFile = "new.lineage";
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

#my @suffixes = (".seqID",".txt");
#my $Seq = basename($SeqID, ".seqID");

## From a list of sequences get the full taxonomic lineage from the source file.
my $taxa;
my $lineage;

my %source = read2colTbl($sourceLineage);
my %tx = read2colTbl($txFile);
my %done;

my $newLineage = $outFile;
open (OUT, ">$newLineage") or die "Cannot open $newLineage for reading: $OS_ERROR\n";

for my $seqID (keys %tx) 
{
  my $taxa = $tx{$seqID};
  my $lineage = $source{$taxa};
  
  next if exists ($done{$taxa});

  $done{$lineage} = $taxa;

  print(OUT "$seqID\t$lineage\n");

}

close OUT;
#print "\nOutgroup lineages written to $seqLineage.\n";


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

print "\nLineages for ". scalar(keys %tx) . " sequences from $txFile written to $newLineage.\n\n";



####################################################################
##                               SUBS
####################################################################


sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR: $file does not exist\n\n\n";
    exit;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_,;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# read two column table
sub read2colTbl{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

exit;
