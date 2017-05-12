#!/usr/bin/env perl

=head1 NAME

  get_lineage.pl

=head1 DESCRIPTION

  Given a taxonomy file (for example, PECAN output columns 1 & 2), and a source lineage file, 
  use the taxon to search for the matching full taxonomic lineage, and print out a new lineage 
  file with sequence ID in column 1, and lineage in column 2.
  
=head1 SYNOPSIS

  get_lineage.pl -a <taxonomic annotation file> -s <source lineage file> [Options]

=head1 OPTIONS

=over

=item B<--taxonomic-annotation file, -a>
  List of selected outgroup sequence IDs. 

=item B<--source-lineage-file, -s>
  The starting, source file containing lineage information for all sequences.

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
  "annotation-file|a=s"   => \my $tx,
  "source-lineage|s=s"    => \my $sourceLineage,
  "igs"                => \my $igs,
  "verbose|v"            => \my $verbose,
  "debug"               => \my $debug,
  "debug2"              => \my $debug2,## file name debug
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if ( !$tx )
{
  print "\n\nERROR: Missing taxonomic annotation file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$sourceLineage)
{
  print "\n\nERROR: Missing source lineage file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
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

#my @suffixes = (".seqID",".txt");
#my $Seq = basename($SeqID, ".seqID");

## From a list of sequences get the full taxonomic lineage from the source file.
my @source;
open (SOURCE, "<$sourceLineage") or die "Cannot open $sourceLineage for reading: $OS_ERROR\n";
while (<SOURCE>)
{
    @source = split /[\t;]/, $_;
}
close SOURCE;

my @outLineage;
open (IN, "<$tx") or die "Cannot open $tx for reading: $OS_ERROR\n";
while (<IN>)
{
    my @t = split /[\t]/, $_;
    
    #if (exists ($source{$t[1]}) )
	for ($t[1] =~ $source[8])
    {
		push @outLineage, $t[0]."\t".$source[1].";".$source[2].";".$source[3].";".$source[4].";".$source[5].";".$source[6].";".$source[7].";".$t[1]."\n";
		#print "$_ found in $sourceLineage and lineage is: ".$source{$_}."\n";
		}
}
close IN;


my $newLineage = "new.lineage";
open (OUT, ">$newLineage") or die "Cannot open $newLineage for reading: $OS_ERROR\n";
print OUT @outLineage;
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

print "Lineages for sequences in $tx written to $newLineage.\n\n";



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
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

exit;
