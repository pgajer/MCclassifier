#!/usr/bin/env perl

=head1 NAME

  microcontax_to_banfield.pl

=head1 DESCRIPTION

  Align microcontax sequences to the Banfield tree of life (Bacteria and Archaea only).

  1. Start from Banfield's alignemtn of Archaea and Bacterial seq's.

  2. For a randomly selected genus of microcontax do align.seqs() with Banfield
  as a template.

  3. Merge the alignment and the template (make sure they have the same lengths)
  and iterate the process for each genus.


=head1 SYNOPSIS

  microcontax_to_banfield.pl [Options]

=head1 OPTIONS

=over

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

  microcontax_to_banfield.pl -o contax_trim_banfield_itr_mothur_algn2.fa

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Temp qw/ tempfile /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $nProc = 8;

GetOptions(
  "out-file|o=s"    => \my $outFile,
  "igs"             => \my $igs,
  "quiet"           => \my $quiet,
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

if ( !$outFile )
{
  print "ERROR: Missing output file name\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

# my $bfdFullAlgnFile = "/Users/pgajer/projects/PECAN/data/banfield/banfield_16S_algn_no_Ekta.fa";
# my $mxxFaFile       = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim.fa";
# my $mxxGeFaDir      = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_genus_dir";
# my $mxxTxFile       = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_genus.tx";

my $bfdFullAlgnFile = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_banfield_itr_mothur_algn.fa";
my $mxxFaFile       = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim.fa";
my $mxxGeFaDir      = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_genus_dir";
my $mxxTxFile       = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_genus_not_processed_yet.tx";

# my $bfdFullAlgnFile = "/Users/pgajer/projects/PECAN/data/banfield/banfield_16S_algn.fa";
# my $mxxFaFile       = "/Users/pgajer/projects/PECAN/data/microcontax/microcontax.trim.align_V3V4_algn_nr.fa";
# my $mxxGeFaDir      = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_V3V4_genus_dir";
# my $mxxTxFile       = "/Users/pgajer/projects/PECAN/data/microcontax/contax_trim_genus.tx";

my $mothur          = "mothur";

if ( $igs )
{
  $bfdFullAlgnFile = "/usr/local/projects/pgajer/projects/PECAN/data/banfield/banfield_16S_algn_no_Ekta.fa";
  $mxxFaFile       = "/usr/local/projects/pgajer/projects/PECAN/data/microcontax/contax_trim.fa";
  $mxxTxFile       = "/usr/local/projects/pgajer/projects/PECAN/data/microcontax/contax_trim_genus.tx";
  $mxxGeFaDir      = "/usr/local/projects/pgajer/projects/PECAN/data/microcontax/contax_trim_genus_dir";
}

if ( ! -e $bfdFullAlgnFile )
{
  warn "\n\n\tERROR: $bfdFullAlgnFile does not exist";
  print "\n\n";
  exit 1;
}

if ( ! -e $mxxFaFile )
{
  warn "\n\n\tERROR: $mxxFaFile does not exist";
  print "\n\n";
  exit 1;
}

if ( ! -e $mxxTxFile )
{
  warn "\n\n\tERROR: $mxxTxFile does not exist";
  print "\n\n";
  exit 1;
}

# creating output directory
# my $cmd = "mkdir -p $outDir";
# print "cmd=$cmd\n" if $dryRun;
# system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $tmpDir = "temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my $quietStr = "";
if ( $quiet )
{
  $quietStr = "--quiet";
}

my $debugStr = "";
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

my $nProcStr = "";
if ( $nProc )
{
  $nProcStr = "--thread $nProc";
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


## testig
# my $blen = extract_1st_seq_length_from_fasta_file( $bfdFullAlgnFile );
# print "\n\nBanfield's algn, $bfdFullAlgnFile, length: $blen\n\n"; exit;

print "--- Parsing microcontax contax_trim_genus.tx table\n";
my %geTbl = parse_tx_tbl( $mxxTxFile ); # ge => seq id of that genus


print "--- Aligning microcontax genera to Banfield alignment\n";

#my $algnFile = $outDir . "/microcontax_and_banfield_itr_mothur_algn.fa";
my $algnFile = $outFile;
my $cmd = "cp $bfdFullAlgnFile $algnFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

my $nGenera = keys %geTbl;
print "Number of genera to be aligned to Banfield alignment: $nGenera\n";

my $counter = 1;
for my $ge ( keys %geTbl )
{
  my $perc = sprintf("%.1f%%", 100 * $counter / $nGenera);
  print "\r$counter [$perc] Processing $ge                                           ";
  $counter++;
  my $candidateFile = $mxxGeFaDir . "/$ge" . ".fa";
  warn "\n\n\tERROR: $candidateFile not found" if ! -e $candidateFile;
  my ($seqCountBefore, $seqCountAfter) = mothur_align_and_add( $candidateFile, $algnFile, $nProc );
}

print "\r                                                           \n";

print "\n\n\tThe combined alignment written to $algnFile\n\n";

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

####################################################################
##                               SUBS
####################################################################


sub mothur_align_and_add
{
  my ($candidateFile, $templateFile, $nProc) = @_;

  my $seqCountBefore = seq_count( $templateFile );
  my $templateAlgnLen = get_algn_length( $templateFile );

  my @tmp;
  push (@tmp,"align.seqs(candidate=$candidateFile, template=$templateFile, flip=T, processors=$nProc)");

  print_array( \@tmp, "mothur commands" ) if ( $debug || $verbose );

  my $scriptFile = create_mothur_script( \@tmp );
  my $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  my $mothurAlgnFile = $candidateFile;
  $mothurAlgnFile =~ s/fa$/align/;
  if ( ! -e $mothurAlgnFile )
  {
    warn "\n\n\tERROR: Count not find $mothurAlgnFile";
    print "\n\n";
    exit 1;
  }

  my $mAlgnLen = get_algn_length( $mothurAlgnFile );

  if ( $mAlgnLen == $templateAlgnLen )
  {
    $cmd = "cat $mothurAlgnFile >> $templateFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }
  else
  {
    warn "\n\n\tERROR: mothur candidate alignment length, $mAlgnLen\nis not equal to the length of template alignment: $templateAlgnLen";
    print "\n\n";
    exit 1;
  }

  my $seqCountAfter = seq_count( $templateFile );
  if ( !$quiet )
  {
    print "\nAlignment file line count BEFORE: $seqCountBefore\n";
    print   "Alignment file line count AFTER:  $seqCountAfter\n\n";
  }

  return ($seqCountBefore, $seqCountAfter);
}

sub create_mothur_script
{
    my (@arr) = @{$_[0]};

    my ($fh, $inFile) = tempfile("mothur.XXXX", SUFFIX => '', OPEN => 1, DIR => $tmpDir);
    foreach my $c (@arr)
    {
        print $fh $c . "\n";
    }
    print $fh "quit()\n";
    close $fh;

    return $inFile;
}


sub parse_tx_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    push @{$tbl{$t}}, $id;
  }
  close IN;

  return %tbl;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
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

sub get_algn_length
{
  my $file = shift;

  my $length = extract_1st_seq_length_from_fasta_file( $file );

  return $length;
}

# extract sequence from a FASTA file
sub extract_1st_seq_length_from_fasta_file
{
    my($file) = @_;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";

    my $line;

    while ( ($line = <IN>) )
    {
	if ( $line !~ /^\s*$|^\s*\#|^>/ )
	{
	    last;
	}
    }

    chomp($line);
    my $seq = $line;
    foreach my $line (<IN>)
    {
	if ( $line =~ /^\s*$|^\s*\#|^>/ )
	{
	    last;
	}
	chomp $line;
	$seq .= $line;
    }

    close IN;

    return length($seq);
}

# extract sequence count of a fasta file
sub seq_count
{
  my $file = shift;

  my $wcline = qx/ grep -c '>' $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
}

exit 0;
