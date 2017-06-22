#!/usr/bin/env perl

=head1 NAME

  build_tree.pl

=head1 DESCRIPTION

  This script

  - builds an alignment
  - trims it
  - dereplicates trimmed sequences
  - builds a tree using FastTree
  - reroots the tree

  Note: This used to be called build_tree.pl.

=head1 SYNOPSIS

  build_tree.pl -i <phylo-group name> -p <number of processors to use for ginsi>

=head1 OPTIONS

=over

=item B<--input-group-name, -i>
  Prefix of input group. For example, if Firmicutes_group_6_dir is a directory
  Firmicutes_group_6 group, then the input group name is "Firmicutes_group_6".

=item B<--num-proc, -p>
  Number of processors to be used. Default value: 8.

=item B<--mafft-auto-algn>
  Generate MAFFT alignment using --auto flag. For fasta files > 200 seq's (?) it will run

  FFT-NS-2 (Fast but rough)
  Progressive method (guide trees were built 2 times.)

=item B<--linsi-algn>
  Generate MAFFT L-INS-i alignment.

  L-INS-i (Probably most accurate, very slow)
  Iterative refinement method (<16) with LOCAL pairwise alignment information

=item B<--ginsi-algn>
  Generate MAFFT G-INS-i alignment.

  G-INS-i (Suitable for sequences of similar lengths, very slow)
  Iterative refinement method (<16) with GLOBAL pairwise alignment information

=item B<--mothur-algn>
  Generate mothur aligner is used with mothur's SILVA SEED db.

=item B<--criteria, -c>
  An integer between 0 and 100. The alignment will be trimmed at the positions
  where percentage of sequences that start and end is at least the criteria value.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--quiet>
  Do not print progress messages.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  cd /Users/pgajer/projects/PECAN/data/phylo_groups/v0.3/cx_hb_ssubRDP_FL_5k_phGr_dir

  build_tree.pl --debug --mafft-auto-algn -i phGr50 -p 4

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

my $criteria = 95;

GetOptions(
  "input-group|i=s" 	=> \my $grPrefix,
  "num-proc|p=i"        => \my $nProc,
  "mafft-auto-algn"     => \my $runMAFFTauto,
  "linsi-algn"          => \my $runLinsi,
  "ginsi-algn"          => \my $runGinsi,
  "mothur-algn"         => \my $runMothurAlgn,
  "criteria|c=i"        => \$criteria,
  "igs"                 => \my $igs,
  "johanna"             => \my $johanna,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "debug2"              => \my $debug2,## file name debug
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$grPrefix )
{
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$runMAFFTauto && !$runLinsi && !$runGinsi && !$runMothurAlgn )
{
  warn "\n\n\tERROR: Missing alignment type";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}


my $johannaStr = "";
if ( defined $johanna )
{
  $johannaStr = "--johanna";
}

my @suffixes;

my $debugStr = "";
my $quietStr = "--quiet";
if ($debug)
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ($verbose)
{
  $verboseStr = "--verbose";
}

my $nProcStr = "";
if ($nProc)
{
  $nProcStr = "--thread $nProc";
}

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "\n\n\tERROR: $grDir does not exist";
  print "\n\n";
  exit 1;
}

my @p = split "/", $grPrefix;
$grPrefix = pop @p;

my $tmpDir                = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.3/tmpDir";
my $silvaSEED             = "/Users/pgajer/projects/PECAN/data/SILVA/v128/silva_seed_v128_algn.fa";

my $nw_labels             = "nw_labels";
my $nw_order              = "nw_order";
my $nw_condense           = "nw_condense";
my $nw_rename             = "nw_rename";
my $nw_prune              = "nw_prune";
my $nw_reroot             = "nw_reroot";
my $uc2clstr2             = "uc2clstr2.pl";
my $extract_seq_IDs       = "extract_seq_IDs.pl";
my $select_seqs           = "select_seqs.pl";
my $rmGaps                = "rmGaps";
my $ginsi                 = "/usr/local/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)
my $linsi                 = "/usr/local/bin/mafft --maxiterate 1000 --localpair"; # MAFFT v7.310 (2017/Mar/17)
my $mafft                 = "/usr/local/bin/mafft --auto"; # MAFFT v7.310 (2017/Mar/17)
my $FastTree              = "FastTree";
my $R                     = "R";
my $fix_fasta_headers     = "fix_fasta_headers.pl";
my $mothur                = "/Users/pgajer/bin/mothur";
my $usearch6              = "/Users/pgajer/bin/usearch6.0.203_i86osx32";
my $vicut                 = "/Users/pgajer/devel/vicut/bin/vicut";
my $readNewickFile        = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";
my $trim_align            = "trim_align.pl";

my $vsearchSORT;
my $vsearch;

my $igsStr = "";
if ( defined $igs )
{
  $tmpDir                = "/home/pgajer/tmpDir";
  $silvaSEED             = "/home/pgajer/projects/PECAN/data/SILVA/v128/silva_seed_v128_algn.fa";

  $fix_fasta_headers     = "/home/pgajer/devel/MCclassifier/perl/fix_fasta_headers.pl";
  $nw_labels             = "/usr/local/projects/pgajer/bin/nw_labels";
  $nw_order              = "/usr/local/projects/pgajer/bin/nw_order";
  $nw_condense           = "/usr/local/projects/pgajer/bin/nw_condense";
  $nw_rename             = "/usr/local/projects/pgajer/bin/nw_rename";
  $nw_prune              = "/usr/local/projects/pgajer/bin/nw_prune";
  $nw_reroot             = "/usr/local/projects/pgajer/bin/nw_reroot";
  $trim_align            = "/home/pgajer/devel/MCclassifier/perl/trim_align.pl";
  $uc2clstr2             = "/home/pgajer/devel/MCclassifier/perl/uc2clstr2.pl";
  $extract_seq_IDs       = "/home/pgajer/devel/MCclassifier/perl/extract_seq_IDs.pl";
  $select_seqs           = "/home/pgajer/devel/MCclassifier/perl/select_seqs.pl";
  $rmGaps                = "/usr/local/projects/pgajer/bin/rmGaps";
  $ginsi                 = "/home/pgajer/bin/mafft --maxiterate 1000 --globalpair"; # MAFFT v7.310 (2017/Mar/17)
  $linsi                 = "/home/pgajer/bin/mafft --maxiterate 1000 --localpair"; # MAFFT v7.310 (2017/Mar/17)
  $mafft                 = "/home/pgajer/bin/mafft --auto"; # MAFFT v7.310 (2017/Mar/17)
  $FastTree              = "/home/pgajer/bin/FastTree_no_openMP";
  $R                     = "/home/pgajer/bin/R";
  $mothur                = "/usr/local/projects/pgajer/bin/mothur";
  $usearch6              = "/local/projects/pgajer/bin/usearch6.0.203_i86linux32";
  $vicut                 = "/usr/local/projects/pgajer/bin/vicut";
  $readNewickFile        = "/local/projects/pgajer/devel/MCclassifier/perl/read.newick.R";
  $vsearchSORT           = "/usr/local/packages/vsearch/bin/vsearch";
  $vsearch               = "/usr/local/bin/vsearch";

  $quietStr              = "";
  $igsStr                = "--igs";
}

####################################################################
##                               MAIN
####################################################################

chdir $grDir;
print "--- Changed dir to $grDir\n";



my $faFile	     = $grPrefix . ".fa";
my $algnFile	 = $grPrefix . "_algn.fa";
my $trAlgnFile   = $grPrefix . "_algn_trimmed.fa";
my $treeFile	 = $grPrefix . ".tree";
my $outgroupFile = $grPrefix . "_outgroup.seqIDs";

## Gathering outgroup data
print "--- Parsing $outgroupFile\n";
my @ogSeqIDs = read_array($outgroupFile);
my %ogInd = map{$_ =>1} @ogSeqIDs; # outgroup elements indicator table

print "--- Extracting seq IDs from alignment fasta file\n";
my @seqIDs = get_seqIDs_from_fa($faFile);

print "--- Testing if outgroup sequences are part of seqIDs\n";
my @ogDiff = diff( \@ogSeqIDs, \@seqIDs );
@ogSeqIDs = comm( \@ogSeqIDs, \@seqIDs );

if ( scalar(@ogDiff) != 0 )
{
  warn "\n\tERROR the following outgroup seq IDs are not in the trimmed alignment file:\n\n";
  print_array( \@ogDiff );
  exit 1;
}

if ( scalar( @ogSeqIDs ) == 0 )
{
  warn "\n\n\tERROR: All outgroup seq's were lost";
  print "\n\n";
  exit 1;
}

if ( $debug )
{
  print "\n\tNumber of seq's in the trimmed alignment/lineage files: " . @seqIDs . "\n";
  print   "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n\n";
}

if ( $runMAFFTauto )
{
  print "--- Generating MAFFT --auto alignment\n";
  my $mafftAlgnFile = abs_path( $grPrefix . "_mafft_auto_algn.fa" );
  my $cmd = "rm -f $mafftAlgnFile; time $mafft --inputorder $quietStr $nProcStr $faFile > $mafftAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $algnFile; ln -s $mafftAlgnFile $algnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
elsif ( $runLinsi )
{
  print "--- Generating MAFFT L-INS-i alignment\n";
  my $linsiAlgnFile = abs_path( $grPrefix . "_linsi_algn.fa" );
  my $cmd = "rm -f $algnFile; time $linsi --inputorder $quietStr $nProcStr $faFile > $linsiAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $algnFile; ln -s $linsiAlgnFile $algnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
elsif ( $runGinsi )
{
  print "--- Generating MAFFT G-INS-i alignment\n";
  my $ginsiAlgnFile = abs_path( $grPrefix . "_ginsi_algn.fa" );
  my $cmd = "rm -f $algnFile; time $ginsi --inputorder $quietStr $nProcStr $faFile > $ginsiAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $algnFile; ln -s $ginsiAlgnFile $algnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
elsif ( $runMothurAlgn )
{
  print "--- Generating mothur alignment using mothur's SILVA SEED db\n";
  my $mothurAlgnFile = abs_path( $grPrefix . "_mothur_algn.fa" );
  mothur_align( $faFile, $nProc, $mothurAlgnFile );

  my $cmd = "rm -f $algnFile; ln -s $mothurAlgnFile $algnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

my $urTreeFile = $grPrefix . "_unrooted.tree";

if ( ! $runMothurAlgn ) # mothur aligner removes bases and so we cannot use its
                        # result to then filter out redundant sequnences (after
                        # trimming)
{
  print "--- Generating sequence summary of the alignment file\n";
  my @tmp;
  push (@tmp,"summary.seqs(fasta=$algnFile)");
  print_array(\@tmp, "mothur commands") if ($debug || $verbose);

  my $scriptFile = create_mothur_script( \@tmp );
  my $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  print "--- Trimming alignment to length of $criteria% of sequences\n";
  my $summaryFile     = $grPrefix . "_algn.summary";
  if ( ! -e $summaryFile )
  {
    warn "\n\n\tERROR: cannot find $summaryFile";
    print "\n\n";
    exit 1;
  }
  $cmd = "$trim_align -c $criteria -j $summaryFile -i $algnFile -o $trAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Dereplicating trimmed alignemnt file\n";
  my $trFaFile = $grPrefix . "_trimmed.fa";
  $cmd = "$rmGaps -i $trAlgnFile -o $trFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $ucFile = $grPrefix . "_trimmed.uc";
  my $nrFile = $grPrefix . "_trimmed_nr.fa";
  $cmd = "$usearch6 -cluster_fast $trFaFile -id 1.0 -uc $ucFile -centroids $nrFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  print "--- Restricting alignment to non-redundant seq's\n";
  my $nrSeqsFile = $grPrefix . "_trimmed_nr.seqIDs";
  $cmd = "$extract_seq_IDs -i $nrFile -o $nrSeqsFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $nrTrAlgnFile   = $grPrefix . "_algn_trimmed_nr.fa";
  $cmd = "$select_seqs -s $nrSeqsFile -i $trAlgnFile -o $nrTrAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "mv $nrTrAlgnFile $trAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Producing tree\n";
  $cmd = "rm -f $urTreeFile; $FastTree -nt $trAlgnFile > $urTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;# && !$skipFastTree;
}
else
{
  print "--- Producing tree\n";
  my $cmd = "rm -f $urTreeFile; $FastTree -nt $algnFile > $urTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;# && !$skipFastTree;
}

print "--- Rooting the tree\n";
my $cmd = "rm -f $treeFile; $nw_reroot $urTreeFile @ogSeqIDs | $nw_order -  > $treeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


####################################################################
##                               SUBS
####################################################################

# extract line count of a file
sub line_count
{
  my $file = shift;

  my $wcline = qx/ wc -l $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
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

sub mothur_align
{
  my ($candidateFile, $nProc, $algnFile ) = @_;

  my @tmp;
  push (@tmp,"align.seqs(candidate=$candidateFile, template=$silvaSEED, flip=T, processors=$nProc)");

  my $mothurAlgnFile = $candidateFile;
  $mothurAlgnFile =~ s/fa$/align/;

  push (@tmp,"filter.seqs(fasta=$mothurAlgnFile)");

  my $mothurFilteredFile = $candidateFile;
  $mothurFilteredFile =~ s/fa$/filter.fasta/;

  print_array( \@tmp, "mothur commands" ) if ( $debug || $verbose );

  my $scriptFile = create_mothur_script( \@tmp );
  my $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  if ( ! -e $mothurFilteredFile )
  {
    warn "\n\n\tERROR: Count not find $mothurFilteredFile";
    print "\n\n";
    exit 1;
  }

  $cmd = "mv $mothurFilteredFile $algnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

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

# common part of two arrays
sub comm{

  my ($a1, $a2) = @_;

  my @c; # common array
  my %count;

  foreach my $e (@{$a1}, @{$a2}){ $count{$e}++ }

  foreach my $e (keys %count)
  {
    push @c, $e if $count{$e} == 2;
  }

  return @c;
}

# difference of two arrays
sub diff{

  my ($a1, $a2) = @_;

  my (%aa1, %aa2);

  foreach my $e (@{$a1}){ $aa1{$e} = 1; }
  foreach my $e (@{$a2}){ $aa2{$e} = 1; }

  my @d; # dfference array

  foreach my $e (keys %aa1, keys %aa2)
  {
    push @d, $e if exists $aa1{$e} && !exists $aa2{$e};
  }

  return @d;
}

# read table with one column
sub read_array{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in read_array(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
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

sub get_seqIDs_from_fa
{
  my $file = shift;

  my $quiet = 1;
  my $startRun = time();
  my $endRun = time();

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR";
  $/ = ">";
  my $junkFirstOne = <IN>;
  my $count = 1;
  my $timeStr = "";
  my @seqIDs;
  while (<IN>)
  {
    if ( !$quiet && ($count % 500 == 0) )
    {
      $endRun = time();
      my $runTime = $endRun - $startRun;
      if ( $runTime > 60 )
      {
	my $timeMin = int($runTime / 60);
	my $timeSec = sprintf("%02d", $runTime % 60);
	$timeStr = "$timeMin:$timeSec";
      }
      else
      {
	my $runTime = sprintf("%02d", $runTime);
	$timeStr = "0:$runTime";
      }
      print "\r$timeStr";
    }

    chomp;
    my ($id,@seqlines) = split /\n/, $_;
    push @seqIDs, $id;
    $count++;
  }
  close IN;
  $/ = "\n";

  return @seqIDs;
}


# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

exit 0;
