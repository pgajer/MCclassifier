#!/usr/bin/env perl

=head1 NAME

  add_to_phylogroups.pl


=head1 DESCRIPTION
  Summary:
  add_to_phylogroups.pl will take as input an unaligned fasta file 
  of sequences truncated to the indicated variable regions and, using
  PECAN classification with the current database, will combine them
  with those of the correct phylogroup and perform dereplication. 
  Additionally, other files required for taxonomy_cleanup.pl for each 
  phylogroup will be updated. Currently, the updated files will be 
  produced in the same directory as the input fasta file. 

  More details:
  The script will run the input sequence file through PECAN classifier
  (currently only V3V4 is working). The taxonomy in the MC_output file 
  from this run will be used to extract the phylogroup from the sp_to_phGr.txt 
  file. Using the phylogroup information, the lineages from the current, 
  phylogroup PECAN lineage files. The species level annotation will be replaced
  with the annotation given through the source database of the input
  sequence file. This lineage file will be concatenated with that of the
  corresponding phylogroup. The first and last columns of this file
  will be used to produce the taxonomy file for the input sequences. This
  file will be concatenated with the phylogroup taxonomy file.  

  The mothur command align.seqs will be used to align the input sequences
  with those of the corresponding phylogroup and concatenation of these 
  and files will produce the concatenated alignment file. The sequences will 
  then be dereplicated. The resulting sequence IDs of the _nr.fa file will used
  to extract the file lineage and taxonomy files. A tree will be produced
  from the final _nr.fa alignment.  

=head1 SYNOPSIS

  add_to_phylogroups.pl -i <sequence file> -v <variable region>

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Input sequence file. Must be unaligned.

=item B<--variable-region, -v>
  16S rRNA gene variable region.
  The current options are V3V4 or V4

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

  cd ~/local/scratch/archaea_pecan/gg_silva_rdpArchaea/V4

  trim_seqs.pl -i SILVA_128_SSURef_Nr99_tax_silva_16S.fasta -v V3V4


=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Statistics::Descriptive;
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-sequences|i=s"   => \my $seqFile,
  "variable-region|v=s"   => \my $varReg,
  "verbose"               => \my $verbose,
  "debug"                 => \my $debug,
  "dry-run"               => \my $dryRun,
  "help|h!"               => \my $help,
  "igs"                   => \my $igs,
  "johanna"               => \my $johanna,

  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$varReg)
{
  print "\n\n ERROR: Must provide a variable region with (-v)\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

elsif (!$seqFile)
{
  print "\n\nERROR: Missing input sequence file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $mothur   = "/Users/pgajer/bin/mothur";
my $dB       = "/Users/pgajer/devel/MCextras/data/RDP/";
my $usearch6 = "/Users/pgajer/bin/usearch6.0.203_i86osx32";

if ( defined $igs )
{
  $mothur   = "/usr/local/packages/mothur-1.36.1/mothur";
  $dB       = "/usr/local/projects/pgajer/devel/MCextras/data/RDP/";
  $usearch6 = "/local/projects/pgajer/bin/usearch6.0.203_i86linux32";
}

if ( defined $johanna )
{
  $mothur = "/Users/jholm/bin/mothur";
  $dB     = "/Users/jholm/RDP/";
  $usearch6 = "/Users/jholm/bin/usearch9.2.64_i86osx32";
}

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

####################################################################
##                               MAIN
####################################################################

my $debugStr = "";
$debugStr = "--debug" if $debug;

if ( ! -e $seqFile )
{
  warn "ERROR: $seqFile does not exist";
  exit 1;
}

my $pecanFile;
my $sp_to_phylogroup;

## Set the PECAN db to the correct location based on the variable region.
if ($varReg =~ 'V3V4')
{
  $pecanFile = "/local/projects/pgajer/devel/MCextras/data/all_bacteria_V3V4_MC_models_dir/";
  print "--- The V3V4 PECAN database is located at $pecanFile.\n"
  $sp_to_phylogroup = "/local/projects/pgajer/devel/MCextras/data/RDP/V3V4/sp_to_phGr.txt";
}
elsif ($varReg =~ 'V4')
{
  $pecanFile = "/local/projects/pgajer/devel/MCextras/data/all_bacteria_V4_MC_models_dir/";
  print "--- The V4 PECAN database does not yet exist.\n"
  $sp_to_phylogroup = "/local/projects/pgajer/devel/MCextras/data/RDP/V4/sp_to_phGr.txt";
  exit 1;
}
else
{
  warn "\n\n\tERROR: Invalid variable region $varReg";
  print "\n\n";
  exit 1;
}

my @suffixes = (".fasta",".fa",".fna");
my $inputBasename = basename($seqFile, @suffixes);
my %sp_to_phGr = readTbl($sp_to_phylogroup)

## Classify the input sequences to the current PECAN database. 
print "--- Classifying $seqFile file the $varReg PECAN database.\n" if !$quiet;
my $classifyDir = $inputBasename . "_" . $varReg . "classify_dir";
$cmd = "classify -v --skip-err-thld -d $pecanFile -i $seqFile -o $classifyDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

## Using the taxonomy resulting from PECAN, determine the phylogroup for each sequence. 
my $classifyTxFile = $classifyDir . "/MC_order7_results.txt";
my @classifyTx = read2colTbl($classifyTxFile)

my %seqID_to_phyGr;
foreach my $i (@classify)
  {
    my ($seqID, $tx) = split /\t/, $_;
    chomp $tx;

    if ($tx =~ keys %sp_to_phGr)
    {
      ## make a hash table with the keys as phylogroups and the seqIDs as values. phGrp => {S1, S2, ... , S3}
      $seqID_to_phyGr{$sp_to_phGr{$tx}} = $seqID;
      push @phylogroups, $sp_to_phGr{$tx};
    }
    else 
    {
      print "--- $tx not found in $sp_to_phylogroup \n";
      print "--- Make sure $sp_to_phylogroup contains all data \n";
      print "--- for current $varReg PECAN database. \n\n";
      exit;
    }
  }

## Using the seqID to phylogroup file, we have to begin to divide the input sequences
## by the phylogroup, and join the .fa, lineage, and .tx files, as well as make another
## tree. This could all be put into a loop, which iterates over the phylogroups to be
## updated. 

foreach my $phGr (keys %seqID_to_phyGr)
  {
    ## use the phyGr to pull out all seqIDs matching it. then print this array as the
    ## phylogroup.seqID, and extract the seq's of this group. Then use the phylogroup
    ## to obtain the oldPECAN phylogroup .fa file to perform align.seqs.

  }

my @tmp;
push (@tmp,"align.seqs(candidate=$trRefFile, template=$seqFile, processors=4, flip=T)");
printArray(\@tmp, "mothur commands") if ($debug || $verbose);
my $scriptFile = createCommandTxt(\@tmp);

$cmd = "$mothur < $scriptFile; rm -f $scriptFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


my @suffixes = (".fasta",".fa",".fna");
my $candBasename = basename($trRefFile, @suffixes); ## This may have to change ($trRefFileBasename to $trRefFile, depending on where mothur writes it)
my $candAlgn = $candBasename . ".align";
my $candFile = "/local/projects/pgajer/devel/MCextras/data/RDP/" . $candAlgn;


## removing $trRefFile as it is not needed anymore
#$cmd = "rm -f $trRefFile";
#print "\tcmd=$cmd\n" if $dryRun || $debug;
#system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Calculating alignment range for $varReg region of $candAlgn\n" if !$quiet;

my $startStats = Statistics::Descriptive::Full->new();
my $endStats = Statistics::Descriptive::Full->new();

my %startTbl;
my %endTbl;

open (IN, "<$candFile") or die "Cannot open $candFile for reading: $OS_ERROR";
$/ = ">";
my $junkFirstOne = <IN>;
while (<IN>)
{
  chomp;
  my ($def,@seqlines) = split /\n/, $_;
  my $seq = join '', @seqlines;
  my ($id) = split /\s+/, $def;
  $id =~ s/^>//;

  my @seq = split "", $seq;

  # find first non-gapped position
  my $startPos = 0;
  while ( $seq[$startPos] !~ /\w/ )
  {
    $startPos++;
  }

  # find last non-gapped position
  my $endPos = $#seq;
  while ( $seq[$endPos] !~ /\w/ )
  {
    $endPos--;
  }

  $startStats->add_data($startPos);
  $endStats->add_data($endPos);
}
$/ = "\n";
close IN;

print "\r                                              ";
print "\nNumber of sequences:\t" . $startStats->count() . "\n\n";

print "Start summary stats [0-based positions]\n";
print "Min:\t" . $startStats->min() . "\n";
print "Max:\t" . $startStats->max() . "\n";
print "Mode:\t" . $startStats->mode() . "\n";
print "Median:\t" . $startStats->median() . "\n";
print "IQR:\t" . $startStats->percentile(25) . "-" . $startStats->percentile(75) . "\n";
print "Mean:\t" . $startStats->mean() . "\n";

print "\nEnd summary stats [0-based positions]\n";
print "Min:\t" . $endStats->min() . "\n";
print "Max:\t" . $endStats->max() . "\n";
print "Mode:\t" . $endStats->mode() . "\n";
print "Median:\t" . $endStats->median() . "\n";
print "IQR:\t" . $endStats->percentile(25) . "-" . $endStats->percentile(75) . "\n";
print "Mean:\t" . $endStats->mean() . "\n\n";


if (defined $manual)
{
  print "--- Please examine above table and indicate the base to START truncation: \n";
  chomp ($s = <STDIN>);
  print "--- Please examine above table and indicate the base to END truncation: \n";
  chomp ($e = <STDIN>);
}
else
{
  print "--- Automatically using mode positions for trimming. \n";
  $s = $startStats->mode();
  $e = $endStats->mode();
  print "--- Start position = $s. \n";
  print "--- End position = $e. \n";
}

if ($start && $end)
{
  print "--- Trimming start position provided. Trimming at position $start. \n";
  chomp ($s = $start);
  print "--- Trimming start position provided. Trimming at position $end. \n";
  chomp ($e = $end);
}

#print "s: $s\te: $e\n";

my @suffixes = (".fasta",".fa",".fna");
my $trPrefix = basename($seqFile, @suffixes); 

print "--- Trimming alignment to $s and $e\n";
my $trAlgnFile = $trPrefix . "_" . $varReg . "_algn.fa";
$cmd = "trimAlign -i $seqFile -o $trAlgnFile -s $s -e $e --min-seq-len $minLen";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $trFaFile = $trPrefix . "_" . $varReg . ".fa";
print "--- Creating corresonding gap free fasta file $trFaFile\n";
$cmd = "rmGaps -i $trAlgnFile -o $trFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Dereplicating $trFaFile\n" if !$quiet;
my $trUCfile = $trPrefix . "_" . $varReg . ".uc";
my $trNRfile = $trPrefix . "_" . $varReg . "_nr.fa";
my $trUCfilelog = $trPrefix . "_" . $varReg . "_uc.log";
$cmd = "$usearch6 -cluster_fast $trFaFile -id 1.0 -uc $trUCfile -centroids $trNRfile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $nrSeqIDs = $trPrefix . "_" . $varReg . "_nr.seqIDs";
$cmd = "extract_seq_IDs.pl -i $trNRfile -o $nrSeqIDs";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $trAlgnFileNR = $trPrefix . "_" . $varReg . "_algn_nr.fa";
print "--- Dereplicating truncated alignment file\n" if !$quiet;
$cmd = "rm -f $trAlgnFileNR; select_seqs.pl -s $nrSeqIDs -i $trAlgnFile -o $trAlgnFileNR";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


####################################################################
##                               SUBS
####################################################################

sub createCommandTxt{

    my (@arr) = @{$_[0]};
    my $file = "mothur_script.txt"; ##tmpnam();
    open OUT, ">$file" or die "Cannot open file $file to write: $!\n";
    foreach my $c (@arr){
        print OUT $c . "\n";
    }
    print OUT "quit()\n";

    return $file;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR: $file does not exist\n\n\n";
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

sub read2colTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
    exit 1;
  }

  my @tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my @tbl = split /\s+/,$_;
  }
  close IN;

  return %tbl;
}

exit;

