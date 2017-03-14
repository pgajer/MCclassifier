#!/usr/bin/env perl

=head1 NAME

  get_curation_files.pl

=head1 DESCRIPTION

  Given a fasta file of a group of sequences and a file of the corresponding
  outgroups sequence IDs, the script generates

  - fasta file of input and outgroup sequences
  - corresponding alignment file
  - lineage file
  - taxon file

  The naming of the above files is so that the correctness of outgroup selection
  can be tested with outgroup_rectifier.pl script.

=head1 SYNOPSIS

  get_curation_files.pl -i <fasta file> -g <outgroup seqIDs file> -o <out dir> [Options]

=head1 OPTIONS

=over

=item B<--in-tx-file, -i>
  Input fasta file.

=item B<--outgroup-seqID-file, -g>
  List of selected outgroup sequence IDs.

=item B<--output-dir, -o>
  Output directory. Optional.

=item B<--lineage-file, -l>
 Input lineage file. Deprecated.

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

  cd ~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir

  get_curation_files.pl -g Verrucomicrobia_outgroup.seqIDs -i Verrucomicrobia_nr.fa -o Verrucomicrobia_nr_wOG_dir

  or

  get_curation_files.pl -g rdp_Bacteria_phylum_OGs_dir/p_Verrucomicrobia_outgroup.seqID -i Verrucomicrobia_nr.fa
  
  or 
  
  for f in *_nr.fa; do get_curation_files.pl -i $f -g ${f/_nr.fa}_outgroup.seqIDs -l ${f/.fa}.lineage --johanna; done

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
  "outgroup-file|g=s"   => \my $ogSeqIDsFile,
  "output-dir|o=s"      => \my $outDir,
  "seq-file|i=s"        => \my $seqFile,
  "lineage-file|l=s"    => \my $lineageFile,
  "source-lineage|s=s"	=> \my $sourceLineage,
  "igs"                 => \my $igs,
  "verbose|v"           => \my $verbose,
  "quiet"               => \my $quiet,
  "debug"               => \my $debug,
  "debug2"              => \my $debug2,## file name debug
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  "johanna"				=> \my $johanna,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if ( !$ogSeqIDsFile )
{
  print "\n\nERROR: Missing outgroup sequence ID file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$seqFile)
{
  print "\n\nERROR: Missing input sequence file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my $masterLineageFile = "rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis.lineage";
my $masterFaFile      = "rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_nr.fa";

my $mothur = "/Users/pgajer/bin/mothur";

if ( defined $igs )
{
  $mothur = "/usr/local/packages/mothur-1.36.1/mothur";
  $masterFaFile      = "/home/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_nr.fa";
  $masterLineageFile = "/home/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis.lineage";
}

if ( defined $johanna )
{
  $mothur = "/Users/jholm/bin/mothur";
  $masterFaFile      = "/Users/jholm/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_nr.fa";
  $masterLineageFile = "/Users/jholm/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis.lineage";
}

## export LD_LIBRARY_PATH=/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

if ( ! -f $masterFaFile )
{
  print "\n\nERROR: $masterFaFile and does not exist\n\n\n";
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

my @suffixes = ("_nr.fasta","_nr.fa","_nr.fna");
my $faPrefix = basename($seqFile, @suffixes);

my $grPrefix = $faPrefix;

if (!$outDir)
{
  $outDir = $grPrefix . "_dir";
  print ""
}

my $cmd = "mkdir -p $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
##chdir $outDir;

$grPrefix = "$outDir/$grPrefix";

my $newFaFile         = $grPrefix . ".fa";
my $newAlgnFile       = $grPrefix . "_algn.fa";
my $newLineageFile    = $grPrefix . ".lineage";
my $newTxFile         = $grPrefix . ".tx";
my $newOutgroupFile   = $grPrefix . "_outgroup.seqIDs";
my $ogFaFile          = $grPrefix."_outgroup.fa";

my $urTreeFile        = $grPrefix."_unrooted.tree";
my $treeFile          = $grPrefix.".tree";
#my $sppTreeFile       = $grPrefix."_spp.tree";
#my $sppCondTreeFile   = $grPrefix."_sppCondensed.tree";
my $sppSeqIDsTreeFile = $grPrefix."_sppSeqIDs.tree";

if ($igs)
{
  print "\n--- Extracting outgroup sequences from master fasta file\n";
  my $ogFaFile1 = $grPrefix."_tmp_og.fa";
  my $cmd = "select_seqs.pl -i $masterFaFile -s $ogSeqIDsFile -o $ogFaFile1";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Fixing headers of fasta file from outgroup sequences\n";
  $cmd = q(awk '{print$1}')." $ogFaFile1 > $ogFaFile; rm -f $ogFaFile1";
  print "$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
else
{
  print "\n--- Extracting outgroup sequences from master fasta file\n";
  my $cmd = "select_seqs.pl -i $masterFaFile -s $ogSeqIDsFile -o $ogFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## Combine this with the input fasta file
print "--- Concatenating outgroup & input sequences\n";
$cmd = "cat $ogFaFile $seqFile > $newFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## Align, make output the
print "--- Aligning sequences\n";
$cmd = "mafft --auto --inputorder $newFaFile > $newAlgnFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Parsing master lineage table\n";
my %mLineageTbl = read2colTbl($masterLineageFile);

print "--- Extracting sequence IDs from fasta file\n";
my @seqIDs = get_seqIDs_from_fa($newFaFile);

print "--- Creating lineage file\n";
my %lineageTbl;
@lineageTbl{@seqIDs} = @mLineageTbl{@seqIDs};

writeTbl(\%lineageTbl, $newLineageFile);

print "--- Reading outgroup sequence IDs\n";
my @ogSeqIDs = readArray($ogSeqIDsFile);
my %ogSeqIDsTbl = map{$_ =>1} @ogSeqIDs;

print "--- Creating taxonomy file including OG's\n";
open OUT, ">$newTxFile" or die "Cannot open $newTxFile for writing: $OS_ERROR\n";
my $ogCount = 1;
for ( @seqIDs )
{
  my $lineage = $lineageTbl{$_};
  my @f = split ";", $lineage;
  my $sp = pop @f;

  if (exists $ogSeqIDsTbl{$_})
  {
    print OUT "$_\t$sp" . "_OG\n";
    ##print OUT "$_\tOUTGROUP_$ogCount\n";
    $ogCount++;
  }
  else
  {
    print OUT "$_\t$sp\n";
  }
}
close OUT;

#if ( !$treeFile )
if (1)
{
  print "--- Generating sequence summary of the alignment file\n";
  my @tmp;
  push (@tmp,"summary.seqs(fasta=$newAlgnFile)");
  printArray(\@tmp, "mothur commands") if ($debug || $verbose);

  my $scriptFile = createCommandTxt(\@tmp);
  $cmd = "rm -f mothur.log; $mothur < $scriptFile > mothur.log 2>&1; rm -f $scriptFile";
  #$cmd = "$mothur < $scriptFile; rm -f $scriptFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Trimming alignment to length of 95% of sequences\n";
  my $summaryFile     = $grPrefix . "_algn.summary";
  my $trimmedAlgnFile = $grPrefix . "_algn_trimmed.fa";
  $cmd = "trim_align.pl -c 95 -j $summaryFile -i $newAlgnFile -o $trimmedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Generating a phylogenetic tree based on trimmed alignment\n";
  $cmd = "FastTree -nt $trimmedAlgnFile > $urTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;# && !$skipFastTree;

  print "--- Rerooting the tree using outgroup sequences\n";
  $cmd = "rm -f $treeFile; nw_reroot $urTreeFile @ogSeqIDs | nw_order -  > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # print "--- Generating a tree with species names at leaves\n";
  # my $cmd = "rm -f $sppTreeFile; nw_rename $treeFile $newTxFile | nw_order - > $sppTreeFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
  my $sppSeqIDsFile     = $grPrefix.".sppSeqIDs";
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $newTxFile > $sppSeqIDsFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $sppSeqIDsTreeFile; nw_rename $treeFile $sppSeqIDsFile | nw_order -  > $sppSeqIDsTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # print "--- Generating a condensed tree with species clades collapsed to a single node \n";
  # $cmd = "rm -f $sppCondTreeFile; nw_condense $sppTreeFile > $sppCondTreeFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;# && !$skipCondensed;
}


print "\nCreated the following files\n";
print "\t$newFaFile\n";
print "\t$newAlgnFile\n";
print "\t$newLineageFile\n";
print "\t$newTxFile\n";

if ( $ogSeqIDsFile ne $newOutgroupFile )
{
  print "Producing symbolic link for $newOutgroupFile to $ogSeqIDsFile";
  my $cmd = "rm -f $newOutgroupFile; ln -s $ogSeqIDsFile $newOutgroupFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  print "\t$newOutgroupFile\n";
}

print "\t$urTreeFile\n";
print "\t$treeFile\n";
#print "\t$sppTreeFile\n";
#print "\t$sppCondTreeFile\n";
print "\t$sppSeqIDsTreeFile\n";

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
print "\n";


####################################################################
##                               SUBS
####################################################################

sub get_seqIDs_from_fa
{
  my $file = shift;

  my $quiet = 1;
  my $startRun = time();
  my $endRun = time();

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <IN>;
  my $count = 1;
  my $timeStr = "";
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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
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

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# read table with one column
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    print "\n\nERROR in readArray() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
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

exit;
