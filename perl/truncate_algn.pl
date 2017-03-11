#!/usr/bin/env perl

=head1 NAME

  truncate_algn.pl

=head1 VERSION

  Versionn 0.92

=head1 CHANGE LOG

  Mar 6, 2017:
  - Fixed the issue with inability of select_seqs.pl not handling symbolic links properly.

  Mar 5, 2017:
  - clstr2_spp.pl use for generating updated taxonomy has been removed as it was
    introducing c<cluster number> type taxonomies. Now the updated taxonomy is
    derived directly from the full length seq's taxonomy


=head1 DESCRIPTION

  truncate_algn.pl will truncate the aligned, full-length 16S rRNA database of
  a taxonomic group of interest to the variable region of choice and place this
  into a new directory. New taxonomy, lineage, and seqID files will be made.
  Additionally, a phylogenetic tree will be produced, re-rooted, and outgroup
  sequences tested for correct placement on the new tree with outgroup_rectifier.pl.

=head1 SYNOPSIS

  truncate_algn.pl -i <group prefix> -v <variable region>

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Group prefix.

=item B<--variable-region, -v>
  16S rRNA gene variable region.
  The current options are V3V4 or V4

=item B<--min-seq-len, -l>
  Sequence minimal length

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

  cd /local/scratch/jbh_pecan_new_tx/nr_files

  truncate_algn.pl -i Chloroflexi -v V3V4 --igs



  cd ~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Proteobacteria_dir

  truncate_algn.pl -i Proteobacteria_group_0 -v V3V4 --min-seq-len 350

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Statistics::Descriptive;
use Cwd qw(abs_path);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "group-prefix|i=s" 	  => \my $grPrefix,
  "variable-region|v=s"   => \my $varReg,
  "min-seq-len|l=i"       => \my $minLen,
  "quiet"           	  => \my $quiet,
  "verbose"           	  => \my $verbose,
  "debug"            	  => \my $debug,
  "dry-run"          	  => \my $dryRun,
  "help|h!"          	  => \my $help,
  "igs"                   => \my $igs,
  "johanna"               => \my $johanna,

  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$grPrefix)
{
  print "\n\nERROR: Missing group prefix\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$varReg)
{
  print "\n\nERROR: Missing variable region.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$minLen)
{
  print "ERROR: Missing sequence min length threshold\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
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
  $usearch6 = "/Users/jholm/bin/usearch6.0.203_i86linux32";
}

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

####################################################################
##                               MAIN
####################################################################


my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  exit;
}

my $trDir = $grPrefix . "_" . $varReg . "_dir";
print "--- Generating truncated data directory $trDir\n";
my $cmd = "mkdir -p $trDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $trPrefix = "$trDir/$grPrefix";

my $origGrPrefix = $grPrefix;
$grPrefix = "$grDir/$grPrefix";

my $txFile       = $grPrefix . "_final.tx";
my $algnFile     = $grPrefix . "_algn_trimmed_final.fa";
my $ogSeqIDsFile = $grPrefix . "_outgroup.seqIDs"; # and here _final_outgroup.seqIDs
my $lineageFile  = $grPrefix . "_final_no_tGTs.lineage";

if ( ! -e $algnFile )
{
  warn "ERROR: $algnFile does not exist";
  exit;
}

if ( ! -e $ogSeqIDsFile )
{
  warn "ERROR: $ogSeqIDsFile does not exist";
  exit;
}

if ( ! -e $txFile )
{
  warn "ERROR: $txFile does not exist";
  exit;
}

if ( ! -e $lineageFile )
{
  warn "ERROR: $lineageFile does not exist";
  exit;
}

# if ( -l $ogSeqIDsFile )
# {
#   $ogSeqIDsFile = readlink($ogSeqIDsFile);
# }

if ( -l $txFile )
{
  $txFile = readlink($txFile);
  if ( ! -s $txFile )
  {
    $txFile = $grDir . "/" . $txFile;
    $txFile = abs_path( $txFile );

    if ( ! -s $txFile )
    {
      warn "ERROR: txaxonomy file, $txFile, does not exist";
      exit;
    }
  }
}

# test consistency of tx and algn_trimmed_final

## also check consistency of OG in og file and sppSeqIDs tree

#print "-- Debug end\n"; exit;

my $trRefFile;
my $trRefFileBasename;
if ($varReg =~ 'V3V4')
{
  $trRefFileBasename = "V400.unique.subsampled.fa";
  #$trRefFile = $dB . $trRefFileBasename;
  print "\n--- Trimming to V3V4 variable region\n"
}
elsif ($varReg =~ 'V4')
{
  $trRefFileBasename = "BEAM.unique.subsampled.fa";
  print "\n--- Trimming to V4 variable region\n"
}
else
{
  print "\nVariable region $varReg not defined or invalid\n\n";
  exit;
}

# check existence
my $dbTrRefFileBasename = $dB . $trRefFileBasename;
if ( ! -e $dbTrRefFileBasename )
{
  warn "ERROR: $dbTrRefFileBasename cannot be found";
  exit;
}

$trRefFile = "$trDir/" . $trRefFileBasename;
$cmd = "cp $dbTrRefFileBasename $trRefFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


print "--- Aligning $trRefFileBasename to $algnFile file\n" if !$quiet;
my @tmp;
push (@tmp,"align.seqs(candidate=$trRefFile, template=$algnFile, processors=4, flip=T)");
printArray(\@tmp, "mothur commands") if ($debug || $verbose);
my $scriptFile = createCommandTxt(\@tmp);

$cmd = "$mothur < $scriptFile; rm -f $scriptFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my @suffixes = (".fasta",".fa",".fna");
my $candBasename = basename($trRefFileBasename, @suffixes); ## This may have to change ($trRefFileBasename to $trRefFile, depending on where mothur writes it)
my $candAlgn = "$trDir/" . $candBasename . ".align";
my $candFile = $candBasename . ".align";

print "--- Calculating alignment range of $candFile\n" if !$quiet;

my $startStats = Statistics::Descriptive::Full->new();
my $endStats = Statistics::Descriptive::Full->new();

my %startTbl;
my %endTbl;

open (IN, "<$candAlgn") or die "Cannot open $candAlgn for reading: $OS_ERROR";
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

my $s = $startStats->mode();
my $e = $endStats->mode();

#print "s: $s\te: $e\n";

print "--- Trimming alignment to $s and $e\n";
my $trAlgnFile = $trPrefix . "_". $varReg . "_algn.fa";
$cmd = "trimAlign -i $algnFile -o $trAlgnFile -s $s -e $e --min-seq-len $minLen";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $trFaFile = $trPrefix . "_". $varReg . ".fa";
print "--- Creating corresonding gap free fasta file $trFaFile\n";
$cmd = "rmGaps -i $trAlgnFile -o $trFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Dereplicatin $trFaFile\n" if !$quiet;
my $trUCfile = $trPrefix . "_". $varReg . ".uc";
my $trNRfile = $trPrefix . "_". $varReg . "_nr.fa";
my $trUCfilelog = $trPrefix . "_". $varReg . "_uc.log";
$cmd = "$usearch6 -cluster_fast $trFaFile -id 1.0 -uc $trUCfile -centroids $trNRfile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

$cmd = "mv $trNRfile $trFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


print "--- Creating non-redundatn seq's taxonomy file\n" if !$quiet;
## extracting seq IDs from the alignment file and selecting those IDs from the taxon file
my $trSeqIDs = $trPrefix  . "_". $varReg . "_nr.seqIDs";
$cmd = "extract_seq_IDs.pl -i $trFaFile -o $trSeqIDs";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $trTxFile = $trPrefix . "_". $varReg . ".tx";
$cmd = "select_tx.pl -s $trSeqIDs -i $txFile -o $trTxFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Creating non-redundatn seq's lineage file\n" if !$quiet;
my $trLineageFile = $trPrefix . "_". $varReg . ".lineage";
$cmd = "select_tx.pl -s $trSeqIDs -i $lineageFile -o $trLineageFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


my $trAlgnFileNR = $trPrefix . "_". $varReg . "_algn_nr.fa";
print "--- Dereplicating truncated alignment file\n" if !$quiet;
$cmd = "rm -f $trAlgnFileNR; select_seqs.pl -s $trSeqIDs -i $trAlgnFile -o $trAlgnFileNR";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Moving $trAlgnFileNR to $trAlgnFile\n";
my $ap = abs_path( $trAlgnFile );
$cmd = "mv $trAlgnFileNR $trAlgnFile; ln -s $ap $trAlgnFileNR";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


my $unrootedTreeFile = $trPrefix . "_". $varReg . "_unrooted.tree";
print "--- Producing phylogenetic tree from $unrootedTreeFile\n";
$cmd = "FastTree -nt $trAlgnFileNR > $unrootedTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Reading outgroup IDs for full length seq's\n" if !$quiet;
my @ogSeqIDs = readArray($ogSeqIDsFile);

printArrayByRow(\@ogSeqIDs, "ogSeqIDs") if $debug;

my @trSeqIDs = readArray($trSeqIDs);
my %trSeqIDsTbl = map { $_ => 1 } @trSeqIDs;

my @trOGs;
for (@ogSeqIDs)
{
  if (exists $trSeqIDsTbl{$_})
  {
    push @trOGs, $_;
  }
}

printArrayByRow(\@trOGs, "trOGs") if $debug;


my $trOGseqIDsFile = $trPrefix . "_". $varReg . "_outgroup.seqIDs";
print "--- Writing update OG seq IDs to $trOGseqIDsFile";
writeArray(\@trOGs, $trOGseqIDsFile);

# print "--- Copying $ogSeqIDsFile to $trOGseqIDsFile";
# $cmd = "cp $ogSeqIDsFile $trOGseqIDsFile";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $rrTreeFile = $trPrefix . "_". $varReg . ".tree";
print "--- Rerooting the tree using outgroup seq's\n" if !$quiet;
$cmd = "rm -f $rrTreeFile; nw_reroot $unrootedTreeFile @trOGs > $rrTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
my $sppSeqIDsFile = $trPrefix . "_". $varReg . "_spp.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $trTxFile > $sppSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $sppSeqIdTreeFile = $trPrefix . "_". $varReg . "_sppSeqIDs.tree";
$cmd = "rm -f $sppSeqIdTreeFile; nw_rename $rrTreeFile $sppSeqIDsFile | nw_order -  > $sppSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Rectifying outgroup sequences - making of them a monophyletic clade if they do not form one\n";
my $truncGr = $origGrPrefix . "_". $varReg;
$cmd = "outgroup_rectifier.pl -i $truncGr";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


print "\n--- Generating a tree with species names at leaves\n";
my $sppTreeFile = $trPrefix . "_" . $varReg . "_spp.tree";
$cmd = "rm -f $sppTreeFile; nw_rename $rrTreeFile $trTxFile | nw_order - > $sppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Generating a condensed tree with species clades collapsed to a single node \n";
my $condSppTreeFile = $trPrefix . "_" . $varReg . "_sppCondensed.tree";
$cmd = "rm -f $condSppTreeFile; nw_condense $sppTreeFile > $condSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;



####################################################################
##                               SUBS
####################################################################

sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file && ! -l $file )
  {
     print "\n\nERROR in readArray() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit;
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
    exit;
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

sub printArrayByRow
{
  my ($a, $header) = @_;

  local $" = '\n ';
  ##local $, = ',';
  print "$header:\n" if $header;
  map { print "$_\n" } @{$a};
  print "\n\n";
}

# write hash table to a file
sub writeTbl{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} keys %tbl;
  close OUT;
  print "Table written to $outFile\n";
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

exit;
