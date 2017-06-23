#!/usr/bin/env perl

=head1 NAME

  truncate_algn.pl


=head1 DESCRIPTION

  truncate_algn.pl will truncate the aligned, full-length 16S rRNA database of
  a taxonomic group of interest to the variable region of choice and place this
  into a new directory. New lineage, alignment, outgroup files will be made.
  Additionally, a phylogenetic tree will be produced, re-rooted, and outgroup
  sequences tested for correct placement on the new tree with outgroup_rectifier.pl.

  The multiple sequence alignment produced by taxonomy_cleanup.pl (which is run
  before this script) may leave it out of synch with the lineage file. Therefore,
  the first order of business is to make sure we use the same sequence IDs for
  both and that we still have at least one outgroup sequence there.


=head1 SYNOPSIS

  truncate_algn.pl -i <group prefix> -v <variable region>

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Group prefix.

=item B<--variable-region, -v>
  16S rRNA gene variable region.
  The current options are V3V4 or V4

=item B<--min-seq-len, -l>
  Sequence minimal length

=item B<--min-size, -m>
  Minimal number of elements in fasta/lineage files for this script to run.

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

  cd ~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Proteobacteria_dir

  truncate_algn.pl -i Firmicutes_group_0 -v V3V4 --min-seq-len 350

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

my $minSize       = 10;

GetOptions(
  "group-prefix|i=s" 	  => \my $grPrefix,
  "variable-region|v=s"   => \my $varReg,
  "min-seq-len|l=i"       => \my $minLen,
  "min-size|m=i"          => \$minSize,
  "quiet"           	  => \my $quiet,
  "verbose"           	  => \my $verbose,
  "debug"            	  => \my $debug,
  "dry-run"          	  => \my $dryRun,
  "help|h!"          	  => \my $help,
  "igs"                   => \my $igs,
  "johanna"               => \my $johanna,
  "manual"		  => \my $manual,

  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$grPrefix)
{
  print "\n\nERROR: Missing group prefix\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$varReg)
{
  print "\n\nERROR: Missing variable region.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$minLen)
{
  print "ERROR: Missing sequence min length threshold\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $debugStr = "";
$debugStr = "--debug" if $debug;

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  exit 1;
}

my $mothur   = "/Users/pgajer/bin/mothur";
my $dB       = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
my $usearch6 = "/Users/pgajer/bin/usearch6.0.203_i86osx32";

if ( defined $igs )
{
  $mothur   = "/usr/local/packages/mothur-1.36.1/mothur";
  $dB       = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
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
##                           I/O CONFIGURATION
####################################################################

# Generating truncated data directory
my $trDir = $grPrefix . "_" . $varReg . "_dir";
my $cmd = "rm -rf $trDir; mkdir -p $trDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $trPrefix = "$trDir/$grPrefix" . "_". $varReg;

# Generating truncated temp directory
# This is where all but the final files will be stored
my $trTmpDir = $trDir . "/trTempDir";
$cmd = " mkdir -p $trTmpDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $trTmpPrefix = "$trTmpDir/$grPrefix" . "_". $varReg;

my $origGrPrefix = $grPrefix;
$grPrefix = "$grDir/$grPrefix";

my $trRefFile;
my $trRefFileBasename;
if ($varReg =~ 'V3V4')
{
  $trRefFileBasename = "V400.unique.subsampled.fa";
  #$trRefFile = $dB . $trRefFileBasename;
  print "--- Detected ref db for the V3V4 variable region\n"
}
elsif ($varReg =~ 'V4')
{
  $trRefFileBasename = "BEAM.unique.subsampled.fa";
  print "--- Detected ref db for the V4 variable region\n"
}
else
{
  warn "\n\n\tERROR: Invalid variable region $varReg";
  print "\n\n";
  exit 1;
}

# check existence
my $dbTrRefFileBasename = $dB . $trRefFileBasename;
if ( ! -e $dbTrRefFileBasename )
{
  warn "ERROR: $dbTrRefFileBasename cannot be found";
  exit 1;
}

$trRefFile = "$trDir/" . $trRefFileBasename;
$cmd = "cp $dbTrRefFileBasename $trRefFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

##
## main input files
##
my $liFile   = $grPrefix . "_spp.lineage";
my $algnFile = $grPrefix . "_algn_trimmed_final.fa";
my $ogFile   = $grPrefix . "_outgroup.seqIDs"; # and here _final_outgroup.seqIDs

my $errorFile = $trTmpPrefix . "_outgroup_rectifier_ERROR"; # this is a file that indicates that outgroup_rectifier.pl finished with an error
if ( -e $errorFile )
{
  unlink($errorFile);
}

if ( ! -e $algnFile )
{
  warn "ERROR: $algnFile does not exist";
  exit 1;
}

if ( -l $algnFile )
{
  $algnFile = readlink($algnFile);
}

if ( ! -e $ogFile )
{
  warn "ERROR: $ogFile does not exist";
  exit 1;
}

if ( ! -e $liFile )
{
  warn "ERROR: $liFile does not exist";
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

# The multiple sequence alignment produced by taxonomy_cleanup.pl (which is run
# before this script) may leave it out of synch with the lineage file. Therefore,
# the first order of business is to make sure we use the same sequence IDs for
# both and that we still have at least one outgroup sequence there.

print "--- Extracting seq IDs from trimmed alignment fasta file\n";
my @seqIDs = get_seqIDs_from_fa( $algnFile );

if ( @seqIDs < $minSize )
{
  print "\n\nWARNING: $grPrefix has less than $minSize elements; Exiting\n\n";
  exit 0;
}

print "--- Parsing lineage table\n";
my %liTbl = read_tbl( $liFile );

## Testing if lineage and fa files have the same seq IDs
print "--- Checking if seqIDs of $algnFile and $liFile are the same\n";
my @lSeqIDs = keys %liTbl;
my @commSeqIDs = comm(\@seqIDs, \@lSeqIDs);
if (@commSeqIDs != @seqIDs || @commSeqIDs != @lSeqIDs)
{
  warn "\n\tWARNING: seq IDs of trimmed alignment fasta file and lineage file do not match";
  print "\tNumber of elements in the trimmed alignment file: " . @seqIDs . "\n";
  print "\tNumber of elements in the lineage file: " . @lSeqIDs . "\n";
  print "\tNumber of elements common to the alignment and lineage files: " . @commSeqIDs . "\n";
  print "\tNOTE: the sequences common to both will be used to construct truncated alignment\n\n";

  if ( @lSeqIDs > @commSeqIDs )
  {
    my @d = diff(\@lSeqIDs, \@commSeqIDs);
    delete @liTbl{@d};
  }

  if (@commSeqIDs != @seqIDs)
  {
    my $commSeqIDsFile = $trTmpPrefix . "_comm.seqIDs";
    writeArray(\@commSeqIDs, $commSeqIDsFile);

    print "---  Selecting seq IDs common to the lineage and the alignment\n";
    my $tmpAlgnFile = $trTmpPrefix . "_pruned_algn.fa";
    my $cmd = "select_seqs.pl -s $commSeqIDsFile -i $algnFile -o $tmpAlgnFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

    $algnFile = $tmpAlgnFile;
  }
}

print "--- Parsing $ogFile\n";
my @ogSeqIDs = read_array( $ogFile );

print "--- Testing if outgroup sequences are part of seqIDs\n";
my @ogDiff = diff( \@ogSeqIDs, \@commSeqIDs );
@ogSeqIDs = comm( \@ogSeqIDs, \@commSeqIDs );

if ( scalar(@ogDiff) != 0 )
{
  warn "\n\tWARNING: the following outgroup seq IDs are not in the trimmed alignment file:\n\n";
  print_array(\@ogDiff);
}

if ( scalar(@ogSeqIDs) == 0 )
{
  warn "\n\tERROR: All outgroup seq's were lost";
  print "\n\n";
  exit 1;
}

print "\n\tNumber of seq's present in the trimmed alignment and lineage files: " . @commSeqIDs . "\n";
print "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n\n";



print "--- Aligning $trRefFileBasename to $algnFile file\n" if !$quiet;
my @tmp;
push (@tmp,"align.seqs(candidate=$trRefFile, template=$algnFile, processors=4, flip=T)");
print_array(\@tmp, "mothur commands") if ($debug || $verbose);
my $scriptFile = createCommandTxt(\@tmp);

$cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my @suffixes = (".fasta",".fa",".fna");
my $candBasename = basename($trRefFileBasename, @suffixes); ## This may have to change ($trRefFileBasename to $trRefFile, depending on where mothur writes it)
my $candAlgn = "$trDir/" . $candBasename . ".align";
my $candFile = $candBasename . ".align";


print "--- Removing $trRefFile as it is not needed anymore\n";
$cmd = "rm -f $trRefFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


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

my $s;
my $e;
if (defined $manual)
{
  print "Please examine above table and indicate the base to START truncation: \n";
  chomp ($s = <STDIN>);
  print "Please examine above table and indicate the base to END truncation: \n";
  chomp ($e = <STDIN>);
}
else
{
  $s = $startStats->mode();
  $e = $endStats->mode();
}

#print "s: $s\te: $e\n";

print "--- Trimming alignment to $s and $e\n";
my $trAlgnFile = $trPrefix . "_algn.fa";
$cmd = "trimAlign -i $algnFile -o $trAlgnFile -s $s -e $e --min-seq-len $minLen";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $trFaFile = $trPrefix . ".fa";
print "--- Creating corresonding gap free fasta file $trFaFile\n";
$cmd = "rmGaps -i $trAlgnFile -o $trFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Dereplicating $trFaFile\n" if !$quiet;
my $trUCfile = $trTmpPrefix . ".uc";
my $trNRfile = $trTmpPrefix . "_nr.fa";
my $trUCfilelog = $trTmpPrefix . "_uc.log";
$cmd = "$usearch6 -cluster_fast $trFaFile -id 1.0 -uc $trUCfile -centroids $trNRfile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

$cmd = "mv $trNRfile $trFaFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


print "--- Creating non-redundant seq's taxonomy file\n" if !$quiet;
## extracting seq IDs from the alignment file and selecting those IDs from the taxon file
my $nrSeqIDs = $trTmpPrefix  . "_nr.seqIDs";
$cmd = "extract_seq_IDs.pl -i $trFaFile -o $nrSeqIDs";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Creating non-redundant seq's lineage file\n" if !$quiet;
my $trLineageFile = $trPrefix . ".lineage";
$cmd = "select_tx.pl -s $nrSeqIDs -i $liFile -o $trLineageFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Dereplicating truncated alignment file\n" if !$quiet;
my $trAlgnFileNR = $trTmpPrefix . "_algn_nr.fa";
$cmd = "rm -f $trAlgnFileNR; select_seqs.pl -s $nrSeqIDs -i $trAlgnFile -o $trAlgnFileNR";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


print "--- Moving $trAlgnFileNR to $trAlgnFile\n";
my $ap = abs_path( $trAlgnFile );
$cmd = "mv $trAlgnFileNR $trAlgnFile; ln -s $ap $trAlgnFileNR";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


my $unrootedTreeFile = $trTmpPrefix . "_unrooted.tree";
print "--- Producing phylogenetic tree from $unrootedTreeFile\n";
$cmd = "FastTree -nt $trAlgnFileNR > $unrootedTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Identifying outgroup seq's in the non-redundant set of truncated seq's\n" if !$quiet;
print_array( \@ogSeqIDs, "\nOG seqIDs BEFORE selection of non-redundant OG seq's") if $debug;

my @trSeqIDs = read_array($nrSeqIDs);
my @trOGs = comm(\@trSeqIDs, \@ogSeqIDs);

if (@trOGs == 0)
{
  warn "\n\n\tERROR: All outgroup sequences has been lost";
  print "\n\n\n";
  exit 1;
}

print_array(\@trOGs, "OG seqIDs AFTER selection of non-redundant OG seq's") if $debug;


my $trOGseqIDsFile = $trPrefix . "_outgroup.seqIDs";
print "--- Writing update OG seq IDs to $trOGseqIDsFile";
writeArray(\@trOGs, $trOGseqIDsFile);

my %ogInd = map{$_ =>1} @trOGs; # outgroup elements indicator table

my $rrTreeFile = $trPrefix . ".tree";
print "--- Rerooting the tree using outgroup seq's\n" if !$quiet;
$cmd = "rm -f $rrTreeFile; nw_reroot $unrootedTreeFile @trOGs > $rrTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Extracting species taxonomy from the lineage data\n";

my %spp;
my $nOG = 0;
for my $id ( keys %liTbl )
{
  my $lineage = $liTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  if ( exists $ogInd{$id} )
  {
    $sp .= "_OG";
    $nOG++;
  }
  $spp{$id} = $sp;
}

if ( $nOG != @trOGs )
{
  warn "";
  print "\n\n";
  exit 1;
}

my $trTxFile = $trPrefix . ".tx";
print "--- Writing species taxonomy to $trTxFile\n";
writeTbl(\%spp, $trTxFile);

print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
my $sppSeqIDsFile = $trTmpPrefix . "_spp.seqIDs";
$cmd = "rm -f $sppSeqIDsFile; awk '{print \$1\"\\t\"\$2\"__\"\$1}' $trTxFile > $sppSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

my $sppSeqIdTreeFile = $trTmpPrefix . "_sppSeqIDs.tree";
$cmd = "rm -f $sppSeqIdTreeFile; nw_rename $rrTreeFile $sppSeqIDsFile | nw_order -  > $sppSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

print "--- Generating a tree with species names at leaves\n";
my $sppTreeFile = $trTmpPrefix . "_spp.tree";
$cmd = "rm -f $sppTreeFile; nw_rename $rrTreeFile $trTxFile | nw_order - > $sppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

if (@trOGs > 1)
{
  print "--- Rectifying outgroup sequences - making a monophyletic clade if they do not form one\n";
  my $truncGr = $origGrPrefix  . "_" . $varReg;
  $cmd = "outgroup_rectifier.pl $debugStr -i $truncGr";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  ## In principle we should test if outgroup_rectifier.pl made any changes and
  ## only then recreate the following trees. Yet, since its so fast we can as
  ## well done it without any checks.
  print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
  my $sppSeqIDsFile = $trTmpPrefix . "_spp.seqIDs";
  $cmd = "rm -f $sppSeqIDsFile; awk '{print \$1\"\\t\"\$2\"__\"\$1}' $trTxFile > $sppSeqIDsFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  my $sppSeqIdTreeFile = $trTmpPrefix . "_sppSeqIDs.tree";
  $cmd = "rm -f $sppSeqIdTreeFile; nw_rename $rrTreeFile $sppSeqIDsFile | nw_order -  > $sppSeqIdTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  print "--- Generating a tree with species names at leaves\n";
  my $sppTreeFile = $trTmpPrefix . "_spp.tree";
  $cmd = "rm -f $sppTreeFile; nw_rename $rrTreeFile $trTxFile | nw_order - > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

}

print "--- Generating a condensed tree with species clades collapsed to a single node \n";
my $condSppTreeFile = $trTmpPrefix . "_spp_condensed.tree";
$cmd = "rm -f $condSppTreeFile; nw_condense $sppTreeFile > $condSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "\n\n";

## Checking if outgroup_rectifier produced an error
$errorFile = $trTmpPrefix . "_outgroup_rectifier_ERROR"; # this is a file that indicates that outgroup_rectifier.pl finished with an error
if ( -e $errorFile )
{
  print "\n\tERROR: outgroup_rectifier.pl finished with an error!\n";
  print "\tPlease fix the source of the error\n\n";
}

## Getting rid of truncation process intermediate files
$cmd =  "rm -f $trDir/$candBasename" . ".*";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

####################################################################
##                               SUBS
####################################################################

sub read_array
{
  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file && ! -l $file )
  {
     print "\n\nERROR in read_array() at line " . __LINE__ . ": $file does not exist\n\n\n";
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
sub print_array
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


# write hash table to a file
sub writeTbl{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} keys %tbl;
  close OUT;
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl{

  my $file = shift;

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

exit 0;
