#!/usr/bin/env perl

=head1 NAME

  trim_seqs.pl


=head1 DESCRIPTION

  trim_seqs.pl will take a full-length, aligned set of 16S rRNA gene 
  sequences, and produce an alignment truncated to either the variable region
  of choice (currently V3V4 or V4) or the indicated start and end positions. 
  Dereplication of sequences will also be performed via usearch6.
  Three files will be produced: a truncated alignment file, a dereplicated
  truncated alignment file, a truncated ungapped sequence file, and a 
  dereplicated truncated ungapped sequence file.

=head1 SYNOPSIS

  trim_seqs.pl -i <sequence file> -l <minimum sequence length> -v <variable region> | -s <start position> -e <end position>

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Input sequence file. Must be alignment!!

=item B<--variable-region, -v>
  16S rRNA gene variable region.
  The current options are V3V4 or V4

=item B<--min-seq-len, -l>
  Sequence minimal length

=item B<--start-position, -s>
Start position to trim alignment

=item B<--start-position, -e>
End position to trim alignment

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

  cd ~/local/scratch/archaea_pecan/gg_silva_rdpArchaea

  trim_seqs.pl -i SILVA_128_SSURef_Nr99_tax_silva_16S.fasta -v V3V4 -l 350

  or

  trim_seqs.pl -i SILVA_128_SSURef_Nr99_tax_silva_16S.fasta -s 1200 -e 1600 -l 350

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
  "input-sequences|i=s" 	=> \my $seqFile,
  "variable-region|v=s"   => \my $varReg,
  "min-seq-len|l=i"       => \my $minLen,
  "quiet"           	    => \my $quiet,
  "verbose"           	  => \my $verbose,
  "debug"            	    => \my $debug,
  "dry-run"          	    => \my $dryRun,
  "help|h!"          	    => \my $help,
  "igs"                   => \my $igs,
  "johanna"               => \my $johanna,
  "manual"		            => \my $manual,
  "start|s=i"             => \my $start,
  "end|e=i"               => \my $end,

  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$varReg && !$start && !$end)
{
  print "\n\n ERROR: Must provide either a variable region with (-v)\n";
  print "\n\n or start (-s) and end (-e) trimming positions.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$varReg && !$start && $end)
{
  print "\n\n ERROR: End position provided, but start\n";
  print "\n\n position (-s) missing.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$varReg && $start && !$end)
{
  print "\n\n ERROR: Start position provided, but end\n";
  print "\n\n position (-e) missing.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$seqFile)
{
  print "\n\nERROR: Missing input sequence file.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$minLen)
{
  print "ERROR: Missing sequence min length threshold\n\n";
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


my $s;
my $e;
my $cmd;

if (defined $varReg)
{
    
    my $trRefFile;
    my $trRefFileBasename;
    if ($varReg =~ 'V3V4')
    {
      $trRefFile = "/local/projects/pgajer/devel/MCextras/data/RDP/V400.unique.subsampled.fa";
      print "--- Detected ref db for the V3V4 variable region\n"
    }
    elsif ($varReg =~ 'V4')
    {
      $trRefFile = "/local/projects/pgajer/devel/MCextras/data/RDP/BEAM.unique.subsampled.fa";
      print "--- Detected ref db for the V4 variable region\n"
    }
    else
    {
      warn "\n\n\tERROR: Invalid variable region $varReg";
      print "\n\n";
      exit 1;
    }

    # check existence
    if ( ! -e $trRefFile )
    {
      warn "ERROR: $trRefFile cannot be found";
      exit 1;
    }
    ## The unaligned sequence database must be aligned to 

    print "--- Aligning $trRefFile to $seqFile file\n" if !$quiet;
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
my $cmd = "trimAlign -i $seqFile -o $trAlgnFile -s $s -e $e --min-seq-len $minLen";
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
$cmd = "extract_seq_IDs.pl -i $trNRFile -o $nrSeqIDs";
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

exit;

