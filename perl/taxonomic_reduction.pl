#!/usr/bin/env perl

=head1 NAME

  taxonomic_reduction.pl


=head1 DESCRIPTION

  Given a fasta file of sequences from some taxonomic rank, subsample species so
  that none of species has more than maxSp sequences.

  Groups families or higher taxonomic ranks into groups of size maxGr.

  First, check if after subsampling all families are within maxGr size. If not,
  one needs to modify maxSp or maxGr.

=head1 SYNOPSIS

  taxonomic_reduction.pl --max-group-size <maxGr> --max-sp-size <maxSp> -i <fasta file> -t <tx rank> -l <lineage file> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--lineage-file, -l>
  Lineage file.

=item B<--fasta-file, -i>
  Fasta file.

=item B<--tx-rank, -t>
  Taxonomic rank. Possible values: domain, phylum, class, order, family, genus, species.

=item B<--max-sp-size, -s>

  Maximal allowed number of sequences within a species. If a species has more
  sequences, its sequences are aligned, trimmed at the 95% majority start/end
  positions, dereplicated. Then phylo tree is produced and its structure
  investigated by eye.

=item B<--output-dir, -o>
  Output dir.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--quiet>
  Do not print progress messages.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir

  taxonomic_reduction.pl --max-group-size 10000 --max-sp-size 1000 -i Actinobacteria_nr.fa -l ../rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.lineage -o Actinobacteria_dir

  OR

  taxonomic_reduction.pl --max-group-size 10000 --max-sp-size 1000 -i Actinobacteria_nr.fa -o Actinobacteria_dir

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

# my $masterLineageFile = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_nr2.lineage";
# my $masterSeqLenFile  = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_nr2.seqLen";
# my $masterFaFile      = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_nr2.fa";

my $masterLineageFile = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis.lineage";
my $masterSeqLenFile  = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.seqLen";
my $masterFaFile      = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB.fa";

##my $thld = 10000;
GetOptions(
  "master-lineage-file|l=s"   => \$masterLineageFile,
  "master-seqLen-file|m=s"    => \$masterSeqLenFile,
  "master-fasta-file|f=s"     => \$masterFaFile,
  "fasta-file|i=s"            => \my $faFile,
  "tx-rank|t=s"               => \my $txRank,
  "max-sp-size|s=i"           => \my $maxSpSize,
  "max-group-size|g=i"        => \my $thld,
  "output-dir|o=s"            => \my $outDir,
  "skip-FastTree"             => \my $skipFastTree,
  "skip-mafft-sp"             => \my $skipMafftSp,
  "skip-usearch-sp"           => \my $skipUsearchSp,
  "igs"                       => \my $igs,
  "verbose|v"                 => \my $verbose,
  "quiet"                     => \my $quiet,
  "debug"                     => \my $debug,
  "dry-run"                   => \my $dryRun,
  "help|h!"                   => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$faFile)
{
  print "\n\nERROR: Missing fasta file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$maxSpSize)
{
  print "\n\nERROR: Missing sample size\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$outDir)
{
  print "\n\nERROR: Missing output dir\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -f $faFile )
{
  print "\n\nERROR: $faFile does not exist\n\n\n";
  exit 1;
}

my $mothur = "/Users/pgajer/bin/mothur";

if ( defined $igs )
{
  $mothur = "/usr/local/packages/mothur-1.36.1/mothur";
}

my $tmpDir = "temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

####################################################################
##                               MAIN
####################################################################

my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

my $perc = 0;
my $fileSize = -s $faFile; # file size in bytes

print "--- Extracting seq IDs from fasta file\n";
#my %seqIDs; # input fasta file's seq IDs
my @seqIDs; # input fasta file's seq IDs
my %seqLen; # sequence length table
open (IN, "<$faFile") or die "Cannot open $faFile for reading: $OS_ERROR\n";
$/ = ">";
my $junkFirstOne = <IN>;
my $count = 1;
while (<IN>)
{
  if ( !$quiet && ($count % 500 == 0) )
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

    my $perc = sprintf("%.1f%%", 100 * (tell IN) / $fileSize);
    print "\r$timeStr [$perc]";
  }

  chomp;
  my ($id,@seqlines) = split /\n/, $_;
  my $seq = join '', @seqlines;
  $seqLen{$id} = length($seq);
  ##my @seqA = split '', $seq;
  push @seqIDs, $id;
}
close IN;
$/ = "\n";

print "--- Parsing master seqLen table\n";
my %masterSeqLen = read2colTbl($masterSeqLenFile);


print "--- Parsing master lineage table\n";
my %masterLineageTbl = read2colTbl($masterLineageFile);

my %masterSpIdxTbl;
my %masterChildren;
my %masterSpTbl;
my %masterPhTbl;
for my $id ( keys %masterLineageTbl )
{
  my $lineage = $masterLineageTbl{$id};

  my @f = split ";", $lineage;
  ##shift @f;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  my $or = pop @f;
  my $cl = pop @f;
  my $ph = pop @f;

  #$sp = "s_$sp";
  $ge = "g_$ge";
  $fa = "f_$fa";
  $or = "o_$or";
  $cl = "c_$cl";
  $ph = "p_$ph";

  $masterSpIdxTbl{$id} = $sp;

  $masterChildren{"d_Bacteria"}{$ph}++;
  $masterChildren{$ph}{$cl}++;
  $masterChildren{$cl}{$or}++;
  $masterChildren{$or}{$fa}++;
  $masterChildren{$fa}{$ge}++;
  $masterChildren{$ge}{$sp}++;

  push @{$masterSpTbl{$sp}}, $id;
  push @{$masterPhTbl{$ph}}, $id;
}


my %lineageTbl;
@lineageTbl{@seqIDs} = @masterLineageTbl{@seqIDs};


print "\r--- Extracting sizes of  full taxonomy table\n";

# print "Size of seqIDs: " . @seqIDs . "\n";
# print "Size of masterLineageTbl: " . (values %masterLineageTbl) . "\n";
# print "Size of lineageTbl: " . (keys %lineageTbl) . "\n";

# writeTbl(\%lineageTbl, "lineageTbl.txt");
# print "Wrote lineageTbl to lineageTbl.txt\n";
# printTbl(\%masterLineageTbl,"masterLineageTbl");
# writeTbl(\%masterLineageTbl, "masterLineageTbl.txt");

my %spTbl;
my %spIdxTbl;
my %geTbl;
my %faTbl;
my %orTbl;
my %clTbl;
my %phTbl;
my %children;
my %parent;
for my $id ( keys %lineageTbl )
{
  my $lineage = $lineageTbl{$id};

  my @f = split ";", $lineage;
  ##shift @f;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  my $or = pop @f;
  my $cl = pop @f;
  my $ph = pop @f;

  #$sp = "s_$sp";
  $ge = "g_$ge";
  $fa = "f_$fa";
  $or = "o_$or";
  $cl = "c_$cl";
  $ph = "p_$ph";

  $children{"d_Bacteria"}{$ph}++;
  $children{$ph}{$cl}++;
  $children{$cl}{$or}++;
  $children{$or}{$fa}++;
  $children{$fa}{$ge}++;
  $children{$ge}{$sp}++;

  $spIdxTbl{$id} = $sp;

  $parent{$sp} = $ge;
  $parent{$ge} = $fa;
  $parent{$fa} = $or;
  $parent{$or} = $cl;
  $parent{$cl} = $ph;
  $parent{$ph} = "d_Bacteria";

  push @{$spTbl{$sp}}, $id;
  push @{$geTbl{$ge}}, $id;
  push @{$faTbl{$fa}}, $id;
  push @{$orTbl{$or}}, $id;
  push @{$clTbl{$cl}}, $id;
  push @{$phTbl{$ph}}, $id;
}

# print "Sum over species: ";
# my $spSum = 0;
# map{ $spSum += @{$spTbl{$_}} } keys %spTbl;
# print "$spSum\n\n";

# print "Sum over genera: ";
# my $geSum = 0;
# map{ $geSum += @{$geTbl{$_}} } keys %geTbl;
# print "$geSum\n\n";

# Checking the number of sequences within each species. If a species has more
# than $maxSpSize sequences, its sequences are aligned, trimmed at the 95%
# majority start/end positions, dereplicated. Then phylo tree is produced and its
# structure investigated by eye.

##my $cmd = "rm -rf $outDir; mkdir $outDir";
my $cmd = "mkdir -p $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

chdir $outDir;
$faFile = "../$faFile";

my $logFile = "logFile.txt";
open LOGOUT, ">$logFile" or die "Cannot open $logFile for writing: $OS_ERROR\n";

print "--- Subsampling species with more than $maxSpSize sequences\n";

##my %ssSpTbl = %spTbl;
my %spDiff;
my %geDiff;
my %faDiff;
my %orDiff;
my %clDiff;
my %phDiff;

my @ssSpSize; # sizes of subsampled species

for my $sp ( sort{ scalar(@{$spTbl{$a}}) <=> scalar(@{$spTbl{$b}}) } keys %spTbl )
{
  my $nSpSeqs = scalar(@{$spTbl{$sp}});

  if ( $nSpSeqs > $maxSpSize )
  {
    print "\nProcessing $sp\tsize: $nSpSeqs\n";

    print LOGOUT "$sp size: $nSpSeqs\n";

    ## create $sp dir
    my $spDir = $sp . "_dir";
    # if ( -e $spDir )
    my $cmd = "mkdir -p $spDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    # produce hash table of sequences of the given species
    my %selSp;
    @selSp{@{$spTbl{$sp}}} = (1) x $nSpSeqs;
    #print "\nselSp size:" . (keys %selSp) . "\n";

    my $spFile = "$spDir/$sp" . ".fa";
    print "\t--- Producing $spFile of $sp seq's\n";
    if (!$skipMafftSp)
    {
      open OUT, ">$spFile" or die "Cannot open $spFile for writing: $OS_ERROR\n";
      open (IN, "<$faFile") or die "Cannot open $faFile for reading: $OS_ERROR\n";
      $/ = ">";
      my $junkFirstOne = <IN>;
      my $count = 1;
      while (<IN>)
      {
	if ( !$quiet && ($count % 500 == 0) )
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
	  print "\r$timeStr";
	}

	chomp;
	my ($id,@seqlines) = split /\n/, $_;
	my $seq = join '', @seqlines;

	if ( defined $selSp{$id} )
	{
	  print OUT ">$id\n";
	  print OUT "$seq\n";
	}
	$count++;
      }
      close IN;
      $/ = "\n";
      close OUT;
    }
    #print "\r                     ";
    print "\r\t--- Screen $spFile\n";

    my @tmp;
    push (@tmp,"summary.seqs(fasta=$spFile)");

    my $spSummaryFile = "$spDir/$sp" . ".summary";

    my $criteria = "90";
    # if ( 0.1 * $nSpSeqs > $maxSpSize )
    # {
    #   $criteria = "9";
    # }
    # elsif ( 0.5 * $nSpSeqs < $maxSpSize )
    # {
    #   $criteria = "49";
    # }
    # else
    # {
    #   $criteria = sprintf("%d",int(10 * $maxSpSize / $nSpSeqs));
    # }
    # print "\n===> criteria: $criteria\n";

    push (@tmp,"screen.seqs(fasta=$spFile, summary=$spSummaryFile, maxambig=0, optimize=minlength, criteria=$criteria)");

    my $spGoodFile = "$spDir/$sp" . ".good.fa";;

    push (@tmp,"summary.seqs(fasta=$spGoodFile)");

    printArray(\@tmp, "mothur commands") if ($debug || $verbose);

    my $scriptFile = create_mothur_script( \@tmp );
    $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
    #$cmd = "rm -f mothur.log; $mothur < $scriptFile > mothur.log 2>&1; rm -f $scriptFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

    print "\t--- Reading $spGoodFile to update array of seq IDs from $sp\n";
    my @nrSp;
    open (IN, "<$spGoodFile") or die "Cannot open $spGoodFile for reading: $OS_ERROR\n";
    $/ = ">";
    $junkFirstOne = <IN>;
    #$count = 1;
    while (<IN>)
    {
      # if ( !$quiet && ($count % 500 == 0) )
      # {
      # 	$endRun = time();
      # 	$runTime = $endRun - $startRun;
      # 	if ( $runTime > 60 )
      # 	{
      # 	  $timeMin = int($runTime / 60);
      # 	  $timeSec = sprintf("%02d", $runTime % 60);
      # 	  $timeStr = "$timeMin:$timeSec";
      # 	}
      # 	else
      # 	{
      # 	  $runTime = sprintf("%02d", $runTime);
      # 	  $timeStr = "$timeMin:$runTime";
      # 	}
      # 	print "\r$timeStr";
      # }

      chomp;
      my ($id,@seqlines) = split /\n/, $_;
      push @nrSp, $id;
      #$count++;
    }
    close IN;
    $/ = "\n";
    #print "\r                                ";

    print LOGOUT "size(nrSp): " . @nrSp ."\n";


    ##@{$ssSpTbl{$sp}} = @nrSp;
    ##@{$spTbl{$sp}} = @nrSp;

    my @ssSp;
    if ( @nrSp > $maxSpSize )
    {
      fisher_yates_shuffle(\@nrSp);
      @nrSp = @nrSp[0 .. ($maxSpSize-1)];
      #@{$spTbl{$sp}} = @{$spTbl{$sp}}[0 .. ($maxSpSize-1)];

      print "\t--- Selecting seq's of nrSp from $spGoodFile\n";
      my $selSeqsFile = "$spDir/$sp" . "_ss.SeqIDs";
      open OUT, ">$selSeqsFile" or die "Cannot open $selSeqsFile for writing: $OS_ERROR\n";
      for (@nrSp)
      {
	print OUT "$_\n";
      }
      close OUT;

      my $spSelFile = "$spDir/$sp" . "_ss.fa";
      $cmd = "select_seqs.pl --quiet -s $selSeqsFile -i $spGoodFile -o $spSelFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;


      print "\t--- Generating MSA of $sp seq's\n";
      my $spAlgnFile = "$spDir/$sp" . "_algn.fa";
      $cmd = "mafft --auto --inputorder --quiet --thread 4 $spSelFile > $spAlgnFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      print "\t--- QCing and trimming alignment file\n";
      my @tmp;
      push (@tmp,"summary.seqs(fasta=$spAlgnFile)");
      printArray(\@tmp, "mothur commands") if ($debug || $verbose);

      my $scriptFile = create_mothur_script( \@tmp );
      $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      print "\t--- Trimming $spAlgnFile\n";
      my $spAlgnSummaryFile = "$spDir/$sp" . "_algn.summary";
      my $spGoodAlgnFile = "$spDir/$sp" . "_algn.good.fa";;
      my $spTrAlgnFile = "$spDir/$sp" . "_algn_trimmed.fa";
      $cmd = "trim_align.pl --quiet -c 95 -j $spAlgnSummaryFile -i $spAlgnFile -o $spTrAlgnFile";
      ##$cmd = "trim_align.pl -c 95 -j $spAlgnSummaryFile -i $spGoodAlgnFile -o $spTrAlgnFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      print "\t--- Removing gaps from $spTrAlgnFile\n";
      my $spTrFile = "$spDir/$sp" . "_trimmed.fa";
      $cmd = "rmGaps --quiet -i $spTrAlgnFile -o $spTrFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      print "\t--- Dereplicating $spTrFile\n";
      my $spNrFile = "$spDir/$sp" . "_nr.fa";
      if ($igs)
      {
	$cmd = "usearch7 -quiet -cluster_fast $spTrFile -id 1.0 -centroids $spNrFile";
      }
      else
      {
	$cmd = "usearch8 -quiet -cluster_fast $spTrFile -id 1.0 -centroids $spNrFile";
      }
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      print "\t--- Selecting nr seq's from the alignment\n";

      my $seqIDsNrFile = "$spDir/$sp" . "_nr.seqIDs";
      open OUT, ">$seqIDsNrFile" or die "Cannot open $seqIDsNrFile for writing: $OS_ERROR\n";
      open (IN, "<$spNrFile") or die "Cannot open $spNrFile for reading: $OS_ERROR\n";
      $/ = ">";
      my $junkFirstOne = <IN>;
      #my $count = 1;
      while (<IN>)
      {
	# if ( !$quiet && ($count % 500 == 0) )
	# {
	#   $endRun = time();
	#   $runTime = $endRun - $startRun;
	#   if ( $runTime > 60 )
	#   {
	#     $timeMin = int($runTime / 60);
	#     $timeSec = sprintf("%02d", $runTime % 60);
	#     $timeStr = "$timeMin:$timeSec";
	#   }
	#   else
	#   {
	#     $runTime = sprintf("%02d", $runTime);
	#     $timeStr = "$timeMin:$runTime";
	#   }
	#   my $perc = sprintf("%.1f%%", 100 * (tell IN) / $fileSize);
	#   print "\r$timeStr [$perc]";
	# }

	chomp;
	my ($id,@seqlines) = split /\n/, $_;
	print OUT "$id\n";
      }
      close IN;
      $/ = "\n";
      close OUT;
      # if ( !$quiet )
      # {
      # 	print "\r            ";
      # }

      # my $spSeqIDsFile = "$spDir/$sp" . "_nr.seqIDs";
      # open OUT, ">$spSeqIDsFile" or die "Cannot open $spSeqIDsFile for writing: $OS_ERROR\n";
      # for (@nrSp)
      # {
      # 	print OUT "$_\n";
      # }
      # close OUT;

      my $spTrNrAlgnFile = "$spDir/$sp" . "_algn_nr.fa";
      $cmd = "select_seqs.pl --quiet -s $seqIDsNrFile -i $spTrAlgnFile -o $spTrNrAlgnFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;


      #####  Computing the distance between sequences
      push (@tmp,"dist.seqs(fasta=$spTrNrAlgnFile, calc=onegap, output=lt)");
      ##push (@tmp,"dist.seqs(fasta=$spTrNrAlgnFile, calc=onegap)");
      printArray(\@tmp, "mothur commands") if ($debug || $verbose);

      my $scriptFile = create_mothur_script( \@tmp );
      $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;


      ##### Translate phylip dist. matrix into a full dist. matrix
      print "\t--- Translating into full distance matrix\n";
      my $phylipDistFile = "$spDir/$sp" . "_algn_nr.phylip.dist";
      my $fullDistFile = "$spDir/$sp" . "_algn_nr.dist";
      $cmd = "lowerTriaMatToFull.pl --quiet -i $phylipDistFile -o $fullDistFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      ##### Generate Ward linkage hierarchical clustering
      ##### Use Silhouette Width to determine the number of clusters
      ##### Remove small clusters
      print "\t--- Extracting seqIDs of largest Ward linkage hierarchical clustering clusters comprising 90% of data\n";
      my $selSSseqsFile = "$spDir/$sp" . "_algn_nr_ss.seqIDs";
      $cmd = "distMat_to_selSeqs.pl --quiet -m ward.D2 -p 90 -i $fullDistFile -o $selSSseqsFile";
      print "cmd=$cmd\n" if $dryRun;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$skipMafftSp;

      print "\t--- Extracting selected seq IDs from $selSSseqsFile\n";
      open IN, "<$selSSseqsFile" or die "Cannot open $selSSseqsFile for reading: $OS_ERROR\n";
      for (<IN>)
      {
	chomp;
	push @ssSp, $_;
      }
      close IN;

      if (0)
      {
	print "Test ssSp:\n";
	my $i = 0;
	for (@ssSp)
	{
	  print "$_\n";
	  $i++;
	  if ($i>10)
	  {
	    exit 1;
	  }
	}
      }
    }
    else
    {
      @ssSp = @nrSp;
    }

    push @ssSpSize, scalar(@ssSp);
    print "\nSize after subsampling: " . scalar(@ssSp) . "\n";
    ##print "\nsize(ssSp): " . scalar(@ssSp) ."\n";
    print LOGOUT "size(ssSp): " . scalar(@ssSp) ."\n";

    ## updating family and other taxon levels
    my @spDiff = diff($spTbl{$sp}, \@ssSp);
    @{$spTbl{$sp}} = @ssSp;
    $spDiff{$sp} += @spDiff;

    my $ge = $parent{$sp};
    @{$geTbl{$ge}} = diff($geTbl{$ge}, \@spDiff);
    $geDiff{$ge} += @spDiff;

    my $fa = $parent{$ge};
    @{$faTbl{$fa}} = diff($faTbl{$fa}, \@spDiff);
    $faDiff{$fa} += @spDiff;

    my $or = $parent{$fa};
    @{$orTbl{$or}} = diff($orTbl{$or}, \@spDiff);
    $orDiff{$or} += @spDiff;

    my $cl = $parent{$or};
    @{$clTbl{$cl}} = diff($clTbl{$cl}, \@spDiff);
    $clDiff{$cl} += @spDiff;

    my $ph = $parent{$cl};
    @{$phTbl{$ph}} = diff($phTbl{$ph}, \@spDiff);
    $phDiff{$ph} += @spDiff;

  } # end of  if ( $nSpSeqs > $maxSpSize )
}


close LOGOUT;

if (0)
{
  print "Phyla\n";
  my @phs = keys %phTbl;
  printFormatedTbl(\%phTbl, \@phs);

  print "Classes\n";
  my @cls = sort{scalar(@{$clTbl{$b}}) <=> scalar(@{$clTbl{$a}})} keys %clTbl;
  printFormatedTbl(\%clTbl, \@cls);

  print "Orders\n";
  my @ors = sort{scalar(@{$orTbl{$b}}) <=> scalar(@{$orTbl{$a}})} keys %orTbl;
  printFormatedTbl(\%orTbl, \@ors);

  print "Families\n";
  my @fas = sort{scalar(@{$faTbl{$b}}) <=> scalar(@{$faTbl{$a}})} keys %faTbl;
  printFormatedTbl(\%faTbl, \@fas);

  print "\n\n";
}


# Print all taxonomic ranks up to family in the nested way with size of each
# Report if any family has size greater than the group size threshold !
# Stop if it is the case.

print "\nTaxonomic ranks present in the data\n";
my %bigFa; # table of families with size greater than maxGr
my @phs = keys %phTbl;
#printFormatedTbl(\%phTbl, \@phs);
my $phKVL = getKeyValStrLengths(\%phTbl);
for my $ph (@phs)
{
  printF(0, $ph, scalar(@{$phTbl{$ph}}), $phKVL);
  my @cls = keys %{$children{$ph}};
  my $clKVL = getKeyValStrLengths(\%clTbl, \@cls);
  for my $cl ( sort{scalar(@{$clTbl{$a}}) <=> scalar(@{$clTbl{$b}})} @cls)
  {
    printF(1, $cl, scalar(@{$clTbl{$cl}}), $clKVL);
    my @ors = keys %{$children{$cl}};
    my $orKVL = getKeyValStrLengths(\%orTbl, \@ors);
    for my $or ( sort{scalar(@{$orTbl{$a}}) <=> scalar(@{$orTbl{$b}})} @ors)
    {
      printF(2, $or, scalar(@{$orTbl{$or}}), $orKVL);
      my @fas = keys %{$children{$or}};
      my $faKVL = getKeyValStrLengths(\%faTbl, \@fas);
      for my $fa ( sort{scalar(@{$faTbl{$a}}) <=> scalar(@{$faTbl{$b}})} @fas)
      {
	printF(3, $fa, scalar(@{$faTbl{$fa}}), $faKVL);
	if ( scalar(@{$faTbl{$fa}}) > $thld)
	{
	  $bigFa{$fa} = scalar(@{$faTbl{$fa}});
	}
      }
    }
  }
}

if ( keys %bigFa > 0 )
{
  print "WARNING: Detected families of size exceeding $thld\n";
  my @fas =  keys %bigFa;
  my $faKVL = getKeyValStrLengths(\%faTbl, \@fas);
  for my $fa ( sort{ $bigFa{$b} <=> $bigFa{$a} } @fas)
  {
    printF(1, $fa, $bigFa{$fa}, $faKVL);
    print "Applying hard-wired upper limit for species size of 240 (see ssSpsizes.R)\n";

    my @spp = get_spp_of_fa($fa);
    ##my @sppSize; # sizes of species of the given family
    $maxSpSize = 240;
    for my $sp (@spp)
    {
      my $nSpSeqs = scalar(@{$spTbl{$sp}});
      if ( $nSpSeqs > $maxSpSize )
      {
	fisher_yates_shuffle($spTbl{$sp});
	my @newIDs = @{$spTbl{$sp}}[0 .. ($maxSpSize-1)];
	my @spDiff = diff($spTbl{$sp}, \@newIDs);
	@{$spTbl{$sp}} = @newIDs;

	my $ge = $parent{$sp};
	@{$geTbl{$ge}} = diff($geTbl{$ge}, \@spDiff);
	$geDiff{$ge} += @spDiff;

	my $fa = $parent{$ge};
	@{$faTbl{$fa}} = diff($faTbl{$fa}, \@spDiff);
	$faDiff{$fa} += @spDiff;

	my $or = $parent{$fa};
	@{$orTbl{$or}} = diff($orTbl{$or}, \@spDiff);
	$orDiff{$or} += @spDiff;

	my $cl = $parent{$or};
	@{$clTbl{$cl}} = diff($clTbl{$cl}, \@spDiff);
	$clDiff{$cl} += @spDiff;

	my $ph = $parent{$cl};
	@{$phTbl{$ph}} = diff($phTbl{$ph}, \@spDiff);
	$phDiff{$ph} += @spDiff;
      }
    }

    if ( scalar(@{$faTbl{$fa}}) > $thld )
    {
      print "\nERROR: after applying 240 species size limit $fa is still too big with size: " . scalar(@{$faTbl{$fa}}) . "\n";
      exit 1;
    }
    else
    {
     print "\nAfter applying 240 species size limit to $fa its size is " . scalar(@{$faTbl{$fa}}) . "\n";
    }
  }

  # when you run it with the default settings its going to stop after species
  # subsampling as the group size is too small. The R script I want to write
  # suppose predict max family size given different values of max species
  # size. This script should be run automatically when group size criterion is
  # violated and a plot should be generated with max family size vs. max species
  # size. Looking at the plot, one can choose both the group size and max species
  # sizes

  # my $sppSizesFile = "ssSppSizes.txt";
  # print "\rWriting species sizes after subsampling to $sppSizesFile";
  # open OUT, ">$sppSizesFile" or die "Cannot open $sppSizesFile for writing: $OS_ERROR\n";
  # for my $sp ( sort{scalar(@{$spTbl{$b}}) <=> scalar(@{$spTbl{$a}})} keys %spTbl)
  # {
  #   print OUT "$sp\t" . scalar(@{$spTbl{$sp}}) . "\n";
  # }
  # close OUT;

  # my $ssOnlySppSizeFile = "ssOnlySppSizes.txt";
  # print "\rWriting subsampled only species sizes to $ssOnlySppSizeFile";
  # writeArray(\@ssSpSize, $ssOnlySppSizeFile);

  # my $parentFile = "parentTbl.txt";
  # print "\rWriting parent table to $parentFile          ";
  # writeTbl(\%parent, $parentFile);

  # my $childrenFile = "childrenTbl.txt";
  # print "\rWriting children table to $childrenFile          ";
  # writeChildrenTbl($childrenFile);

  # print "\r                                                                       \n";
  # print "\tSpecies sizes after subsampling written to $outDir/$sppSizesFile\n";
  # print "\tSubsampled only species sizes written to $outDir/$ssOnlySppSizeFile\n";
  # print "\tParent table written to $outDir/$parentFile\n\n";

  # exit 1;
}


if (0)
{
  my $spTblFile = "spTbl.csv";
  writeArrayTbl(\%spTbl, $spTblFile);

  my $seqLenFile = "seqLen.csv";
  writeTblCSV(\%seqLen, $seqLenFile);
}

my $sppSizesFile = "ssSppSizes.txt";
print "\rWriting species sizes after subsampling to $sppSizesFile";
open OUT, ">$sppSizesFile" or die "Cannot open $sppSizesFile for writing: $OS_ERROR\n";
for my $sp ( sort{scalar(@{$spTbl{$b}}) <=> scalar(@{$spTbl{$a}})} keys %spTbl)
{
  print OUT "$sp\t" . scalar(@{$spTbl{$sp}}) . "\n";
}
close OUT;

my $ssOnlySppSizeFile = "ssOnlySppSizes.txt";
print "\rWriting subsampled only species sizes to $ssOnlySppSizeFile";
writeArray(\@ssSpSize, $ssOnlySppSizeFile);

my $parentFile = "parentTbl.txt";
print "\rWriting parent table to $parentFile          ";
writeTbl(\%parent, $parentFile);

my $childrenFile = "childrenTbl.csv";
print "\rWriting children table to $childrenFile          ";
writeChildrenTbl($childrenFile);

print "\r                                                                       \n";
print "\tSpecies sizes after subsampling written to $outDir/$sppSizesFile\n";
print "\tSubsampled only species sizes written to $outDir/$ssOnlySppSizeFile\n";
print "\tParent table written to $outDir/$parentFile\n\n";


print "\n\nSplitting subsequenced data into groups\n";
print "Reporting taxon size and group to which it belongs\n";

my %group;   # @{$group{$idx}}: array of seqID's of group $idx
my %groupTx; # taxonomic ranks of the given group
my $idx = 0;
## here we assume that the input fasta file consists of sequences of a single phylum
my @cls = sort{scalar(@{$clTbl{$b}}) <=> scalar(@{$clTbl{$a}})} keys %clTbl;
for my $cl ( reverse @cls)
{
  ##print "Processing $cl\t" . scalar(@{$clTbl{$cl}}) . "\n";
  print "$cl\t" . scalar(@{$clTbl{$cl}}) . "\t";
  if ( scalar(@{$clTbl{$cl}}) < $thld )
  {
    if (defined $group{$idx} )
    {
      if ( (scalar(@{$group{$idx}}) + scalar(@{$clTbl{$cl}})) < $thld )
      {
	push @{$group{$idx}}, @{$clTbl{$cl}};
      }
      else
      {
	$idx++;
	push @{$group{$idx}}, @{$clTbl{$cl}};
      }
    }
    else
    {
      push @{$group{$idx}}, @{$clTbl{$cl}};
    }
    if ($debug)
    {
      print "Added $cl to group $idx\n";
      print "Added $parent{$cl} to group $idx\n";
      printFormatedTbl(\%group);
    }
    push @{$groupTx{$idx}}, $cl;
    push @{$groupTx{$idx}}, $parent{$cl};
    print "$idx\n";
  }
  else
  {
    print "\n";
    my @ors = sort{scalar(@{$orTbl{$a}}) <=> scalar(@{$orTbl{$b}})} keys %{$children{$cl}};
    if ($debug)
    {
      print "\tchildren{$cl}:\n";
      map{ print "\t\t$_\t" . scalar(@{$orTbl{$_}}) . "\n" } @ors;
      print "\n";
    }
    for my $or ( @ors )
    {
      ##print "\tProcessing $or\t" . scalar(@{$orTbl{$or}}) . "\n";
      print "\t$or\t" . scalar(@{$orTbl{$or}}) . "\t";

      if ( scalar(@{$orTbl{$or}}) < $thld )
      {
	if (defined $group{$idx})
	{
	  if ( (scalar(@{$group{$idx}}) + scalar(@{$orTbl{$or}})) < $thld )
	  {
	    push @{$group{$idx}}, @{$orTbl{$or}};
	  }
	  else
	  {
	    $idx++;
	    push @{$group{$idx}}, @{$orTbl{$or}};
	  }
	}
	else
	{
	  push @{$group{$idx}}, @{$orTbl{$or}};
	}
	if ($debug)
	{
	  print "Added $or to group $idx\n";
	  print "group\n";
	  printFormatedTbl(\%group);
	}
	push @{$groupTx{$idx}}, $or;
	push @{$groupTx{$idx}}, $parent{$or};
	push @{$groupTx{$idx}}, $parent{$parent{$or}};
	print "$idx\n";
      }
      else
      {
	print "\n";
	my @fas = sort{scalar(@{$faTbl{$a}}) <=> scalar(@{$faTbl{$b}})} keys %{$children{$or}};
	if ($debug)
	{
	  print "\t\tchildren{$or}:\n";
	  map{ print "\t\t\t$_\t" . scalar(@{$faTbl{$_}}) . "\n" } @fas;
	  print "\n";
	}

	for my $fa ( @fas )
	##for my $fa ( keys %{$children{$or}} )
	{
	  print "\t\t$fa\t" . scalar(@{$faTbl{$fa}}) . "\t";

	  if ( scalar(@{$faTbl{$fa}}) < $thld )
	  {
	    if (defined $group{$idx} )
	    {
	      if ( (scalar(@{$group{$idx}}) + scalar(@{$faTbl{$fa}})) < $thld )
	      {
		push @{$group{$idx}}, @{$faTbl{$fa}};
	      }
	      else
	      {
		$idx++;
		push @{$group{$idx}}, @{$faTbl{$fa}};
	      }
	    }
	    else
	    {
	      push @{$group{$idx}}, @{$faTbl{$fa}};
	    }
	    if ($debug)
	    {
	      print "Added $fa to group $idx\n";
	      print "group\n";
	      printFormatedTbl(\%group);
	    }
	    push @{$groupTx{$idx}}, $fa;
	    push @{$groupTx{$idx}}, $parent{$fa};
	    push @{$groupTx{$idx}}, $parent{$parent{$fa}};
	    push @{$groupTx{$idx}}, $parent{$parent{$parent{$fa}}};
	    print "$idx\n";
	  }
	  else
	  {
	    print "ERROR: $fa exceeds $thld size: " . scalar(@{$faTbl{$fa}}) . "\n";
	    exit 1;
	  }
	}
      }
    }
  }
}

print "\n--- Removing redundant taxons from groupTx{idx}\n";
for my $idx ( keys %groupTx )
{
  my @u = unique($groupTx{$idx});
  @{$groupTx{$idx}} = @u;
}

## print "\n\nFinal idx: $idx\n";
print "\n\nFinal groups:\n";
my @grIdx = sort { $a <=> $b } keys %group;
printFormatedTbl(\%group, \@grIdx);


print "\n\nTaxonomic ranks of each group\n";
for my $idx ( sort {$a <=> $b} keys %groupTx)
{
  print "\nGroup $idx\n";
  for my $tx ( sort @{$groupTx{$idx}} )
  {
    my $tabs = ""; # no tabs for classes
    my @f = split "_", $tx;
    my $prefix = shift @f;
    if ( $prefix eq "o" )
    {
      $tabs = "\t";
    }
    elsif ( $prefix eq "f" )
    {
      $tabs = "\t\t";
    }
    print "$tabs$tx\n";
  }
}

my $grTxFile = "groupTx.cvs";
writeArrayTbl(\%groupTx, $grTxFile);
print "\n\tGroup taxonomic ranks written to $grTxFile\n";

## Finding outgroup seqIDs for each group
## They will be added to fa and tx files of each group

print "\n\n--- Finding outgroup seq's for each group\n";

my %useMasterFasta;

my %outgroup;
my %txInt;
$txInt{'p'} = 0;
$txInt{'c'} = 1;
$txInt{'o'} = 2;
$txInt{'f'} = 3;
$txInt{'g'} = 4;
$txInt{'s'} = 5;
for my $idx ( sort {$a <=> $b} keys %groupTx)
{
  print "\rGroup $idx";

  my $debugOG = 1;
  my %grTx; # group taxons
  my @taxons = @{$groupTx{$idx}};
  for my $tx ( @taxons )
  {
    my $c = get_tx_prefix($tx);
    my $i = $txInt{$c};
    push @{$grTx{$c}}, $tx;
  }

  printArrayTbl(\%grTx, "\ngrTx") if $debugOG;

  ## pick the lowest taxon with the single element
  my @cc = ('s','g','f','o','c','p');
  my $cmin;
  for my $c (@cc)
  {
    if ( exists $grTx{$c} && @{$grTx{$c}}==1 )
    {
      $cmin = $c;
      last;
    }
  }
  my @t;
  my $gp; # grandparent
  my $p;  # parent
  if ( exists $grTx{$cmin} )
  {
    @t  = @{$grTx{$cmin}}; # extracting highest level singleton taxon
    $p  = $t[0];
    $gp = $parent{$t[0]};
  }
  else
  {
    print "\nERROR: grTx{$cmin} does not exist\n";
    exit 1;
  }

  if ($gp eq "d_Bacteria")
  {
    $useMasterFasta{$idx} = 1;
  }
  else
  {
    $useMasterFasta{$idx} = 0;
  }

  printArray(\@t, "t=grTx{cmin}") if $debugOG;

  my @ch = keys %{$masterChildren{$gp}};

  printArray(\@ch, "ch") if $debugOG;

  @t = ($p);
  my @sibs = diff(\@ch, \@t);

  printArray(\@sibs, "sibs") if $debugOG;
  print "\n\n" if $debugOG;

  if (@sibs == 0)
  {
    print "\nERROR: Difference between ch and t is zero\n";
    print "ch: @ch\n";
    print "sibs: @sibs\n";
    exit 1;
  }

  my @seqIDs;
  my $perc = 90; # randomly selected seq's should have lengths in the upper 90% percentile
  for my $sib (@sibs)
  {
    ##my $id = get_one_rand_sp_seq($tx, \%spTbl, \%seqLen, $perc);
    ## selecting only one sequences for each outgroup taxon
    my @a;
    if ($cmin eq 'p')
    {
      @a = @{$masterPhTbl{$sib}};
    }
    elsif ($cmin eq 'c')
    {
      @a = @{$clTbl{$sib}};
    }
    elsif ($cmin eq 'o')
    {
      @a = @{$orTbl{$sib}};
    }
    elsif ($cmin eq 'f')
    {
      @a = @{$faTbl{$sib}};
    }
    my $id = $a[rand @a];
    push @seqIDs, $id;
  }

  push @{$outgroup{$idx}}, @seqIDs;
}
print "\r                                                                           \n";

print "\n\n--- Producing fa, algn, tree, tx and outgroup files of each group\n";
my ($dirPrefix, $str) = split "_", $outDir;
for my $idx (keys %group)
{
  print "\n\nProcessing Group $idx\n";

  my $grPrefix = "$dirPrefix" . "_group_$idx";
  my $grDir = $grPrefix . "_dir";

  my $cmd = "mkdir -p $grDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  ##chdir $outDir;

  $grPrefix = "$grDir/$grPrefix";

  my $grFaFile          = $grPrefix . ".fa";
  my $grAlgnFile        = $grPrefix . "_algn.fa";
  my $grLineageFile     = $grPrefix . ".lineage";
  my $grTxFile          = $grPrefix . ".tx";
  my $grTxSummaryFile   = $grPrefix . ".txSummary";
  my $grOutgroupFile    = $grPrefix . "_outgroup.seqIDs";
  my $grOgFaFile        = $grPrefix . "_og.fa";
  my $grUrTreeFile      = $grPrefix . "_unrooted.tree";
  my $grTreeFile        = $grPrefix . ".tree";
  my $sppSeqIDsTreeFile = $grPrefix . "_sppSeqIDs.tree";

  my @ogSeqIDs = @{$outgroup{$idx}};

  ##print "\rCreating seqIDs file of the give group's seqIDs plus its outgroup sequences";
  my $seqIDsGrFile = "$grPrefix" . ".seqIDs";
  open OUT, ">$seqIDsGrFile" or die "Cannot open $seqIDsGrFile for writing: $OS_ERROR\n";
  for ( @{$group{$idx}} )
  {
    print OUT "$_\n";
  }
  for (@ogSeqIDs)
  {
    print OUT "$_\n";
  }
  close OUT;

  # Creating $grTxFile
  open OUT, ">$grTxFile" or die "Cannot open $grTxFile for writing: $OS_ERROR\n";
  for ( @{$group{$idx}} )
  {
    my $sp = $spIdxTbl{$_};
    # my @f = split "_", $sp;
    # shift @f;
    # $sp = join "_", @f;
    print OUT "$_\t$sp\n";
    ## if outgroup than use OUTGROUP_i label
  }
  ##my $oCount = 1;
  for (@{$outgroup{$idx}})
  {
    ##print OUT "$_\tOUTGROUP_$oCount\n";
    my $sp = $masterSpIdxTbl{$_};
    print OUT "$_\t$sp" . "_OG\n";
    ##$oCount++;
  }
  close OUT;

  # Creating taxon summary group file
  open OUT, ">$grTxSummaryFile" or die "Cannot open $grTxSummaryFile for writing: $OS_ERROR\n";
  for ( @{$groupTx{$idx}} )
  {
    print OUT "$_\n";
  }
  close OUT;

  # Creating outgrGrFile
  open OUT, ">$grOutgroupFile" or die "Cannot open $grOutgroupFile for writing: $OS_ERROR\n";
  for (@{$outgroup{$idx}})
  {
    print OUT "$_\n";
  }
  close OUT;

  # Creating $faGrFile
  print "--- Creating $grFaFile\n";
  if ($useMasterFasta{$idx})
  {
    $cmd = "select_seqs.pl --quiet -s $seqIDsGrFile -i $masterFaFile -o $grFaFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }
  else
  {
    $cmd = "select_seqs.pl --quiet -s $seqIDsGrFile -i $faFile -o $grFaFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }

  # Creating $algnGrFile
  print "--- Creating $grAlgnFile\n";
  $cmd = "rm -f $grAlgnFile; mafft --auto --inputorder --quiet --thread 4 $grFaFile > $grAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # trim
  print "--- Generating seq's summary of the alignment file\n";
  my @tmp;
  push (@tmp,"summary.seqs(fasta=$grAlgnFile)");
  printArray(\@tmp, "mothur commands") if ($debug || $verbose);

  my $scriptFile = create_mothur_script( \@tmp );
  $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Trimming alignment\n";
  my $summaryFile     = $grPrefix . "_algn.summary";
  my $trimmedAlgnFile = $grPrefix . "_algn_trimmed.fa";
  $cmd = "trim_align.pl -c 95 -j $summaryFile -i $grAlgnFile -o $trimmedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Generating a phylo tree based on trimmed alignment\n";
  $cmd = "rm -f $grUrTreeFile; FastTree -nt $trimmedAlgnFile > $grUrTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;# && !$skipFastTree;

  print "--- Rerooting the tree using outgroup seq's\n";
  $cmd = "rm -f $grTreeFile; nw_reroot $grUrTreeFile @ogSeqIDs | nw_order -  > $grTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
  my $sppSeqIDsFile  = $grPrefix . ".sppSeqIDs";
  $cmd = "rm -f $sppSeqIDsFile; awk '{print \$1\"\\t\"\$2\"__\"\$1}' $grTxFile > $sppSeqIDsFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $sppSeqIDsTreeFile; nw_rename $grTreeFile $sppSeqIDsFile | nw_order -  > $sppSeqIDsTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


  # Creating lineage file
  open OUT, ">$grLineageFile" or die "Cannot open $grLineageFile for writing: $OS_ERROR\n";
  for ( @{$group{$idx}} )
  {
    print OUT "$_\t" . $lineageTbl{$_} . "\n";
  }
  if ($useMasterFasta{$idx})
  {
    for (@{$outgroup{$idx}})
    {
      print OUT "$_\t" . $masterLineageTbl{$_} . "\n";
    }
  }
  else
  {
    for (@{$outgroup{$idx}})
    {
      print OUT "$_\t" . $lineageTbl{$_} . "\n";
    }
  }
  close OUT;
}
print "\rfa, algn, lineage, tx, tree and outgroup files of each group written to $outDir\n";



my @subsampledIDs;
for (keys %group)
{
  push @subsampledIDs, @{$group{$_}};
}

print "\n\tNumber of subsampled IDs: " . @subsampledIDs . "\n";

## subsampled spp tbl
my %ssSpFreqTbl;
for (@subsampledIDs)
{
  $ssSpFreqTbl{$spIdxTbl{$_}}++;
}

print "\tNumber of species: " . (keys %spTbl) . "\n";
print "\tNumber of subsampled species: " . (keys %ssSpFreqTbl) . "\n";

if (0)
{
  print "\nComparing frequencies of species before and after subsampling\n";
  ##for my $sp ( sort{ $ssSpFreqTbl{$b} <=> $ssSpFreqTbl{$a} } keys %ssSpFreqTbl )
  for my $sp ( keys %spTbl )
  {
    if ( $ssSpFreqTbl{$sp} != scalar(@{$spTbl{$sp}}) )
    {
      print "$sp\t" . @{$spTbl{$sp}} . "\t" . $ssSpFreqTbl{$sp} . "\n";
    }
    else
    {
      print "$sp\t" . @{$spTbl{$sp}} . "\t" . $ssSpFreqTbl{$sp} . "\n";
    }
  }
  print "\n";
}


if (0)
{
  print "Genera\n";
  for (keys %geTbl)
  {
    print "\t$_\t" . @{$geTbl{$_}} . "\n";
  }
  print "\n";

  print "Species\n";
  for (keys %spTbl)
  {
    print "\t$_\t" . @{$spTbl{$_}} . "\n";
  }
  print "\n";
}

print "--- Removing mothur log files\n" if !$quiet;
my $dir = '.';
opendir(DIR, $dir) or die $!;
while (my $file = readdir(DIR))
{
  next if ($file !~ /logfile/);
  unlink $file;

}
closedir(DIR);

print "\tOutput written to $outDir\n\n";

## report timing
$endRun = time();
$runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  $timeMin = int($runTime / 60);
  $timeSec = sprintf("%02d", $runTime % 60);
  print "\r\tCompleted in $timeMin:$timeSec\n"
}
else
{
  print "\r\tCompleted in $runTime seconds\n"
}


####################################################################
##                               SUBS
####################################################################

## get species names of a given family
sub get_spp_of_fa
{
  my $tx = shift;

  my $c = get_tx_prefix($tx);

  my @spp;
  if ($c eq 'f')
  {
    my @ges = keys %{$children{$tx}};
    for my $ge (@ges)
    {
      my @sps = keys %{$children{$ge}};
      push @spp, @sps;
    }
  }
  else
  {
    print "ERROR in get_spp_of_fa(): Input taxon not a family\n";
    exit 1;
  }

  return @spp;
}

# write hash table to a CSV file
sub writeTblCSV
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "," . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# write hash table to a file
sub writeTbl{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# write array table to a csv file
sub writeArrayTbl
{
  my ($rTbl, $outFile) = @_;

  local $, = ',';
  local $" = ',';
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . ",@{$tbl{$_}}\n"} sort keys %tbl;
  close OUT;
}

# print array table to a csv file
sub printArrayTbl
{
  my ($rTbl, $header) = @_;

  local $, = ',';
  local $" = ',';
  my %tbl = %{$rTbl};
  if ($header)
  {
    print "$header\n";
  }
  map {print $_ . ",@{$tbl{$_}}\n"} sort keys %tbl;
}

# write %children to a file
sub writeChildrenTbl
{
  my $outFile = shift;

  local $, = ',';
  local $" = ',';
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  for my $p (sort keys %children)
  {
    print OUT "$p";
    for ( keys %{$children{$p}} )
    {
      print OUT ",$_";
    }
    print OUT "\n";
  }
  close OUT;
}

# write %spTbl to a csv file
sub writeSpTbl
{
  my $outFile = shift;

  local $, = ',';
  local $" = ',';
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . ",@{$spTbl{$_}}\n"} sort keys %spTbl;
  close OUT;
}


## get random sequence ID from seq IDs of a given species here the selection is using master lineage of all bacteria
## the selected seq should have length greater than in the perc of the max length
sub get_rand_seqID_from_master_sp
{
  my ($sp, $perc) = @_;

  my $c = get_tx_prefix($sp);
  if ($c ne 's')
  {
    print "ERROR in get_rand_seqID_from_sp(): input taxon $sp is not a species\n";
    exit 1;
  }

  if (!defined $masterSpTbl{$sp})
  {
    print "ERROR in get_rand_seqID_from_sp(): Could not find $sp in masterSpTbl\n";
    exit 1;
  }

  my @s = @{$masterSpTbl{$sp}};
  my @len = @masterSeqLen{@s};
  my $imax = argmax(\@len);
  my $maxLen = $len[$imax];

  my $prop = $perc / 100.0;
  my $goodLen = $prop * $maxLen;
  my @ss = grep{ $masterSeqLen{$_} > $goodLen } @s;
  my $id = $ss[rand @ss];

  return $id;
}


## get random sequence ID from seq IDs of a given species
## the selected seq should have lengths in the upper, perc, percentile of all
## sequence of the given species
sub get_rand_seqID_from_sp
{
  my ($sp, $rSpTbl, $rSeqLen, $perc) = @_;

  my $c = get_tx_prefix($sp);
  if ($c ne 's')
  {
    print "ERROR in get_rand_seqID_from_sp(): input taxon $sp is not a species\n";
    exit 1;
  }

  my $prop = $perc / 100.0;
  my %seqLen = %{$rSeqLen};
  my %spTbl = %{$rSpTbl};

  if (!defined $spTbl{$sp})
  {
    print "ERROR in get_rand_seqID_from_sp(): Could not find $sp in spTbl\n";
    exit 1;
  }
  my @s = @{$spTbl{$sp}};
  my @len = @seqLen{@s};
  my $imax = argmax(\@len);
  my $maxLen = $len[$imax];
  my $goodLen = $prop * $maxLen;
  my @ss = grep{ $seqLen{$_} > $goodLen } @s;
  my $id = $ss[rand @ss];

  return $id;
}

## get random sequence from ONE randomly selected species of the given taxon
sub get_one_rand_sp_seq
{
  my ($tx, $rSpTbl, $rSeqLen, $perc) = @_;

  my $prop = $perc / 100.0;
  my %seqLen = %{$rSeqLen};
  my %spTbl = %{$rSpTbl};

  my $c = get_tx_prefix($tx);
  my $id;
  if ($c eq 's')
  {
    $id = get_rand_seqID_from_sp($tx, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'g')
  {
    my @sps = keys %{$children{$tx}};
    # pick randomly a species and then select from it
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'f')
  {
    my @ges = keys %{$children{$tx}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'o')
  {
    my @fas = keys %{$children{$tx}};
    my $fa = $fas[rand @fas];
    my @ges = keys %{$children{$fa}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'c')
  {
    my @ors = keys %{$children{$tx}};
    my $or = $ors[rand @ors];
    my @fas = keys %{$children{$or}};
    my $fa = $fas[rand @fas];
    my @ges = keys %{$children{$fa}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$children{$ge}};
    my $sp = $sps[rand @sps];
    $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
  }
  elsif ($c eq 'p')
  {
    my @cls = keys %{$masterChildren{$tx}};
    my $cl = $cls[rand @cls];
    my @ors = keys %{$masterChildren{$cl}};
    my $or = $ors[rand @ors];
    my @fas = keys %{$masterChildren{$or}};
    my $fa = $fas[rand @fas];
    my @ges = keys %{$masterChildren{$fa}};
    my $ge = $ges[rand @ges];
    my @sps = keys %{$masterChildren{$ge}};
    my $sp = $sps[rand @sps];
    ##$id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
    $id = get_rand_seqID_from_master_sp($sp, $perc);
  }

  return $id;
}

## get random sequence from each species of the given taxon
## the selected seq's should have lengths in the upper, perc, percentile of all
## sequence of the given species
sub get_rand_sp_seq
{
  my ($tx, $rSpTbl, $rSeqLen, $perc) = @_;

  my $prop = $perc / 100.0;
  my %seqLen = %{$rSeqLen};
  my %spTbl = %{$rSpTbl};

  my $c = get_tx_prefix($tx);
  my @seqIDs;

  if ($c eq 's')
  {
    my $id = get_rand_seqID_from_sp($tx, $rSpTbl, $rSeqLen, $perc);
    push @seqIDs, $id;
  }
  elsif ($c eq 'g')
  {
    my @sps = keys %{$children{$tx}};
    for my $sp (@sps)
    {
      my $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
      push @seqIDs, $id;
    }
  }
  elsif ($c eq 'f')
  {
    my @ges = keys %{$children{$tx}};
    for my $ge (@ges)
    {
      my @sps = keys %{$children{$ge}};
      for my $sp (@sps)
      {
	my $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
	push @seqIDs, $id;
      }
    }
  }
  elsif ($c eq 'o')
  {
    my @fas = keys %{$children{$tx}};
    for my $fa (@fas)
    {
      my @ges = keys %{$children{$fa}};
      for my $ge (@ges)
      {
	my @sps = keys %{$children{$ge}};
	for my $sp (@sps)
	{
	  my $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
	  push @seqIDs, $id;
	}
      }
    }
  }
  elsif ($c eq 'c')
  {
    my @ors = keys %{$children{$tx}};
    for my $or (@ors)
    {
      my @fas = keys %{$children{$or}};
      for my $fa (@fas)
      {
	my @ges = keys %{$children{$fa}};
	for my $ge (@ges)
	{
	  my @sps = keys %{$children{$ge}};
	  for my $sp (@sps)
	  {
	    my $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
	    push @seqIDs, $id;
	  }
	}
      }
    }
  }
  elsif ($c eq 'p')
  {
    my @cls = keys %{$children{$tx}};
    for my $cl (@cls)
    {
      my @ors = keys %{$children{$cl}};
      for my $or (@ors)
      {
	my @fas = keys %{$children{$or}};
	for my $fa (@fas)
	{
	  my @ges = keys %{$children{$fa}};
	  for my $ge (@ges)
	  {
	    my @sps = keys %{$children{$ge}};
	    for my $sp (@sps)
	    {
	      my $id = get_rand_seqID_from_sp($sp, $rSpTbl, $rSeqLen, $perc);
	      push @seqIDs, $id;
	    }
	  }
	}
      }
    }
  }

  return @seqIDs;
}


## given a taxon name of the form f_Lactobacillaceae
## extract the taxon prefix, in this example 'f'
sub get_tx_prefix
{
  my $tx = shift;
  my @f = split "_", $tx;
  my $prefix = shift @f;
  return $prefix;
}

sub get_families
{
  my $tx = shift;
  my $c = get_tx_prefix($tx);
  my $i = $txInt{$c};
  if ($i>3)
  {
    print "ERROR: in get_families() at __LINE__ the input taxon $tx is genus or species\n";
    return "";
  }

  my @fas;

  if ($c eq 'f')
  {
    ##print "WARNING: in get_families() the input taxon is a family\n";
    push @fas, $tx;
  }
  elsif ($c eq 'o')
  {
    @fas = keys %{$children{$tx}};
  }
  elsif ($c eq 'c')
  {
    my @ors = keys %{$children{$tx}};
    for (@ors)
    {
      push @fas, keys %{$children{$_}};
    }
  }
  elsif ($c eq 'p')
  {
    my @cls = keys %{$children{$tx}};
    for my $cl (@cls)
    {
      my @ors = keys %{$children{$cl}};
      for my $or (@ors)
      {
	push @fas, keys %{$children{$_}};
      }
    }
  }

  return @fas;
}

sub argmax
{
  my $r = shift;

  my $index = undef;
  my $max   = undef;

  my $count = 0;
  foreach my $val(@{$r})
  {
    if ( not defined $max or $val > $max )
    {
      $max   = $val;
      $index = $count;
    }
    $count++;
  }
  return $index;
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

# difference of two arrays
sub diff
{
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

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle
{
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; )
    {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify {
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}

# get maxKeyLen and maxValLen
sub getKeyValStrLengths{

  my ($rTbl, $rSub) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxKeyLen = 0;
  map { $maxKeyLen = length($_) if( length($_) > $maxKeyLen )} @args;

  my $maxValLen = 0;
  map { $maxValLen = length(commify(scalar(@{$rTbl->{$_}}))) if( length(commify(scalar(@{$rTbl->{$_}}))) > $maxValLen )} @args;

  my @ret = ($maxKeyLen, $maxValLen);
  return \@ret;
}

## print two strings (key and scalar value of a hash table) formated using $maxKeyLen, $maxValLen generated using previeous routine
## (4, $fa, scalar(@{$faTbl{$fa}}), $faKVL);
sub printF
{
  my ($nTabs, $key, $val, $rKV ) = @_;

  my $tabs = "";
  my $i = 0;
  while ( $i < $nTabs)
  {
    $tabs .= "\t";
    $i++;
  }

  my $maxKeyLen = $rKV->[0];
  my $maxValLen = $rKV->[1];
  my $n = ($maxKeyLen - length($key)) + ($maxValLen - length($val));
  my $pad = ": ";
  for (my $i=0; $i<$n; $i++)
  {
    $pad .= " ";
  }
  print "$tabs$key$pad$val\n";
}


# print elements of a hash table so that arguments are aligned
sub printFormatedTbl{

  my ($rTbl, $rSub) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxKeyLen = 0;
  map { $maxKeyLen = length($_) if( length($_) > $maxKeyLen )} @args;

  my $maxValLen = 0;
  map { $maxValLen = length(commify(scalar(@{$rTbl->{$_}}))) if( length(commify(scalar(@{$rTbl->{$_}}))) > $maxValLen )} @args;

  for (@args)
  {
    ##my $n = $maxKeyLen - length($_);
    my $n = ($maxKeyLen - length($_)) + ($maxValLen - length(commify(scalar(@{$rTbl->{$_}}))));
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "\t$_$pad" . commify(scalar(@{$rTbl->{$_}})) . "\n";
  }
}

# print elements of a hash table
sub printTbl
{
  my ($rTbl, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR: $file does not exist\n\n\n";
    exit 1;
  }

  my %idTbl;
  #my %lineageTbl;

  open IN, "<$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $lineage) = split /\s+/,$_;
    $idTbl{$id} = $lineage;
    #push @{$lineageTbl{$lineage}}, $id;
  }
  close IN;

  return %idTbl;
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

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
  print "\n";
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# read table with one column
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    print "\n\nERROR in readArray() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit 1;
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

exit 0;
