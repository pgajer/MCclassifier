#!/usr/bin/env perl

=head1 NAME

  outgroup_rectifier.pl


=head1 DESCRIPTION

  taxonomic_reduction.pl produces groups with outgroups (OG) from random
  sequences of siblings of the lowest rank _singleton taxon_ of the group. For
  example,

  Group 0
  p_Actinobacteria
    c_Actinobacteria
    c_Thermoleophilia

  The singleton taxon is p_Actinobacteria so the OG will be from other phyla.

  or

  Group 1
  c_Actinobacteria
	o_Actinomycetales
		f_Actinomycetaceae
		f_Brevibacteriaceae
		f_Dermabacteraceae
		f_Dietziaceae
		f_Intrasporangiaceae
		f_Mycobacteriaceae
		f_Nocardiopsaceae
		f_Promicromonosporaceae
		f_Streptosporangiaceae
		f_Thermomonosporaceae

  Here the lowest rank singleton taxon is o_Actinomycetales so the OG of this
  group is from other orders of the c_Actinobacteria class.

  Strangely enought, sometimes OG sequences do not form a homogenous outgroup
  clade and form two or more separate clusters. The aim of this script is to
  identify these clusters and take the larger one as the new outgoup. Once this
  is done and tested, the alignment, lineage, taxon files are updated together
  with trees. Removed OG seqIDs are written to a file.


=head1 SYNOPSIS

  outgroup_rectifier.pl -i <group prefix> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

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

  cd ~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Actinobacteria_dir

  outgroup_rectifier.pl --keep-tmp-files --debug -i Actinobacteria_group_2

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

GetOptions(
  "group-prefix|i=s" => \my $grPrefix,
  "keep-tmp-files"   => \my $keepTmpFiles,
  "report-only"      => \my $reportOnly,
  "quiet"            => \my $quiet,
  "verbatim|v"       => \my $verbatim,
  "debug"            => \my $debug,
  "dry-run"          => \my $dryRun,
  "help|h!"          => \my $help,
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

####################################################################
##                               MAIN
####################################################################

# For a given group prefix Actinobacteria_group_0 the following files/directories
# has been created by taxonomic_reduction.pl and taxonomy_cleanup.pl scripts:

# Actinobacteria_group_0.fa
# Actinobacteria_group_0.lineage
# Actinobacteria_group_0.seqIDs
# Actinobacteria_group_0.tx
# Actinobacteria_group_0.txSummary
# Actinobacteria_group_0_algn.fa
# Actinobacteria_group_0_outgroup.seqIDs
# Actinobacteria_group_0_dir

# this script will generate deleted outgroup seqIDs file

# Actinobacteria_group_0_outgroup_deleted.seqIDs

# and will update the above files.

my $grDir     = $grPrefix . "_dir";
$grPrefix     = "$grDir/$grPrefix";
my $errorFile = $grPrefix . "_outgroup_rectifier_ERROR"; # this is a file that indicates that outgroup_rectifier.pl finished with an error

print "\n\nerrorFile: $errorFile\n" if $debug;

if ( -e $errorFile )
{
  unlink($errorFile);
}

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  touchFile( $errorFile );
  exit;
}

my $faGrFile        = $grPrefix . ".fa";
my $lineageFile     = $grPrefix . ".lineage";
my $txGrFile        = $grPrefix . ".tx";
my $txSummaryGrFile = $grPrefix . ".txSummary";
my $faAlgnGrFile    = $grPrefix . "_algn.fa";
my $trimmedAlgnFile = $grPrefix . "_algn_trimmed.fa";
my $ogSeqIDsFile    = $grPrefix . "_outgroup.seqIDs";
my $delogSeqIDsFile = $grPrefix . "_outgroup_deleted.seqIDs";
my $seqIDsFile      = $grPrefix . "_after_pruning.seqIDs";


print "--- Parsing lineage table\n" if !$quiet;
my %lineageTbl = read2colTbl($lineageFile);

print "--- Reading OG sequence IDs\n" if !$quiet;
my @ogSeqIDs = readArray($ogSeqIDsFile);

printArrayByRow(\@ogSeqIDs, "ogSeqIDs") if $debug;

if ( @ogSeqIDs == 1 )
{
  print "\n\tNumber of outgroup seq's: " . @ogSeqIDs . "\n";
  print "\tNothing to be done\n\n";
  exit;
}
elsif ( @ogSeqIDs == 0 )
{
  warn "\n\n\tERROR: No outgroup sequences found";
  print "\n\n\n";
  touchFile( $errorFile );
  exit;
}

my @allSeqIDs = keys %lineageTbl;
my @noOGseqIDs = diff(\@allSeqIDs, \@ogSeqIDs);

print "\n\tNumber of all seq's: " . @allSeqIDs . "\n" if !$quiet;
print "\tNumber of OG seq's: " . @ogSeqIDs . "\n" if !$quiet;
print "\tNumber of non-OG seq's: " . @noOGseqIDs . "\n\n" if !$quiet;


print "--- Parsing phylo spp_seqIDs tree to test if OG seq's form one or more clusters\n" if !$quiet;

my $sppSeqIdTreeFile = $grPrefix . "_sppSeqIDs.tree";
if ( ! -e $sppSeqIdTreeFile )
{
  warn "ERROR: $sppSeqIdTreeFile does not exist";
  touchFile( $errorFile );
  exit;
}

print "--- Extracting leaves\n" if !$quiet;
my $sppSeqIdTreeLeavesFile = "$grPrefix" . "_sppSeqIDs.leaves";
my $cmd = "rm -f $sppSeqIdTreeLeavesFile; nw_labels -I $sppSeqIdTreeFile > $sppSeqIdTreeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Reading leaves\n" if !$quiet;
my @leaves = readArray($sppSeqIdTreeLeavesFile);

print "--- Checking the number of clusters formed by OG seqs\n" if !$quiet;
my @ogIdx;
for my $i (0..$#leaves)
{
  push @ogIdx, $i if ( $leaves[$i] =~ /_OG_/ );
}

printArray(\@ogIdx, "\nPositions of OG seq's") if ($debug);

## identifying consecutive indices ranges
my @start;
my @end;

push @start, $ogIdx[0];
if (@ogIdx>1)
{
  for my $i (1..$#ogIdx)
  {
    ##if ( !$foundEnd && $ogIdx[$i] != $start[$rIdx] )
    if ( $ogIdx[$i-1]+1 != $ogIdx[$i] )
    {
      #$foundEnd = 1;
      push @end, $ogIdx[$i-1];
      push @start, $ogIdx[$i];
    }
    if ($i==$#ogIdx)
    {
      push @end, $ogIdx[$i];
    }

    if (0 && $debug)
    {
      print "\ni: $i\n";
      printArray(\@start, "start");
      printArray(\@end, "end");
    }
  }
}
else
{
  push @end, $ogIdx[0];
}

my @ogPos1;  # OG positions
for my $i (0..$#start)
{
  push @ogPos1, ($start[$i] .. $end[$i]);
}

my @og = @leaves[@ogPos1];
my @ogBig = @leaves[($start[0] .. $end[$#end])];

printArrayByRow(\@og, "\nOutgroup elements") if ($debug);

if ( scalar(@start) != scalar(@end) )
{
  warn "$grPrefix\nERROR: start and end arrays have different lengths!";
  print "length(start): " . @start . "\n";
  print "length(end): " . @end . "\n";
  touchFile( $errorFile );
  exit;
}

my @rangeSize;
for my $i (0..$#start)
{
  push @rangeSize, ($end[$i] - $start[$i]+1);
}

if ($debug)
{
  print "\nstart\tend\tsize\n";
  for my $i (0..$#start)
  {
    print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
  }
  print "\n";
}

if ($reportOnly)
{
  print "\nstart\tend\tsize\n";
  for my $i (0..$#start)
  {
    print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
  }
  print "\n";

  if (@rangeSize>1)
  {
    print "Detected multiple OG clusters\n";

    my $imax = argmax(\@rangeSize);
    print "imax: $imax\n";
    print "Maximal range size: " . $rangeSize[$imax] . "\n";

    my $minCladeSize = @leaves;
    my $minCladeSizeIdx = $imax;
    print "Clade size of each cluster of maximal range size\n";
    print "\nidx\tstart\tend\trgSize\tcladeSize\n";
    for my $i (0..$#rangeSize)
    {
      #if ($rangeSize[$i] == $rangeSize[$imax])
      if (1)
      {
	my @pos = ($start[$i] .. $end[$i]);
	my @og = @leaves[@pos];

	#print "--- Extracting the clade of OG sequences\n";
	my $ogCladeTreeFile = "$grPrefix" . "_clade.tree";
	$cmd = "rm -f $ogCladeTreeFile; nw_clade $sppSeqIdTreeFile @og > $ogCladeTreeFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

	#print "--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$grPrefix" . "_clade.leaves";
	$cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

	#print "--- Reading the leaves\n" if !$quiet;
	my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

	print "i: $i\t$start[$i]\t$end[$i]\t$rangeSize[$i]\t" . @ogCladeLeaves . "\n";
	if ( @ogCladeLeaves < $minCladeSize )
	{
	  $minCladeSize = @ogCladeLeaves;
	  $minCladeSizeIdx = $i;
	}
      }
    }

    $imax = $minCladeSizeIdx;
    print "\nUpdated imax: $imax\n";
  }
  print "\n";
  exit;
}

## if there is more than one cluster, eliminate all but the largest one
my @badOGs;
my $nCltrs = @rangeSize;
if (@rangeSize>1)
{
  ## from http://stackoverflow.com/questions/18532026/how-to-append-system-date-to-a-filename-in-perl
  my @now = localtime();
  my $timeStamp = sprintf("%04d-%02d-%02d_%02d:%02d:%02d",
			  $now[5]+1900, $now[4]+1, $now[3],
			  $now[2],      $now[1],   $now[0]);

  my $origDataDir = $grDir . "/orig_data_dir_$timeStamp";

  $cmd = "mkdir -p $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug; # || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


  print "Detected multiple OG clusters\n" if $debug;

  ## selecting imax with the smaller clade size
  my $imax = argmax(\@rangeSize);
  print "imax: $imax\n" if $debug;
  print "Maximal range size: " . $rangeSize[$imax] . "\n" if $debug;

  my $minCladeSize = @leaves;
  my $minCladeSizeIdx = $imax;
  print "Clade size of each cluster of maximal range size\n" if $debug;
  print "\nidx\tstart\tend\trgSize\tcladeSize\n" if $debug;
  for my $i (0..$#rangeSize)
  {
    if ($rangeSize[$i] == $rangeSize[$imax])
    {
      my @pos = ($start[$i] .. $end[$i]);
      my @og = @leaves[@pos];

      #print "--- Extracting the clade of OG sequences\n";
      my $ogCladeTreeFile = "$grPrefix" . "_clade.tree";
      $cmd = "rm -f $ogCladeTreeFile; nw_clade $sppSeqIdTreeFile @og > $ogCladeTreeFile";
      #print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

      #print "--- Extracting leaves of the OG clade\n";
      my $ogCladeTreeLeavesFile = "$grPrefix" . "_clade.leaves";
      $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
      #print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

      #print "--- Reading the leaves\n" if !$quiet;
      my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);
      print "i: $i\t$start[$i]\t$end[$i]\t$rangeSize[$i]\t" . @ogCladeLeaves . "\n" if $debug;
      if ( @ogCladeLeaves < $minCladeSize )
      {
	$minCladeSize = @ogCladeLeaves;
	$minCladeSizeIdx = $i;
      }
    }
  }

  $imax = $minCladeSizeIdx;
  print "\nUpdated imax: $imax\n" if $debug;

  my @delPos; # row numbers (0-based) to be deleted
  my @ogPos;  # OG positions (after removing the small clusters)
  for my $i (0..$#start)
  {
    if ($i != $imax)
    {
      push @delPos, ($start[$i] .. $end[$i]);
    }
    else
    {
      push @ogPos, ($start[$i] .. $end[$i]);
    }
  }

  printArray(\@delPos, "\ndelPos") if ($debug);
  print "\n" if ($debug);

  printArray(\@ogPos, "\nogPos") if ($debug);
  print "\n" if ($debug);

  # Corynebacterium_kroppenstedtii__S000266890
  # Corynebacterium_kroppenstedtii__S001416212
  # Metascardovia_criceti_OG__S000616640
  # Sciscionella_marina__S001016157
  # Thermocrispum_agreste__S002161645
  # Thermocrispum_agreste__S002161657

  #my $dx = 2;
  @badOGs = @leaves[@delPos];
  printArrayByRow(\@badOGs, "badOGs") if ($debug);

  ## new outgroup elements
  @og = @leaves[@ogPos];

  printArrayByRow(\@og, "og") if ($debug);

  # Extracting seqIDs from new outgroup elements
  my @ogIDs;
  for (@og)
  {
    my ($s, $id) = split '__', $_;
    push @ogIDs, $id;
    #print "$_\tID: $id\n";
  }
  #print "\n";

  ## Identifying deleted OG seqIDs and writing them to file
  my @delogIDs;
  for (@badOGs)
  {
    my ($s, $id) = split '__', $_;
    push @delogIDs, $id;
    #print "$_\tID: $id\n";
  }
  #print "\n";

  printArrayByRow(\@delogIDs, "deleted OG seqIDs") if ($debug);
  printArrayByRow(\@ogIDs, "new OG seqIDs") if ($debug);

  writeArray(\@delogIDs, $delogSeqIDsFile);

  ## non-outgroup seqIDs with the OG ones that survived pruning
  my @seqIDs = (@noOGseqIDs, @ogIDs);

  writeArray(\@seqIDs, $seqIDsFile);

  print "--- Pruning lost seqIDs from the current sppSeqIDs phylo tree\n" if !$quiet;
  my $sppSeqIdTreeFile = "$grPrefix" . "_sppSeqIDs.tree";
  my $prunedTreeFile = "$grPrefix" . "_sppSeqIDs_pruned.tree";
  $cmd = "rm -f $prunedTreeFile; nw_prune $sppSeqIdTreeFile @badOGs > $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Rerooting the tree using pruned outgroup seq's\n" if !$quiet;
  my $prunedRerootedTreeFile = "$grPrefix" . "_sppSeqIDs_pruned_rr.tree";
  $cmd = "rm -f $prunedRerootedTreeFile; nw_reroot $prunedTreeFile @og | nw_order -  > $prunedRerootedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ##print "\n\tPruned rerooted tree written to $prunedRerootedTreeFile\n\n" if $debug;

  ## save the original tree in $origDataDir
  $cmd = "mv $sppSeqIdTreeFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  $cmd = "mv $prunedRerootedTreeFile $sppSeqIdTreeFile; rm -f $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;


  my $urTreeFile    = "$grPrefix" . "_unrooted.tree";
  my $tmpurTreeFile = "$grPrefix" . "_unrooted_tmp.tree";

  my $treeFile    = "$grPrefix" . ".tree";
  my $tmpTreeFile = "$grPrefix" . "_tmp.tree";

  print "--- Pruning lost seqIDs from the unrooted phylo tree\n" if !$quiet;
  $cmd = "rm -f $tmpurTreeFile; nw_prune $urTreeFile @delogIDs > $tmpurTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Rerooting the tree using pruned outgroup seq's\n" if !$quiet;
  $cmd = "rm -f $tmpTreeFile; nw_reroot $tmpurTreeFile @ogIDs | nw_order -  > $tmpTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ## save the original $urTreeFile in $origDataDir
  $cmd = "mv $urTreeFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  $cmd = "mv $tmpurTreeFile $urTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  ## save the original $treeFile in $origDataDir
  $cmd = "mv $treeFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  $cmd = "mv $tmpTreeFile $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  print "--- Pruning remaining files\n" if !$quiet;

  # print "--- Extracting leaves of the OG clade\n";
  # my $treeLeavesFile = "$grPrefix" . ".leaves";
  # $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # #print "--- Reading the leaves\n" if !$quiet;
  # my @prLeaves = readArray($treeLeavesFile);

  # my $faGrFile        = $grPrefix . ".fa";
  # my $lineageFile     = $grPrefix . ".lineage";
  # my $txGrFile        = $grPrefix . ".tx";
  # my $faAlgnGrFile    = $grPrefix . "_algn.fa";
  # my $ogSeqIDsFile    = $grPrefix . "_outgroup.seqIDs";

  my $tmpfaGrFile        = $grPrefix . "_tmp.fa";
  my $tmplineageFile     = $grPrefix . "_tmp.lineage";
  my $tmptxGrFile        = $grPrefix . "_tmp.tx";
  my $tmpfaAlgnGrFile    = $grPrefix . "_tmp_algn.fa";
  #my $tmpogSeqIDsFile    = $grPrefix . "_tmp_outgroup.seqIDs";
  my $tmpTrimmedAlgnFile = $grPrefix . "_tmp_algn_trimmed.fa";

  ## fa
  print "--- Pruning $faGrFile\n" if !$quiet;
  $cmd = "select_seqs.pl --quiet -s $seqIDsFile -i $faGrFile -o $tmpfaGrFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ## save the original fa in $origDataDir
  $cmd = "mv $faGrFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  $cmd = "mv $tmpfaGrFile $faGrFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;


  ## algn
  if ( -e $faAlgnGrFile)
  {
    print "--- Pruning $faAlgnGrFile\n" if !$quiet;
    $cmd = "select_seqs.pl --quiet -s $seqIDsFile -i $faAlgnGrFile -o $tmpfaAlgnGrFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    ## save the original algn in $origDataDir
    $cmd = "mv $faAlgnGrFile $origDataDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

    $cmd = "mv $tmpfaAlgnGrFile $faAlgnGrFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;
  }

  ## trimmed algn
  if ( -e $trimmedAlgnFile)
  {
    print "--- Pruning $trimmedAlgnFile\n" if !$quiet;
    $cmd = "select_seqs.pl --quiet -s $seqIDsFile -i $trimmedAlgnFile -o $tmpTrimmedAlgnFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    ## save the original trimmed algn in $origDataDir
    $cmd = "mv $trimmedAlgnFile $origDataDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

    $cmd = "mv $tmpTrimmedAlgnFile $trimmedAlgnFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;
  }

  ## tx
  print "--- Pruning $txGrFile\n" if !$quiet;
  $cmd = "select_tx.pl --quiet -s $seqIDsFile -i $txGrFile -o $tmptxGrFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ## save the original tx file in $origDataDir
  $cmd = "mv $txGrFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  $cmd = "mv $tmptxGrFile $txGrFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;


  ## lineage
  print "--- Pruning $lineageFile\n" if !$quiet;
  $cmd = "select_tx.pl --quiet -s $seqIDsFile -i $lineageFile -o $tmplineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ## save the original lineage file in $origDataDir
  $cmd = "mv $lineageFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  $cmd = "mv $tmplineageFile $lineageFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  ## OG
  ## save the original OG file in $origDataDir
  $cmd = "mv $ogSeqIDsFile $origDataDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  writeArray(\@ogIDs, $ogSeqIDsFile);


  print "--- Testing if the resulting tree has OGs in one cluster at the top or bottom of the tree\n" if !$quiet;

  print "--- Extracting leaves\n" if !$quiet;
  # my $sppSeqIdTreeLeavesFileRR = "$grPrefix" . "_sppSeqIDs_pruned_rr.leaves";
  # $cmd = "rm -f $sppSeqIdTreeLeavesFileRR; nw_labels -I $prunedRerootedTreeFile > $sppSeqIdTreeLeavesFileRR";
  my $sppSeqIdTreeLeavesFile = "$grPrefix" . "_sppSeqIDs.leaves";
  $cmd = "rm -f $sppSeqIdTreeLeavesFile; nw_labels -I $sppSeqIdTreeFile > $sppSeqIdTreeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Reading leaves\n" if !$quiet;
  @leaves = readArray($sppSeqIdTreeLeavesFile);

  $cmd = "rm -f $sppSeqIdTreeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun && !$keepTmpFiles;

  print "--- Checking the number of clusters formed by pruned OG seqs\n" if !$quiet;
  @ogIdx = ();
  for my $i (0..$#leaves)
  {
    push @ogIdx, $i if ( $leaves[$i] =~ /_OG_/ );
  }

  printArray(\@ogIdx, "\nogIdx") if $debug;

  ## identifying consecutive indices ranges
  @start = ();
  @end   = ();
  push @start, $ogIdx[0];
  if (@ogIdx>1)
  {
    for my $i (1..$#ogIdx)
    {
      if ( $ogIdx[$i-1]+1 != $ogIdx[$i] )
      {
	push @end, $ogIdx[$i-1];
	push @start, $ogIdx[$i];
      }
      if ($i==$#ogIdx)
      {
	push @end, $ogIdx[$i];
      }

      if (0 && $debug)
      {
	print "\ni: $i\n";
	printArray(\@start, "start");
	printArray(\@end, "end");
      }
    }
  }
  else
  {
    push @end, $ogIdx[0];
  }

  if ( scalar(@start) != scalar(@end) )
  {
    warn "$grPrefix\nERROR: start and end arrays have different lengths!";
    print "length(start): " . @start . "\n";
    print "length(end):   " . @end . "\n";
    touchFile( $errorFile );
    exit;
  }

  @rangeSize = ();
  for my $i (0..$#start)
  {
    push @rangeSize, ($end[$i] - $start[$i]+1);
  }

  if ($debug)
  {
    print "Number of leaves: " . @leaves . "\n";
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";
  }

  if (@start > 1)
  {
    warn "$grPrefix\nERROR: In the pruned tree outgroups still are found in more than one cluster!";
    touchFile( $errorFile );

    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";

    printArrayByRow(\@og, "og");
    print "\n";

    my $maxOGbigSize = 100;
    if ( @ogBig < $maxOGbigSize )
    {
      printArrayByRow(\@ogBig, "Leaves from first to last OG seq");
    }

    #print "--- Extracting the clade of OG sequences\n";
    my $ogCladeTreeFile = "$grPrefix" . "_sppSeqIDs_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; nw_clade $sppSeqIdTreeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_sppSeqIDs_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "--- Reading the leaves\n" if !$quiet;
    my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

    my $maxCladeSize = 100;
    if ( @ogCladeLeaves < $maxCladeSize )
    {
      printArrayByRow(\@ogCladeLeaves, "OG Clade Leaves");
    }

    print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
    print   "\tNumber of OG sequences: " . @og . "\n\n";

    exit;
  }

  if ( !( $start[0] == 0 || $end[0] == $#leaves) )
  {
    warn "$grPrefix\nERROR: In the pruned tree outgroups sequences are not at the top or bottom of the tree!";
    touchFile( $errorFile );

    print "Number of leaves: " . @leaves . "\n";
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";

    printArrayByRow(\@og, "og");
    print "\n";

    my $maxOGbigSize = 100;
    if ( @ogBig < $maxOGbigSize )
    {
      printArrayByRow(\@ogBig, "Leaves from first to last OG seq");
    }

    #print "--- Extracting the clade of OG sequences\n";
    my $ogCladeTreeFile = "$grPrefix" . "_sppSeqIDs_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; nw_clade $sppSeqIdTreeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_sppSeqIDs_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "--- Reading the leaves\n" if !$quiet;
    my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

    my $maxCladeSize = 100;
    if ( @ogCladeLeaves < $maxCladeSize )
    {
      printArrayByRow(\@ogCladeLeaves, "OG Clade Leaves");
    }

    print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
    print   "\tNumber of OG sequences: " . @og . "\n\n";

    exit;
  }
}
else
{
  print "--- Testing if the outgropu seq's are at the top or bottom of the tree\n" if !$quiet;
  printArray(\@ogIdx, "\nogIdx") if $debug;

  ## identifying consecutive indices ranges
  my @start = ();
  my @end   = ();
  push @start, $ogIdx[0];
  if (@ogIdx>1)
  {
    for my $i (1..$#ogIdx)
    {
      if ( $ogIdx[$i-1]+1 != $ogIdx[$i] )
      {
	push @end, $ogIdx[$i-1];
	push @start, $ogIdx[$i];
      }
      if ($i==$#ogIdx)
      {
	push @end, $ogIdx[$i];
      }
    }
  }
  else
  {
    push @end, $ogIdx[0];
  }

  if ( scalar(@start) != scalar(@end) )
  {
    warn "$grPrefix\nERROR: start and end arrays have different lengths!";
    print "length(start): " . @start . "\n";
    print "length(end): " . @end . "\n";
    touchFile( $errorFile );
    exit;
  }

  @rangeSize = ();
  for my $i (0..$#start)
  {
    push @rangeSize, ($end[$i] - $start[$i]+1);
  }

  if ($debug)
  {
    print "Number of leaves: " . @leaves . "\n";
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";
  }

  if ( !( $start[0] == 0 || $end[0] == $#leaves) )
  {
    warn "$grPrefix\nERROR: Outgroups sequences are not at the top or bottom of the tree!";
    touchFile( $errorFile );

    print "Number of leaves: " . @leaves . "\n";
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n\n";

    printArrayByRow(\@og, "og");

    #print "--- Extracting the clade of OG sequences\n";
    my $ogCladeTreeFile = "$grPrefix" . "_sppSeqIDs_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; nw_clade $sppSeqIdTreeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_sppSeqIDs_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "--- Reading the leaves\n" if !$quiet;
    my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

    my $maxCladeSize = 100;
    if (@ogCladeLeaves<$maxCladeSize)
    {
      printArrayByRow(\@ogCladeLeaves, "OG Clade Leaves");
    }

    print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
    print   "\tNumber of OG sequences: " . @og . "\n\n";

    exit;
  }

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

print "\n\tDeleted outgroup seqIDs written to $delogSeqIDsFile\n" if !$quiet;

print "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n";
print "\tNumber of outgroup clusters: " . $nCltrs . "\n";
if ($nCltrs>1)
{
  print "\tNumber of outgroup seq's removed: " . @badOGs . "\n";
  print "\tNumber of remined outgroup seq's: " . @og . "\n\n";
}
else
{
  print "\tNo outgroup seq's has been removed\n\n";

}

####################################################################
##                               SUBS
####################################################################

## get the lowest rank singleton taxon
## the output is the single character
## 's','g','f','o','c','p'
## of the lowest rank singleton taxon
sub get_lowest_rank_singleton_taxon
{
  my $rTaxons = shift;

  my $debug = 1;

  my @taxons = @{$rTaxons};

  my %txInt;
  $txInt{'p'} = 0;
  $txInt{'c'} = 1;
  $txInt{'o'} = 2;
  $txInt{'f'} = 3;
  $txInt{'g'} = 4;
  $txInt{'s'} = 5;

  my %grTx; # group taxons
  for my $tx ( @taxons )
  {
    my $c = get_tx_prefix($tx);
    my $i = $txInt{$c};
    push @{$grTx{$c}}, $tx;
  }

  printArrayTbl(\%grTx, "grTx") if $debug;

  ## pick the lowest taxon with the single element

  my $imin = 5;
  my $cmin;
  my @cc = ('s','g','f','o','c','p');
  for my $c (@cc)
  {
    if ( exists $grTx{$c} && @{$grTx{$c}}==1 )
    {
      $cmin = $c;
      last;
    }
  }

  return $cmin;
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

# read table with one column
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file && ! -l $file )
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

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;

  local $" = ', ';
  ##local $, = ',';
  print "$header: " if $header;
  print "@{$a}\n";
}

# print array to stdout
sub printArrayByRow
{
  my ($a, $header) = @_;

  local $" = '\n ';
  ##local $, = ',';
  print "$header:\n" if $header;
  map { print "$_\n" } @{$a};
  print "\n";
}

sub printLineage
{
  my ($r, $p) = @_;

  my $showSpp = 0;

  my %lineageTbl = %{$r};

  my %spTbl;
  my %geTbl;
  my %faTbl;
  my %orTbl;
  my %clTbl;
  my %phTbl;

  my %children;

  for my $id ( keys %lineageTbl )
  {
    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    ##$sp = "s_$sp";
    ##$sp .= "_OG" if ( exists $ogInd{$id} );
    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $children{$ph}{$cl}++;
    $children{$cl}{$or}++;
    $children{$or}{$fa}++;
    $children{$fa}{$ge}++;
    $children{$ge}{$sp}++;

    $spTbl{$sp}++;
    $geTbl{$ge}++;
    $faTbl{$fa}++;
    $orTbl{$or}++;
    $clTbl{$cl}++;
    $phTbl{$ph}++;
  }


  if (0 && !defined $p)
  {
    print "\nPhyla:\n";
    printTbl(\%phTbl);

    print "Classes:\n";
    printTbl(\%clTbl);

    print "Orders:\n";
    printTbl(\%orTbl);

    print "Families:\n";
    printTbl(\%faTbl);
    print "\n";
  }


  my @phs;
  if (defined $p)
  {
    @phs = @{$p};
  }
  else
  {
   @phs = keys %phTbl;
  }
  my $phKVL = getKeyValStrLengths(\%phTbl);
  for my $ph (@phs)
  {
    printF(0, $ph, $phTbl{$ph}, $phKVL);
    my @cls = keys %{$children{$ph}};
    my $clKVL = getKeyValStrLengths(\%clTbl, \@cls);
    for my $cl ( sort{$clTbl{$a} <=> $clTbl{$b}} @cls)
    {
      printF(1, $cl, $clTbl{$cl}, $clKVL);
      my @ors = keys %{$children{$cl}};
      my $orKVL = getKeyValStrLengths(\%orTbl, \@ors);
      for my $or ( sort{$orTbl{$a} <=> $orTbl{$b}} @ors)
      {
	printF(2, $or, $orTbl{$or}, $orKVL);
	my @fas = keys %{$children{$or}};
	my $faKVL = getKeyValStrLengths(\%faTbl, \@fas);
	for my $fa ( sort{$faTbl{$a} <=> $faTbl{$b}} @fas)
	{
	  printF(3, $fa, $faTbl{$fa}, $faKVL);
	  if ($showSpp)
	  {
	    my @ges = keys %{$children{$fa}};
	    my $geKVL = getKeyValStrLengths(\%geTbl, \@ges);
	    for my $ge ( sort{$geTbl{$a} <=> $geTbl{$b}} @ges)
	    {
	      printF(4, $ge, $geTbl{$ge}, $geKVL);
	      my @sps = keys %{$children{$ge}};
	      my $spKVL = getKeyValStrLengths(\%spTbl, \@sps);
	      for my $sp ( sort{$spTbl{$a} <=> $spTbl{$b}} @sps)
	      {
	        printF(5, $sp, $spTbl{$sp}, $spKVL);
	      }
	    }
	  }
	}
      }
    }
  }

  my @fa = sort keys %faTbl;
  return (\@phs, \@fa);
}

sub printLineage2
{
  my ($r, $showSpp) = @_;

  my %lineageTbl = %{$r};

  my %spTbl;
  my %geTbl;
  my %faTbl;
  my %orTbl;
  my %clTbl;
  my %phTbl;

  my %children;

  for my $id ( keys %lineageTbl )
  {
    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    ##$sp = "s_$sp";
    ##$sp .= "_OG" if ( exists $ogInd{$id} );
    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $children{$ph}{$cl}++;
    $children{$cl}{$or}++;
    $children{$or}{$fa}++;
    $children{$fa}{$ge}++;
    $children{$ge}{$sp}++;

    $spTbl{$sp}++;
    $geTbl{$ge}++;
    $faTbl{$fa}++;
    $orTbl{$or}++;
    $clTbl{$cl}++;
    $phTbl{$ph}++;
  }

  my @phs = keys %phTbl;
  my $phKVL = getKeyValStrLengths(\%phTbl);
  for my $ph (@phs)
  {
    printF(0, $ph, $phTbl{$ph}, $phKVL);
    my @cls = keys %{$children{$ph}};
    my $clKVL = getKeyValStrLengths(\%clTbl, \@cls);
    for my $cl ( sort{$clTbl{$a} <=> $clTbl{$b}} @cls)
    {
      printF(1, $cl, $clTbl{$cl}, $clKVL);
      my @ors = keys %{$children{$cl}};
      my $orKVL = getKeyValStrLengths(\%orTbl, \@ors);
      for my $or ( sort{$orTbl{$a} <=> $orTbl{$b}} @ors)
      {
	printF(2, $or, $orTbl{$or}, $orKVL);
	my @fas = keys %{$children{$or}};
	my $faKVL = getKeyValStrLengths(\%faTbl, \@fas);
	for my $fa ( sort{$faTbl{$a} <=> $faTbl{$b}} @fas)
	{
	  printF(3, $fa, $faTbl{$fa}, $faKVL);
	  if ($showSpp)
	  {
	    my @ges = keys %{$children{$fa}};
	    my $geKVL = getKeyValStrLengths(\%geTbl, \@ges);
	    for my $ge ( sort{$geTbl{$a} <=> $geTbl{$b}} @ges)
	    {
	      printF(4, $ge, $geTbl{$ge}, $geKVL);
	      my @sps = keys %{$children{$ge}};
	      my $spKVL = getKeyValStrLengths(\%spTbl, \@sps);
	      for my $sp ( sort{$spTbl{$a} <=> $spTbl{$b}} @sps)
	      {
	        printF(5, $sp, $spTbl{$sp}, $spKVL);
	      }
	    }
	  }
	}
      }
    }
  }

  my @fa = sort keys %faTbl;
  return (\@phs, \@fa);
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
  map { $maxValLen = length(commify($rTbl->{$_})) if( length(commify($rTbl->{$_})) > $maxValLen )} @args;

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

sub get_taxons
{
  my $rLineage = shift;

  my %lineageTbl = %{$rLineage};

  my @taxons;

  my %spTbl;
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
    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    ##$sp = "s_$sp";
    ##$sp .= "_OG" if ( exists $ogInd{$id} );
    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $parent{$sp} = $ge;
    $parent{$ge} = $fa;
    $parent{$fa} = $or;
    $parent{$or} = $cl;
    $parent{$cl} = $ph;
    $parent{$ph} = "d_Bacteria";

    $children{"d_Bacteria"}{$ph}++;
    $children{$ph}{$cl}++;
    $children{$cl}{$or}++;
    $children{$or}{$fa}++;
    $children{$fa}{$ge}++;
    $children{$ge}{$sp}++;

    $spTbl{$sp}++;
    $geTbl{$ge}++;
    $faTbl{$fa}++;
    $orTbl{$or}++;
    $clTbl{$cl}++;
    $phTbl{$ph}++;
  }


  my @phs = keys %phTbl;
  for my $ph (@phs)
  {
    push @taxons, $ph;
    my @cls = keys %{$children{$ph}};
    for my $cl (@cls)
    {
      push @taxons, $cl;
      my @ors = keys %{$children{$cl}};
      for my $or (@ors)
      {
	push @taxons, $or;
	my @fas = keys %{$children{$or}};
	push @taxons, @fas;
      }
    }
  }

  return (\@taxons, \%parent, \%children);
}

# use a two-ranks-up rule for OG selection: find the highest taxonomic rank (say
# class), go to grandparent (the parent of the parent) (Bacteria in the case of
# the class), take for OG children of the grandparent and subtract from them our
# phylum. Those are taxonomic ranks from which the OGs will be selected. In the
# case of family, we start from genus and go two levels up, so we are in the
# order and so all other families of that order would be the OGs of the given
# family

# For Group 0 of Firmicutes - Find the highest taxonomic rank within Group 0 - it
# is class. Go to the grandparent - Bacteria. Take for OGs the children of the
# grandparent (phyla) that are not the phylum of Group 0.

# For Bacilli, which is separated at lower ranks, say family, outgroups would be
# chosen from sister orders excluding the order containing the family of interest.

sub get_OG_families
{
  my ($rTaxons, $rParent, $rChildren) = @_;

  my @taxons = @{$rTaxons};
  my %parent = %{$rParent};
  my %children = %{$rChildren};

  my %txInt;
  $txInt{'p'} = 0;
  $txInt{'c'} = 1;
  $txInt{'o'} = 2;
  $txInt{'f'} = 3;
  $txInt{'g'} = 4;
  $txInt{'s'} = 5;

  my %grTx; # group taxons
  my %pTbl;
  my @d;
  my @ch;
  my $imin = 5;
  my $cmin;

  my $debug = 1;

  for my $tx ( @taxons )
  {
    my $c = get_tx_prefix($tx);
    my $i = $txInt{$c};
    push @{$grTx{$c}}, $tx;
  }

  printArrayTbl(\%grTx, "grTx") if $debug;

  ## pick the lowest taxon with the single element
  my @cc = ('s','g','f','o','c','p');
  for my $c (@cc)
  {
    if ( exists $grTx{$c} && @{$grTx{$c}}==1 )
    {
      $cmin = $c;
      last;
    }
  }
  my @t;
  my $gp; # parent
  my $p;  # parent
  if ( exists $grTx{$cmin} )
  {
    @t = @{$grTx{$cmin}}; # extracting highest level taxons
    $gp = $parent{$t[0]};
    $p = $t[0];
  }
  else
  {
    print "\nERROR: grTx{$cmin} does not exist\n";
    print "Group Taxons\n";
    for my $tx (@taxons)
    {
      my $c = get_tx_prefix($tx);
      my $i = $txInt{$c};
      print "tx: $tx\tc: $c\ti: $i\n";
      if ( $i > $imin )
      {
	$imin = $i;
	$cmin = $c;
      }
    }
    print "imin: $imin\n";
    print "cmin: $cmin\n";
    printArrayTbl(\%grTx, "grTbl");
    print "\n";
    exit;
  }

  printArray(\@t, "t=grTx{cmin}") if $debug;

  # for (@t)
  # {
  #   $pTbl{$parent{$_}}++ if exists $parent{$_};
  # }
  # ##$p = $parent{$t[0]};# its enough to take the parent of the first taxonomic rank as all taxonomic ranks at the lowest level have to have the same parent

  # printTbl(\%pTbl, "pTbl") if $debug;

  # for my $p (keys %pTbl)
  # {
  #   push @ch, keys %{$children{$p}};
  # }

  @ch = keys %{$children{$gp}};

  printArray(\@ch, "ch") if $debug;

  @t = ($p);
  @d = diff(\@ch, \@t);

  printArray(\@d, "d") if $debug;

  if (@d == 0)
  {
    print "\nERROR: Difference between ch and t is zero\n";
    print "Taxons:\n";
    for my $tx (@taxons)
    {
      print "$tx\n";
    }
    print "\n";

    if (@taxons > 1)
    {
      printTbl(\%pTbl, "pTbl");
      print "\n";
    }
    else
    {
      my $p = $parent{$rTaxons->[0]};
      print "p: $p\n";
    }

    print "ch: @ch\n";
    print "d: @d\n";
    exit;
  }

  return @d;
}

# get families of the given taxonomic rank
sub get_families
{
  my ($tx, $rParent, $rChildren) = @_;

  my %parent = %{$rParent};
  my %children = %{$rChildren};

  my $c = get_tx_prefix($tx);

  my %txInt;
  $txInt{'p'} = 0;
  $txInt{'c'} = 1;
  $txInt{'o'} = 2;
  $txInt{'f'} = 3;
  $txInt{'g'} = 4;
  $txInt{'s'} = 5;

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


## given a taxon name of the form f_Lactobacillaceae
## extract the taxon prefix, in this example 'f'
sub get_tx_prefix
{
  my $tx = shift;
  my @f = split "_", $tx;
  my $prefix = shift @f;
  return $prefix;
}

# print elements of a hash table
sub printTbl
{
  my ($rTbl, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
  print "\n";
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

sub touchFile
{
  my $file = shift;

  my $cmd = "touch $file";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
exit;
