#!/usr/bin/env perl

=head1 NAME

  taxonomy_cleanup.pl

=head1 DESCRIPTION

  Given a phylogenetic tree with outgroup sequences and a lineage data, the
  script runs vicut clustering based on the taxonomic information, using
  different of majority vote taxonomy modification to clean up taxonomy at all
  taxonomic ranks up to the class level.

  The alignment data (especially the trimmed one) is maintained for the purpose
  truncation and generation of variable region reference db. As such MSAs have
  to curry with them outgroup sequences and hence there has to be a version
  of the lineage file that contains OG info.

  One possibility would be to maintain an outgroup lineage file and a separate
  target group lineage file that would not contain any outgroup sequences.

  Outgroup sequences are kept for the species taxonomy cleanup to identify
  sequences of the target group that cluster with outgroup sequences. These
  sequences will be removed from the target db as errors.

=head1 SYNOPSIS

  taxonomy_cleanup.pl -i <input group name>

=head1 OPTIONS

=over

=item B<--input-group-name, -i>
  Prefix of input group. For example, if Firmicutes_group_6_dir is a directory
  Firmicutes_group_6 group, then the input group name is "Firmicutes_group_6".

=item B<--use-long-spp-names>
  Use long species names for sequences from multi-species vicut clusters.

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

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Firmicutes_dir

  taxonomy_cleanup.pl --debug -i Firmicutes_group_6

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Cwd qw(abs_path);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "input-group|i=s" 	=> \my $grPrefix,
  "rm-OGs-from-tx|r"    => \my $rmOGfromTx,
  "skip-FastTree"       => \my $skipFastTree,
  "use-long-spp-names"  => \my $useLongSppNames,
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
  exit;
}

if ( !$grPrefix )
{
  print "\n\nERROR: Missing input group name.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my $mothur = "/Users/pgajer/bin/mothur";

if ( defined $igs )
{
  $mothur = "/usr/local/packages/mothur-1.36.1/mothur";
}

if ( defined $johanna )
{
  $mothur = "/Users/pgajer/bin/mothur";
}

## export LD_LIBRARY_PATH=/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

my $useLongSppNamesStr = "";
if ($useLongSppNames)
{
  $useLongSppNamesStr = "--use-long-spp-names";
}
####################################################################
##                               MAIN
####################################################################
my $debugStr = "";
if ($debug)
{
  $debugStr = "--debug";
}

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  exit;
}

chdir $grDir;
print "--- Changed dir to $grDir\n";

my $lineageFile     = $grPrefix . ".lineage";
my $algnFile	    = $grPrefix . "_algn.fa";
my $trimmedAlgnFile = $grPrefix . "_algn_trimmed.fa";
my $outgroupFile    = $grPrefix . "_outgroup.seqIDs";
my $treeFile	    = $grPrefix . ".tree";
my $txFile          = $grPrefix . ".tx";

if ( ! -f $lineageFile )
{
  warn "ERROR: $lineageFile does not exist\n";
  exit;
}
elsif ( ! -f $algnFile )
{
  warn "ERROR: $algnFile does not exist\n";
  exit;
}
elsif ( ! -f $trimmedAlgnFile )
{
  warn "WARNING: $trimmedAlgnFile does not exist. Creating a symbolic link to $algnFile.\n";
  my $ap = abs_path( $algnFile );
  #print "ap: $ap\n"; exit;
  my $cmd = "ln -s $ap $trimmedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
elsif ( ! -f $treeFile )
{
  warn "ERROR: $treeFile does not exist\n";
  exit;
}
elsif ( ! -f $outgroupFile )
{
  warn "ERROR: $treeFile does not exist\n";
  exit;
}


## Gathering outgroup data
print "--- Parsing $outgroupFile\n";
my @ogSeqIDs = readArray($outgroupFile);
my %ogInd = map{$_ =>1} @ogSeqIDs; # outgroup elements indicator table

print "--- Extracting seq IDs from trimmed alignment fasta file\n";
my @seqIDs = get_seqIDs_from_fa($trimmedAlgnFile);

print "--- Parsing lineage table\n";
my %lineageTbl = read2colTbl($lineageFile);

# if ($debug)
# {
#   print "\n\n--- Testing if S002972785 is in lineage\n";
#   my $bool = exists $lineageTbl{"S002972785"};
#   print "exists lineageTbl{S002972785}: $bool\n\n\n";
#   exit;
# }

## testing if lineage and fa files has the same seq IDs
print "--- Checking if seqIDs of $trimmedAlgnFile and $lineageFile are the same\n";
my @lSeqIDs = keys %lineageTbl;
my @commSeqIDs = comm(\@seqIDs, \@lSeqIDs);
if (@commSeqIDs != @seqIDs || @commSeqIDs != @lSeqIDs)
{
  warn "WARNING: seq IDs of trimmed alignment fasta file and lineage file do not match";
  print "Number of elements in the trimmed alignment file: " . @seqIDs . "\n";
  print "Number of elements in the lineage file: " . @lSeqIDs . "\n";
  print "Number of elements common to the alignment and lineage files: " . @commSeqIDs . "\n";
}


print "--- Testing if outgroup sequences are part of seqIDs\n";
my @ogDiff = diff( \@ogSeqIDs, \@seqIDs );

if ( scalar(@ogDiff) != 0 )
{
  warn "\n\tERROR the following outgroup seq IDs are not in the trimmed alignment file:\n\n";
  printArray(\@ogDiff);
  exit;
}

print "\n\tNumber of seq's in the trimmed alignment/lineage files: " . @seqIDs . "\n";
print "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n\n";


print "--- Testing if $treeFile is consistent with $trimmedAlgnFile\n";
## extracting leaves' IDs
my $treeLeavesFile = "$grPrefix" . "_tree.leaves";
my $cmd = "nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## looking at the difference between leaf IDs and newTx keys
my @treeLeaves = readArray($treeLeavesFile);

my @commTL = comm(\@treeLeaves, \@seqIDs);
if (@commTL != @treeLeaves || @commTL != @seqIDs)
{
  warn "ERROR: seq IDs of $treeFile and $algnFile do not match";
  print "Number of seq IDs of $treeFile: " . @treeLeaves . "\n";
  print "Number of seq IDs of $trimmedAlgnFile: " . @seqIDs . "\n";
  print "Number of common seq IDs elements: " . @commTL . "\n";
  exit;
}

if ($debug)
{
  print "\n\n\tNumber of seq IDs of $treeFile: " . @treeLeaves . "\n";
  print      "\tNumber of seq IDs of $trimmedAlgnFile: " . @seqIDs . "\n";
  print      "\tNumber of common seq IDs elements: " . @commTL . "\n\n";
}


## Keeping outgroup sequences pollutes taxonomy cleanup process while running
## vicut on higher taxonomic ranks. Sometimes it leads to higher taxonomic ranks
## of OG sequences merging with the taxonomic ranks of the target group
## sequences. This obviously is not good, so I am going to split outgroup data
## here so the cleanup is not muddled by the outgroup sequences.

my %ogLineageTbl;
@ogLineageTbl{@ogSeqIDs} = @lineageTbl{@ogSeqIDs};
delete @lineageTbl{@ogSeqIDs};

## lineage for the target group
my %spTbl;
my %geTbl;
my %faTbl;
my %orTbl;
my %clTbl;
my %phTbl;

my %spp;
my %gen;
my %fam;
my %orr;
my %cls;

my %children;
my %parent;

my %spLineage; # species => lineage of the species (recorded as a string)

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

  $spLineage{$sp} = $lineage;

  ##$sp = "s_$sp";
  #$sp .= "_OG" if ( exists $ogInd{$id} );
  $ge = "g_$ge";
  $fa = "f_$fa";
  $or = "o_$or";
  $cl = "c_$cl";
  $ph = "p_$ph";

  $spp{$id} = $sp;
  $gen{$id} = $ge;
  $fam{$id} = $fa;
  $orr{$id} = $or;
  $cls{$id} = $cl;

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

  push @{$spTbl{$sp}}, $id;
  push @{$geTbl{$ge}}, $id;
  push @{$faTbl{$fa}}, $id;
  push @{$orTbl{$or}}, $id;
  push @{$clTbl{$cl}}, $id;
  push @{$phTbl{$ph}}, $id;
}


print "\nTaxonomy BEFORE cleanup\n";
printLineage();

my $summaryStatsFile = "summary_stats.txt";
open my $SRYOUT, ">$summaryStatsFile" or die "Cannot open $summaryStatsFile for writing: $OS_ERROR\n";
printLineageToFile($SRYOUT, "\n\n====== Taxonomy BEFORE cleanup ======\n");


## outgroup seq's lineage
my %ogSpp;
my %ogSpTbl;
for my $id ( keys %ogLineageTbl )
{
  my $lineage = $ogLineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  $sp .= "_OG";
  $ogSpp{$id} = $sp;
  push @{$ogSpTbl{$sp}}, $id;
}


##
## Basic species summary stats
##
print "--- Generating basic species summary stats\n";
my %sppFreq; ## table of number of sequences per species
map { $sppFreq{$_}++ } values %spp;

my %sppFreq2; ## frequency table of species sequence frequencies
map { $sppFreq2{$_}++ } values %sppFreq;

## number of _sp species
my $nSpSpp = 0;
my %genus;
for (keys %sppFreq)
{
  my ($g, $s) = split "_";
  $nSpSpp++ if defined $s && $s eq "sp";
  push @{$genus{$g}}, $_;
}

## Number of genera with only one species and the species being _sp
my $nSpGenera = 0;
my @spGenera;
my %sspSpecies; # singleton _sp species
for my $g (keys %genus)
{
  if ( scalar(@{$genus{$g}})==1 )
  {
    my $sp = $genus{$g}->[0];
    my ($g, $s) = split "_", $sp;
    if (defined $s && $s eq "sp")
    {
      $nSpGenera++;
      push @spGenera, $g;
      $sspSpecies{$sp} = 1;
    }
  }
}

if ( $debug )
{
  print "\n\nGenera and their species:\n";
  for my $g ( sort { scalar(@{$genus{$b}}) <=> scalar(@{$genus{$a}}) } keys %genus)
  {
    print $g . "\n";
    for (@{$genus{$g}})
    {
      print "\t$_\n";
    }
  }

  print "\n\nFrequency of the number of sequences per species:\n";
  my @q = sort { $sppFreq2{$b} <=> $sppFreq2{$a} || $a <=> $b } keys %sppFreq2;
  printFormatedTbl(\%sppFreq2, \@q);

  print "\n\nNumber of sequences per species:\n";
  @q = sort { $sppFreq{$b} <=> $sppFreq{$a} } keys %sppFreq;
  printFormatedTbl(\%sppFreq, \@q);
  print "\n\n";
}

print "--- Generating species ann and query files\n";

# NOTE. We are going to keep outgroup sequences for the species taxonomy cleanup
# to identify sequences of the target group that cluster with outgroup
# sequences. These sequences will be removed from the target db as errors.

my $annFile    = "spp_ann.tx";
my $queryFile  = "spp_query.seqIDs";
my $vicutDir   = "spp_vicut_dir";
my $nQuerySeqs = 0;
my $nAnnSeqs   = 0;
my %querySpp; # species => count of seq's of that species in the query file
open QOUT, ">$queryFile" or die "Cannot open $queryFile for writing: $OS_ERROR\n";
open ANNOUT, ">$annFile" or die "Cannot open $annFile for writing: $OS_ERROR\n";
for my $id ( keys %spp )
{
  my $t = $spp{$id};
  my ($g, $suffix) = split "_", $t;

  if ( defined $suffix && $suffix eq "sp" && !defined $sspSpecies{$t} ) # do not add _sp species coming from genera where they are the only species = this is why !defined $sspSpecies{$t}
  {
    print QOUT "$id\n";
    $nQuerySeqs++;
    $querySpp{$t}++;
  }
  else
  {
    print ANNOUT "$id\t$t\n"; # all OG seq's end up here as their suffices will not be _sp
    $nAnnSeqs++;
  }
}
# Adding outgroup sequences to the annotation file
for my $id ( keys %ogSpp )
{
  print ANNOUT "$id\t" . $ogSpp{$id} . "\n"; # all OG seq's end up here as their suffices will not be _sp
  $nAnnSeqs++;
}
close ANNOUT;
close QOUT;

if ($debug)
{
  print "\n\n\tNumber of annotation seq's: $nAnnSeqs\n";
  print     "\tNumber of query seq's:      $nQuerySeqs\n";
  print     "\tSum:                        " . ($nAnnSeqs + $nQuerySeqs) . "\n";
  print     "\tNumber of leaves: " . @treeLeaves . "\n\n";

  print     "Query species:\n";
  my @q = sort { $querySpp{$b} <=> $querySpp{$a} } keys %querySpp;
  printFormatedTbl(\%querySpp, \@q);
  print "\n\n"
}

print "--- Running vicut\n";
if ($nQuerySeqs)
{
  $cmd = "vicut -t $treeFile -a $annFile -q $queryFile -o $vicutDir";
}
else
{
  $cmd = "vicut -t $treeFile -a $annFile -o $vicutDir";
}
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## update_spp_tx.pl takes vicut clusters and updates taxonomy using majority vote.
## Moreover, if a species is present in more than 1 cluster, it changes
## taxonomies of sequences from the small clusters to <genus>_sp
print "--- Running update_spp_tx.pl\n";
$cmd = "update_spp_tx.pl $debugStr $useLongSppNamesStr -a $txFile -d $vicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## Rerunning vicut with updated taxons giving a chance sequences whose taxonomy
## changed to <genus>_sp to have their taxonomy updated to a proper species
## taxonomy, as now they will be put in the new query file.

my $vicutFinalTx = "$vicutDir/updated.tx";
my $newTxFile = $vicutFinalTx;

my %newSpp = readTbl($vicutFinalTx);

my $cltrFile = "$vicutDir/minNodeCut.cltrs";
if ( ! -f $cltrFile )
{
  warn "\nERROR: $cltrFile does not exist";
  exit;
}

my %id2cltrTbl;  # seqID => cltrID
my %cltrTbl;     # cltrID => ref to array of IDs in the given cluster
open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR\n";
my $header = <IN>;
foreach (<IN>)
{
  chomp;
  my ($id, $cl, $tx) = split /\s+/,$_;
  $id2cltrTbl{$id} = $cl;
  push @{$cltrTbl{$cl}}, $id;
}
close IN;

my %spFreqTbl;
for my $id ( keys %newSpp )
{
  my $sp = $newSpp{$id};
  my $cl = $id2cltrTbl{$id};
  $spFreqTbl{$sp}{$cl}++;
}

my %sspSeqID; # seqID => 1 if seqID is in the largest cluster of a singleton species;
for my $sp (keys %spFreqTbl)
{
  if (defined $sspSpecies{$sp})
  {
    my @cls = sort { $spFreqTbl{$sp}{$b} <=> $spFreqTbl{$sp}{$a} } keys %{$spFreqTbl{$sp}};
    if ( @cls > 1 && $spFreqTbl{$sp}{$cls[0]} > $spFreqTbl{$sp}{$cls[1]} )
    {
      print "\n\nProcessing $sp\tnCltrs: " . @cls . "\n" if $debug;
      print "Cluster sizes: " if $debug;
      map { print $spFreqTbl{$sp}{$_} . ", "} @cls  if $debug;

      my @f = split "_", $sp;
      my $g = shift @f;
      #my $s = shift @f;
      #print "Genus: $g\n";
      my $clmax = shift @cls;
      for ( @{$cltrTbl{$clmax}} )
      {
	$sspSeqID{$_} = 1;
      }
    } # end of if ( @cls > 1
  }
}
print "\n\n" if $debug;


my @query2;
my @ann2;
for my $id ( keys %newSpp )
{
  my $t = $newSpp{$id};
  my @f = split "_", $t;
  my $g = shift @f;
  my $suffix = shift @f;
  my $suffix2 = shift @f;

  if ( defined $suffix && $suffix eq "sp" )
  {
    if ( !defined $sspSeqID{$id} && !defined $suffix2 )
    {
      push @query2, "$id\n";
      #print "Query2: $id\t$t\n" if $debug;
    }
    else
    {
      push @ann2, "$id\t$t\n";
    }
  }
  # elsif ( $g eq "Unclassified" )
  # {
  #   push @query2, "$id\n";
  #   print "\n\n\tWARNING: $g eq Unclassified in Query2b: $id\t$t\n" if $debug;
  # }
  else
  {
    push @ann2, "$id\t$t\n";
  }
}

my $queryFile2 = "spp_query2.seqIDs";
my $annFile2   = "spp_ann2.tx";
## If _sp sequences, print to file, and set vicutDir to 2nd run directory
if (@query2)
{
  $vicutDir  = "spp_vicut_dir2";
  open QOUT, ">$queryFile2" or die "Cannot open $queryFile2 for writing: $OS_ERROR\n";
  print QOUT @query2;
  close QOUT;

  open ANNOUT, ">$annFile2" or die "Cannot open $annFile2 for writing: $OS_ERROR\n";
  print ANNOUT @ann2;
  close ANNOUT;

  print "--- Running vicut again\n";
  $cmd = "vicut -t $treeFile -a $annFile2 -q $queryFile2 -o $vicutDir";
  ##$cmd = "vicut -t $treeFile -a $annFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Running update_spp_tx.pl\n";
  $cmd = "update_spp_tx.pl $debugStr $useLongSppNamesStr -a $vicutFinalTx -d $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $vicutFinalTx = "$vicutDir/updated.tx";
  %newSpp = readTbl($vicutFinalTx);
}


#printTbl(\%ogInd, "ogInd") if $debug;

my @extraOG; # array of seqIDs that had their taxonomy changed to OG - they will be removed
for my $id (keys %newSpp)
{
  next if exists $ogInd{$id};

  my $sp = $newSpp{$id};

  if (!exists $spLineage{$sp})
  {
    if ( $sp =~ /OG$/)
    {
      print "\n\tWARNING: non-OG seq $id was changed to OG species => $sp; Scheduling removal of the sequence.\n";
      print "\tlineageTbl{$id}: " . $lineageTbl{$id} . "\n\n";
      push @extraOG, $id;
      delete $lineageTbl{$id};
      delete $newSpp{$id};
      next;
    }
    else
    {
      warn "\n\n\tERROR: $id => $sp not found in spLineage";
      print "\tlineageTbl{$id}: " . $lineageTbl{$id} . "\n";
      print "\n\n";
      exit;
    }
  }
  $lineageTbl{$id} = $spLineage{$sp};
}

if ( @extraOG )
{
  # removing elements of extraOG from the updated.tx file
  my $newTxFileNoTGTs = "$vicutDir/updated.tx";
  print "--- Removing elements of extraOG from $newTxFileNoTGTs\n";
  my $origNewTxFileNoTGTs = "$vicutDir/updated_orig.tx";
  my $extraOGfile = "$vicutDir/extraOG.seqIDs";
  writeArray(\@extraOG, $extraOGfile);

  $cmd = "mv $newTxFileNoTGTs $origNewTxFileNoTGTs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "select_tx.pl -e $extraOGfile -i $origNewTxFileNoTGTs -o $newTxFileNoTGTs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # pruning alignment
  print "--- Pruning extraOG seq's from $trimmedAlgnFile\n";
  my $origAlgnFile = $grPrefix . "_algn_trimmed_orig.fa";
  $cmd = "mv $trimmedAlgnFile $origAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "select_seqs.pl --quiet -e $extraOGfile -i $origAlgnFile -o $trimmedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # pruning tree
  print "--- Pruning extraOG seq's from $treeFile\n";
  my $origTreeFile = $grPrefix . "_orig.tree";
  $cmd = "mv $treeFile $origTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "nw_prune $origTreeFile @extraOG > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## Removing outgroup sequences from the updated.tx file and the phylogenetic
## tree.

## Deleting outgroup sequences from %newSpp !!!!
delete @newSpp{@ogSeqIDs};

my %newTxNoTGTs = %newSpp; ## readTbl($newTxFileNoTGTs);

# removing OG seq's from the updated.tx file
my $newTxFileNoTGTs = "$vicutDir/updated.tx";
print "--- Removing OG seq's from $newTxFileNoTGTs\n";
my $origNewTxFileNoTGTs = "$vicutDir/updated_orig2.tx";
my $ogFile = "$vicutDir/og.seqIDs";
writeArray(\@ogSeqIDs, $ogFile);

$cmd = "mv $newTxFileNoTGTs $origNewTxFileNoTGTs";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "select_tx.pl -e $ogFile -i $origNewTxFileNoTGTs -o $newTxFileNoTGTs";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# pruning tree
print "--- Pruning OG seq's from $treeFile\n";
my $origTreeFile = $grPrefix . "_orig2.tree";
$cmd = "mv $treeFile $origTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "nw_prune $origTreeFile @ogSeqIDs > $treeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## Testing consistency between the phylo tree and sequences after taxonomy
## cleanup. In particular, the set of leave seqIDs has to be equal to the set of
## seqIDs of the new taxonomy.

## extracting leave IDs
$treeLeavesFile = "$grPrefix" . "_tree.leaves";
$cmd = "nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my %newTx = readTbl($newTxFileNoTGTs);

## looking at the difference between leaf IDs and newTx keys
@treeLeaves = readArray($treeLeavesFile);
my @survivedIDs = keys %newTx;
my @lostLeaves = diff(\@treeLeaves, \@survivedIDs);

## prunning
if (@lostLeaves>0)
{
  $origTreeFile = $grPrefix . "_orig3.tree";
  $cmd = "mv $treeFile $origTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "\n\tSpecies cleanup eliminated " . @lostLeaves . " sequences\n";
  print "--- Pruning lost seqIDs from the current phylo tree\n";
  $cmd = "nw_prune $origTreeFile @lostLeaves > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


print "--- Running genotype_spp.pl\n";
$cmd = "genotype_spp.pl -d $vicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$newTxFile = "$vicutDir/updated2.tx";
%newTx = readTbl($newTxFile);
#delete @newTx{@ogSeqIDs};

print "--- Updating lineage using new species taxonomy\n";
## updating lineage

my @newTxKeys      = keys %newTx;
my @lineageTblKeys = keys %lineageTbl;

my @d1 = diff(\@newTxKeys, \@lineageTblKeys);
if (@d1)
{
  delete @newTx{@d1};
}

my @d2 = diff(\@lineageTblKeys, \@newTxKeys);
if (@d2)
{
  delete @lineageTbl{@d2};
}

## Testing again consistency between the phylo tree and sequences after taxonomy cleanup.

## extracting leave IDs
$cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## looking at the difference between leaf IDs and newTx keys
@treeLeaves = readArray($treeLeavesFile);
@survivedIDs = keys %newTx;
@lostLeaves = diff(\@treeLeaves, \@survivedIDs);

## prunning
if (@lostLeaves>0)
{
  $origTreeFile = $grPrefix . "_orig4.tree";
  $cmd = "mv $treeFile $origTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "\n\tSpecies cleanup eliminated " . @lostLeaves . " sequences\n";
  print "--- Pruning lost seqIDs from the current phylo tree\n";
  $cmd = "nw_prune $origTreeFile @lostLeaves > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


for my $id (keys %lineageTbl)
{
  #next if exists $ogInd{$id};
  if ( exists $newTx{$id} )
  {
    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    # my $ge = pop @f;
    # my $fa = pop @f;
    # my $or = pop @f;
    # my $cl = pop @f;
    # my $ph = pop @f;
    # "d_Bacteria";

    my $newSp = $newTx{$id};
    my @t = (@f, $newSp);

    $lineageTbl{$id} = join ";", @t;
  }
  else
  {
    delete $lineageTbl{$id};
  }
}

## Testing consistency between lineageTbl and newTx keys
@newTxKeys      = keys %newTx;
@lineageTblKeys = keys %lineageTbl;
my @tlComm = comm(\@newTxKeys, \@lineageTblKeys);
if (@tlComm != @newTxKeys || @tlComm != @lineageTblKeys)
{
  warn "\n\nWARNING: seq IDs of new taxonomy table and the new lineage table do not match";
  print "\n\tNumber of elements in the new taxonomy table: " . @newTxKeys . "\n";
  print "\tNumber of elements in the lineage table: " . @lineageTblKeys . "\n";
  print "\tNumber of common elements: " . @tlComm . "\n";

  writeArray(\@newTxKeys, "newTxKeys.txt");
  writeArray(\@lineageTblKeys, "lineageTblKeys.txt");

  print "\n\tNew taxon and lineage keys written to newTxKeys.txt and lineageTblKeys.txt, respectively\n\n";

  if (@newTxKeys > @lineageTblKeys)
  {
    my @d = diff(\@newTxKeys, \@lineageTblKeys);
    print "\nElements in new taxonomy that are not in new lineage:\n";
    for (@d)
    {
      print "\t$_\t" . $newTx{$_} . "\n";
    }
    print "\n\n";
  }
  elsif (@lineageTblKeys > @newTxKeys)
  {
    my @d = diff(\@lineageTblKeys, \@newTxKeys);
    print "\nElements in new lineage that are not in the new taxonomy:\n";
    for (@d)
    {
      print "\t$_\t" . $lineageTbl{$_} . "\n";
    }
    print "\n\n";
  }

  exit;
}


my %sppFreqFinal; ## table of number of sequences per species
map { $sppFreqFinal{$_}++ } values %newTx;

## number of _sp species
my $nSpSppFinal = 0;
my %genusFinal;
for (keys %sppFreqFinal)
{
  my ($g, $s, $s2) = split "_";
  if (!defined $s2 || (defined $s2 && $s2 ne "OG") )
  {
    $nSpSppFinal++ if defined $s && $s eq "sp";
    push @{$genusFinal{$g}}, $_;
  }
}

## Number of genera with only one species and the species being _sp
my $nSpGeneraFinal = 0;
my @spGeneraFinal;
my %sspSpeciesFinal; # singleton _sp species
for my $g (keys %genusFinal)
{
  if ( scalar(@{$genusFinal{$g}})==1 )
  {
    my $sp = $genusFinal{$g}->[0];
    my ($g, $s) = split "_", $sp;
    if (defined $s && $s eq "sp")
    {
      $nSpGeneraFinal++;
      push @spGeneraFinal, $g;
      $sspSpeciesFinal{$sp} = 1;
    }
  }
}

my %sppFreqFinal2; ## frequency table of species sequence frequencies
map { $sppFreqFinal2{$_}++ } values %sppFreqFinal;


my $finalTxFile = $grPrefix . "_final.tx";
# if ($rmOGfromTx) # if OG seq's are to be removed from the final tx file
# {
#   print "   Removing outgroup sequences from final taxonomy file, as requested";
#   open OUT, ">$finalTxFile" or die "Cannot open $finalTxFile for writing: $OS_ERROR\n";
#   for my $id (keys %newTx)
#   {
#     if ( ! exists $ogInd{$id} )
#     {
#       print OUT "$id\t" . $newTx{$id} . "\n";
#     }
#     else
#     {
#       print "Outgroup sequence $id (" . $newTx{$id} . ") excluded from $finalTxFile\n";
#     }
#   }
#   close OUT;
# }
# else
# {
my $ap = abs_path( $newTxFile );
$cmd = "rm -f $finalTxFile; ln -s $ap $finalTxFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
#}

#%newTx = readTbl($finalTxFile);

## final alignment
# my $finalAlgnFile = $grPrefix . "_algn_trimmed_final.fa";

# my $finalnSeqs  = keys %newTx;
# my $nSeqs       = keys %spp;

# if ($nSeqs != $finalnSeqs)
# {
#   print "\nInitial number of sequences: $nSeqs\n";
#   print "Final number of sequences: $finalnSeqs\n\n";

#   print "--- Generating final alignment files (removing sequences lost in taxonomy cleanup).\n";
#   print "    NOTE: No outgroup sequences here !\n\n";

#   $cmd = "cut -f1 $finalTxFile > final.seqIDs";
#   print "\tcmd=$cmd\n" if $dryRun || $debug;
#   system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

#   $cmd = "select_seqs.pl --quiet -s final.seqIDs -i $trimmedAlgnFile -o $finalAlgnFile";
#   print "\tcmd=$cmd\n" if $dryRun || $debug;
#   system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
# }
# else
# {
#   my $ap = abs_path( $trimmedAlgnFile );
#   $cmd = "rm -f $finalAlgnFile; ln -s $ap $finalAlgnFile";
#   print "\tcmd=$cmd\n" if $dryRun || $debug;
#   system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
# }

## comparison between old and new species assignments
print "--- Comparing old and new species assignments\n";
$cmd = "cmp_tx.pl -i $txFile -j $newTxFile -o old_vs_new_spp_cmp";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


##
## phylo tree with final taxonomy
##

print "\n--- Generating a tree with final species names at leaves\n";
my $finalSppTreeFile = "$grPrefix" . "_final_spp.tree";
$cmd = "rm -f $finalSppTreeFile; nw_rename $treeFile $newTxFile | nw_order - > $finalSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Generating tree with <final species name>_<seqID> labels at leaves\n";
my $finalSppSeqIDsFile = "$grPrefix" . "_final_spp.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalTxFile > $finalSppSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalSppSeqIdTreeFile = "$grPrefix" . "_final_sppSeqIDs.tree";
$cmd = "rm -f $finalSppSeqIdTreeFile; nw_rename $treeFile $finalSppSeqIDsFile | nw_order -  > $finalSppSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final species clades collapsed to a single node \n";
my $finalCondSppTreeFile = "$grPrefix" . "_final_spp_condensed.tree";
$cmd = "rm -f $finalCondSppTreeFile; nw_condense $finalSppTreeFile > $finalCondSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


##
## genus-level cleanup
##
print "--- Generating genus ann and query files\n";

my $genusTxFile = "$grPrefix" . ".genus";
print "--- Writing genus assignments to $genusTxFile file\n";
writeTbl(\%gen, $genusTxFile);


my %id2genTbl;  # seqID => genus
my %annGenera;  # genus => count of seq's of that genus
$annFile = "genus_ann.tx";
open ANNOUT, ">$annFile" or die "Cannot open $annFile for writing: $OS_ERROR\n";
for my $id ( keys %newTx )
{
  my ($g, $suffix) = split "_", $newTx{$id};
  $id2genTbl{$id} = $g;
  print ANNOUT "$id\t$g\n";
  $annGenera{$g}++;
  if (!$g)
  {
    warn "WARNING: undefined genus for $newTx{$id}";
  }
}
close ANNOUT;

if ($debug)
{
  print "\nAnnotation Genera:\n";
  my @q = sort { $annGenera{$b} <=> $annGenera{$a} } keys %annGenera;
  printFormatedTbl(\%annGenera, \@q);
  print "\n\n"
}


print "\n--- Running genus vicut on pruned (after spp cleanup) tree\n";
my $genusVicutDir = "genus_vicut_dir";
$cmd = "vicut -t $treeFile -a $annFile -o $genusVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running update_tx.pl for genus level\n";
$cmd = "update_tx.pl $debugStr -a $genusTxFile -d $genusVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running genotype_tx.pl for genus level\n";
$cmd = "genotype_tx.pl -d $genusVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalGenusTx = "$genusVicutDir/updated2.tx";
my %genusTx = readTbl($finalGenusTx);

## updating lineage
for my $id (keys %genusTx)
{
  next if exists $ogInd{$id};

  my $lineage = $lineageTbl{$id};
  if (!defined $lineage)
  {
    warn "\n\n\tERROR: lineage not found for $id";
    print "\n\n";
    exit;
  }
  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $ge = pop @f;
  # my $fa = pop @f;
  # my $or = pop @f;
  # my $cl = pop @f;
  # my $ph = pop @f;
  # "d_Bacteria";

  my $newGenus = $genusTx{$id};
  my @t = (@f, $newGenus, $sp);

  $lineageTbl{$id} = join ";", @t;
}

print "\n--- Generating a tree with final genus names at leaves\n";
my $finalgeTreeFile = "$grPrefix" . "_final_genus.tree";
$cmd = "rm -f $finalgeTreeFile; nw_rename $treeFile $finalGenusTx | nw_order - > $finalgeTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating tree with <final genus name>_<seqID> labels at leaves\n";
my $finalgeSeqIDsFile = "$grPrefix" . "_final_genus.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalGenusTx > $finalgeSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalgeSeqIdTreeFile = "$grPrefix" . "_final_geSeqIDs.tree";
$cmd = "rm -f $finalgeSeqIdTreeFile; nw_rename $treeFile $finalgeSeqIDsFile | nw_order -  > $finalgeSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final genera collapsed to a single node \n";
my $finalCondgeTreeFile = "$grPrefix" . "_final_genus_condensed.tree";
$cmd = "rm -f $finalCondgeTreeFile; nw_condense $finalgeTreeFile > $finalCondgeTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;



##
## family-level cleanup
##
print "--- Generating family ann and query files\n";

my $familyTxFile = "$grPrefix" . ".family";
print "--- Writing family assignments to $familyTxFile file\n";
writeTbl(\%fam, $familyTxFile);

my %annFamilies;  # family => count of seq's of that family
$annFile = "family_ann.tx";
open ANNOUT, ">$annFile" or die "Cannot open $annFile for writing: $OS_ERROR\n";
for my $id ( keys %newTx )
{
  my $f = $fam{$id};
  if (!$f)
  {
    warn "\n\n\tWARNING: undefined family $id";
    print "\n\n";
  }
  $f =~ s/^f_//;
  print ANNOUT "$id\t$f\n";
  $annFamilies{$f}++;
}
close ANNOUT;

if ($debug)
{
  print "\nAnnotation Families:\n";
  my @q = sort { $annFamilies{$b} <=> $annFamilies{$a} } keys %annFamilies;
  printFormatedTbl(\%annFamilies, \@q);
  print "\n\n"
}


print "\n--- Running family vicut on pruned (after spp cleanup) tree\n";
my $familyVicutDir = "family_vicut_dir";
$cmd = "vicut -t $treeFile -a $annFile -o $familyVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running update_tx.pl for family level\n";
$cmd = "update_tx.pl $debugStr -a $familyTxFile -d $familyVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running genotype_tx.pl for family level\n";
$cmd = "genotype_tx.pl -d $familyVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalFamilyTx = "$familyVicutDir/updated2.tx";
my %familyTx = readTbl($finalFamilyTx);

## updating lineage
for my $id (keys %familyTx)
{
  next if exists $ogInd{$id};

  my $lineage = $lineageTbl{$id};
  if (!defined $lineage)
  {
    warn "\n\n\tERROR: lineage not found for $id";
    print "\n\n";
    exit;
  }

  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  # my $or = pop @f;
  # my $cl = pop @f;
  # my $ph = pop @f;
  # "d_Bacteria";

  my $newFamily = $familyTx{$id};
  my @t = (@f, $newFamily, $ge, $sp);

  $lineageTbl{$id} = join ";", @t;
}

print "\n--- Generating a tree with final family names at leaves\n";
my $familyTreeFile = "$grPrefix" . "_final_family.tree";
$cmd = "rm -f $familyTreeFile; nw_rename $treeFile $finalFamilyTx | nw_order - > $familyTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating tree with <final family name>_<seqID> labels at leaves\n";
my $familySeqIDsFile = "$grPrefix" . "_final_family.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalFamilyTx > $familySeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $familySeqIdTreeFile = "$grPrefix" . "_final_family_seqIDs.tree";
$cmd = "rm -f $familySeqIdTreeFile; nw_rename $treeFile $familySeqIDsFile | nw_order -  > $familySeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final families collapsed to a single node \n";
my $familyCondTreeFile = "$grPrefix" . "_final_family_condensed.tree";
$cmd = "rm -f $familyCondTreeFile; nw_condense $familyTreeFile > $familyCondTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;




##
## order-level cleanup
##
print "--- Generating order ann and query files\n";

my $orderTxFile = "$grPrefix" . ".order";
print "--- Writing order assignments to $orderTxFile file\n";
writeTbl(\%orr, $orderTxFile);

my %annOrders;  # order => count of seq's of that order
$annFile = "order_ann.tx";
open ANNOUT, ">$annFile" or die "Cannot open $annFile for writing: $OS_ERROR\n";
for my $id ( keys %newTx )
{
  my $f = $orr{$id};
  if (!$f)
  {
    warn "\n\n\tWARNING: undefined order for $id";
    print "\n\n";
  }
  $f =~ s/^o_//;
  print ANNOUT "$id\t$f\n";
  $annOrders{$f}++;
}
close ANNOUT;

if ($debug)
{
  print "\nAnnotation Orders:\n";
  my @q = sort { $annOrders{$b} <=> $annOrders{$a} } keys %annOrders;
  printFormatedTbl(\%annOrders, \@q);
  print "\n\n"
}


print "\n--- Running order vicut on pruned (after spp cleanup) tree\n";
my $orderVicutDir = "order_vicut_dir";
$cmd = "vicut -t $treeFile -a $annFile -o $orderVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running update_tx.pl for order level\n";
$cmd = "update_tx.pl $debugStr -a $orderTxFile -d $orderVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running genotype_tx.pl for order level\n";
$cmd = "genotype_tx.pl -d $orderVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalOrderTx = "$orderVicutDir/updated2.tx";
my %orderTx = readTbl($finalOrderTx);

## updating lineage
for my $id (keys %orderTx)
{
  next if exists $ogInd{$id};

  my $lineage = $lineageTbl{$id};
  if (!defined $lineage)
  {
    warn "\n\n\tERROR: lineage not found for $id";
    print "\n\n";
    exit;
  }

  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  my $or = pop @f;
  # my $cl = pop @f;
  # my $ph = pop @f;
  # "d_Bacteria";

  my $newOrd = $orderTx{$id};
  my @t = (@f, $newOrd, $fa, $ge, $sp);

  $lineageTbl{$id} = join ";", @t;
}

print "\n--- Generating a tree with final order names at leaves\n";
my $orderTreeFile = "$grPrefix" . "_final_order.tree";
$cmd = "rm -f $orderTreeFile; nw_rename $treeFile $finalOrderTx | nw_order - > $orderTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating tree with <final order name>_<seqID> labels at leaves\n";
my $orderSeqIDsFile = "$grPrefix" . "_final_order.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalOrderTx > $orderSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $orderSeqIdTreeFile = "$grPrefix" . "_final_order_seqIDs.tree";
$cmd = "rm -f $orderSeqIdTreeFile; nw_rename $treeFile $orderSeqIDsFile | nw_order -  > $orderSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final orders collapsed to a single node \n";
my $orderCondTreeFile = "$grPrefix" . "_final_order_condensed.tree";
$cmd = "rm -f $orderCondTreeFile; nw_condense $orderTreeFile > $orderCondTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;



##
## class-level cleanup
##
print "--- Generating class ann and query files\n";

my $classTxFile = "$grPrefix" . ".class";
print "--- Writing class assignments to $classTxFile file\n";
writeTbl(\%orr, $classTxFile);

my %annClasss;  # class => count of seq's of that class
$annFile = "class_ann.tx";
open ANNOUT, ">$annFile" or die "Cannot open $annFile for writing: $OS_ERROR\n";
for my $id ( keys %newTx )
{
  my $f = $cls{$id};
  if (!$f)
  {
    warn "\n\n\tWARNING: undefined class for $id";
    print "\n\n";
  }
  $f =~ s/^c_//;
  print ANNOUT "$id\t$f\n";
  $annClasss{$f}++;
}
close ANNOUT;

if ($debug)
{
  print "\nAnnotation Classes:\n";
  my @q = sort { $annClasss{$b} <=> $annClasss{$a} } keys %annClasss;
  printFormatedTbl(\%annClasss, \@q);
  print "\n\n"
}


print "\n--- Running class vicut on pruned (after spp cleanup) tree\n";
my $classVicutDir = "class_vicut_dir";
$cmd = "vicut -t $treeFile -a $annFile -o $classVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running update_tx.pl for class level\n";
$cmd = "update_tx.pl $debugStr -a $classTxFile -d $classVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running genotype_tx.pl for class level\n";
$cmd = "genotype_tx.pl -d $classVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalClassTx = "$classVicutDir/updated2.tx";
my %classTx = readTbl($finalClassTx);

## updating lineage
for my $id (keys %classTx)
{
  next if exists $ogInd{$id};

  my $lineage = $lineageTbl{$id};
  if (!defined $lineage)
  {
    warn "\n\n\tERROR: lineage not found for $id";
    print "\n\n";
    exit;
  }

  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  my $or = pop @f;
  my $cl = pop @f;
  # my $ph = pop @f;
  # "d_Bacteria";

  my $newCl = $classTx{$id};
  my @t = (@f, $newCl, $or, $fa, $ge, $sp);

  $lineageTbl{$id} = join ";", @t;
}

print "\n--- Generating a tree with final class names at leaves\n";
my $classTreeFile = "$grPrefix" . "_final_class.tree";
$cmd = "rm -f $classTreeFile; nw_rename $treeFile $finalClassTx | nw_order - > $classTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating tree with <final class name>_<seqID> labels at leaves\n";
my $classSeqIDsFile = "$grPrefix" . "_final_class.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalClassTx > $classSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $classSeqIdTreeFile = "$grPrefix" . "_final_class_seqIDs.tree";
$cmd = "rm -f $classSeqIdTreeFile; nw_rename $treeFile $classSeqIDsFile | nw_order -  > $classSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final classs collapsed to a single node \n";
my $classCondTreeFile = "$grPrefix" . "_final_class_condensed.tree";
$cmd = "rm -f $classCondTreeFile; nw_condense $classTreeFile > $classCondTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;



print "--- Generating final lineage files\n";

my $initNumSpecies  = scalar( keys %spTbl );
my $initNumGenera   = scalar( keys %geTbl );
my $initNumFamilies = scalar( keys %faTbl );
my $initNumOrders   = scalar( keys %orTbl );
my $initNumClasses  = scalar( keys %clTbl );
my $initNumPhyla    = scalar( keys %phTbl );

undef %spTbl;
undef %geTbl;
undef %faTbl;
undef %orTbl;
undef %clTbl;
undef %phTbl;
undef %children;
undef %parent;

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
  $sp .= "_OG" if ( exists $ogInd{$id} );
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

  push @{$spTbl{$sp}}, $id;
  push @{$geTbl{$ge}}, $id;
  push @{$faTbl{$fa}}, $id;
  push @{$orTbl{$or}}, $id;
  push @{$clTbl{$cl}}, $id;
  push @{$phTbl{$ph}}, $id;
}



print "\nTaxonomy AFTER cleanup\n";
printLineage();
printLineageToFile($SRYOUT, "\n\n====== Taxonomy AFTER cleanup ======\n");


my $finalLineageFile = $grPrefix . "_final.lineage";
open OUT, ">$finalLineageFile" or die "Cannot open $finalLineageFile for writing: $OS_ERROR\n";
for my $id (keys %newTx)
{
  my $lineage = $lineageTbl{$id};
  print OUT "$id\t$lineage\n";
}
close OUT;


my $finalLineageFile2 = $grPrefix . "_final_no_tGTs.lineage";
open OUT, ">$finalLineageFile2" or die "Cannot open $finalLineageFile2 for writing: $OS_ERROR\n";
for my $id (keys %newTx)
{
  my $lineage = $lineageTbl{$id};

  my @f = split ";", $lineage;
  my $sp = pop @f;
  print OUT "$id\t";
  for (@f)
  {
    print OUT "$_;";
  }
  print OUT $newTxNoTGTs{$id} . "\n";
}
for my $id ( keys %ogSpp )
{
  my $lineage = $ogLineageTbl{$id};
  print OUT "$id\t$lineage\n";
}
close OUT;


## ---------------------------------------
##    Taxonomy Cleanup Summary Stats
## ---------------------------------------

## log file of old and new taxonomic assignments for each sequence
my $oldNewTxFile = "old_new.spp";
open OUT, ">$oldNewTxFile" or die "Cannot open $oldNewTxFile for writing: $OS_ERROR\n";
my %sppTr; ## species transition table
for ( keys %spp )
{
  if ( exists $newTx{$_} )
  {
    print OUT "$_\t$spp{$_}\t$newTx{$_}\n";
    $sppTr{ $spp{$_} }{ $newTx{$_} }++;
  }
  else
  {
    print OUT "$_\t$spp{$_}\tNA\n";
    $sppTr{ $spp{$_} }{ "NA" }++;
  }

}
close OUT;


my $readmeFile = "README";
generateREADMEfile($readmeFile);


print $SRYOUT "\n\n--- Final summary\n";

print $SRYOUT  "\n\tNumber of species (with OG seq's) BEFORE taxonomic cleanup: " . scalar( keys %sppFreq ) . "\n";
print $SRYOUT    "\tNumber of species (with OG seq's) AFTER taxonomic cleanup and tentative ribotype splits:  " . scalar( keys %sppFreqFinal ) . "\n\n";

print $SRYOUT "\tNumber of _sp species BEFORE taxonomic cleanup: $nSpSpp\n";
print $SRYOUT "\tNumber of _sp species AFTER taxonomic cleanup: $nSpSppFinal\n\n";

print $SRYOUT  "\tNumber of singletons species BEFORE taxonomic cleanup: $sppFreq2{1}\n";
print $SRYOUT  "\tNumber of singletons species AFTER taxonomic cleanup: $sppFreqFinal2{1}\n\n";

print $SRYOUT  "\tNumber of genera BEFORE taxonomic cleanup: $initNumGenera\n";
print $SRYOUT  "\tNumber of genera AFTER taxonomic cleanup: " . scalar( keys %geTbl ) . "\n\n";

print $SRYOUT  "\tNumber of genera with only one species and the species being _sp BEFORE taxonomic cleanup: $nSpGenera\n";
print $SRYOUT  "\tNumber of genera with only one species and the species being _sp AFTER taxonomic cleanup: $nSpGeneraFinal\n\n";

print $SRYOUT  "\tNumber of families BEFORE taxonomic cleanup: $initNumFamilies\n";
print $SRYOUT  "\tNumber of families AFTER taxonomic cleanup: " . scalar( keys %faTbl ) . "\n\n";

print $SRYOUT  "\tNumber of orders BEFORE taxonomic cleanup: $initNumOrders\n";
print $SRYOUT  "\tNumber of orders AFTER taxonomic cleanup: " . scalar( keys %orTbl ) . "\n\n";

print $SRYOUT  "\tNumber of classes BEFORE taxonomic cleanup: $initNumClasses\n";
print $SRYOUT  "\tNumber of classes AFTER taxonomic cleanup: " . scalar( keys %clTbl ) . "\n\n";

print $SRYOUT  "\tNumber of phyla BEFORE taxonomic cleanup: $initNumPhyla\n";
print $SRYOUT  "\tNumber of phyla AFTER taxonomic cleanup: " . scalar( keys %phTbl ) . "\n\n";

close $SRYOUT;


print  "\n\n--- Final summary\n";

print   "\n\tNumber of species (with OG seq's) BEFORE taxonomic cleanup: " . scalar( keys %sppFreq ) . "\n";
print     "\tNumber of species (with OG seq's) AFTER taxonomic cleanup and tentative ribotype splits:  " . scalar( keys %sppFreqFinal ) . "\n\n";

print  "\tNumber of _sp species BEFORE taxonomic cleanup: $nSpSpp\n";
print  "\tNumber of _sp species AFTER taxonomic cleanup: $nSpSppFinal\n\n";

print   "\tNumber of singletons species BEFORE taxonomic cleanup: $sppFreq2{1}\n";
print   "\tNumber of singletons species AFTER taxonomic cleanup: $sppFreqFinal2{1}\n\n";

print   "\tNumber of genera BEFORE taxonomic cleanup: $initNumGenera\n";
print   "\tNumber of genera AFTER taxonomic cleanup: " . scalar( keys %geTbl ) . "\n\n";

print   "\tNumber of genera with only one species and the species being _sp BEFORE taxonomic cleanup: $nSpGenera\n";
print   "\tNumber of genera with only one species and the species being _sp AFTER taxonomic cleanup: $nSpGeneraFinal\n\n";

print   "\tNumber of families BEFORE taxonomic cleanup: $initNumFamilies\n";
print   "\tNumber of families AFTER taxonomic cleanup: " . scalar( keys %faTbl ) . "\n\n";

print   "\tNumber of orders BEFORE taxonomic cleanup: $initNumOrders\n";
print   "\tNumber of orders AFTER taxonomic cleanup: " . scalar( keys %orTbl ) . "\n\n";

print   "\tNumber of classes BEFORE taxonomic cleanup: $initNumClasses\n";
print   "\tNumber of classes AFTER taxonomic cleanup: " . scalar( keys %clTbl ) . "\n\n";

print   "\tNumber of phyla BEFORE taxonomic cleanup: $initNumPhyla\n";
print   "\tNumber of phyla AFTER taxonomic cleanup: " . scalar( keys %phTbl ) . "\n\n";


print "\n\n\tSee $grDir/$readmeFile for info about the taxonomy curation algorithm and the content of the output directory.\n";
print     "\tSummary stats written to $grDir/$summaryStatsFile\n\n";

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

sub printF2
{
  my ($nTabs, $key, $val, $rKV, $fh) = @_;

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
  print $fh "$tabs$key$pad$val\n";
}


## print lineage after first round of cleanup (spp level)

sub printLineageAfterSppToFile
{
  my %newChildren;
  my %newSpTbl;
  for my $id (keys %newTx)
  {
    if (!exists $gen{$id})
    {
      warn "\nERROR: gen{$id} does not exist";
      print "newTx{$id}: " . $newTx{$id} . "\n";
    }
    $newChildren{$gen{$id}}{$newTx{$id}}++;
    push @{$newSpTbl{$newTx{$id}}}, $id;
  }

  print $SRYOUT "\nTaxonomy AFTER cleanup\n";
  my @phs = keys %phTbl;
  my $phKVL = getKeyValStrLengths(\%phTbl);
  for my $ph (@phs)
  {
    printF2(0, $ph, scalar(@{$phTbl{$ph}}), $phKVL, $SRYOUT);
    my @cls = keys %{$children{$ph}};
    my $clKVL = getKeyValStrLengths(\%clTbl, \@cls);
    for my $cl ( sort{scalar(@{$clTbl{$a}}) <=> scalar(@{$clTbl{$b}})} @cls)
    {
      printF2(1, $cl, scalar(@{$clTbl{$cl}}), $clKVL, $SRYOUT);
      my @ors = keys %{$children{$cl}};
      my $orKVL = getKeyValStrLengths(\%orTbl, \@ors);
      for my $or ( sort{scalar(@{$orTbl{$a}}) <=> scalar(@{$orTbl{$b}})} @ors)
      {
	printF2(2, $or, scalar(@{$orTbl{$or}}), $orKVL, $SRYOUT);
	my @fas = keys %{$children{$or}};
	my $faKVL = getKeyValStrLengths(\%faTbl, \@fas);
	for my $fa ( sort{scalar(@{$faTbl{$a}}) <=> scalar(@{$faTbl{$b}})} @fas)
	{
	  printF2(3, $fa, scalar(@{$faTbl{$fa}}), $faKVL, $SRYOUT);
	  my @ges = keys %{$children{$fa}};
	  my $geKVL = getKeyValStrLengths(\%geTbl, \@ges);
	  for my $ge ( sort{scalar(@{$geTbl{$a}}) <=> scalar(@{$geTbl{$b}})} @ges)
	  {
	    printF2(4, $ge, scalar(@{$geTbl{$ge}}), $geKVL, $SRYOUT);
	    my @sps = keys %{$newChildren{$ge}};
	    my $spKVL = getKeyValStrLengths(\%newSpTbl, \@sps);
	    for my $sp ( sort @sps)
	    {
	      printF2(5, $sp, $newChildren{$ge}{$sp}, $spKVL, $SRYOUT);
	    }
	  }
	}
      }
    }
  }
}

sub printLineageAfterSpp
{
  my %newChildren;
  my %newSpTbl;
  for my $id (keys %newTx)
  {
    if (!exists $gen{$id})
    {
      warn "\nERROR: gen{$id} does not exist";
      print "newTx{$id}: " . $newTx{$id} . "\n";
    }
    $newChildren{$gen{$id}}{$newTx{$id}}++;
    push @{$newSpTbl{$newTx{$id}}}, $id;
  }

  my @phs = keys %phTbl;
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
	  my @ges = keys %{$children{$fa}};
	  my $geKVL = getKeyValStrLengths(\%geTbl, \@ges);
	  for my $ge ( sort{scalar(@{$geTbl{$a}}) <=> scalar(@{$geTbl{$b}})} @ges)
	  {
	    printF(4, $ge, scalar(@{$geTbl{$ge}}), $geKVL);
	    if (!exists $newChildren{$ge})
	    {
	      warn "ERROR: newChildren{$ge} does not exist";
	      exit;
	    }
	    my @sps = keys %{$newChildren{$ge}};
	    my $spKVL = getKeyValStrLengths(\%newSpTbl, \@sps);
	    for my $sp ( sort @sps)
	    {
	      printF(5, $sp, $newChildren{$ge}{$sp}, $spKVL);
	    }
	  }
	}
      }
    }
  }
}

## print lineage
sub printLineage
{
  my @phs = keys %phTbl;
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
	  my @ges = keys %{$children{$fa}};
	  my $geKVL = getKeyValStrLengths(\%geTbl, \@ges);
	  for my $ge ( sort{scalar(@{$geTbl{$a}}) <=> scalar(@{$geTbl{$b}})} @ges)
	  {
	    printF(4, $ge, scalar(@{$geTbl{$ge}}), $geKVL);
	    my @sps = keys %{$children{$ge}};
	    my $spKVL = getKeyValStrLengths(\%spTbl, \@sps);
	    for my $sp ( sort{scalar(@{$spTbl{$a}}) <=> scalar(@{$spTbl{$b}})} @sps)
	    {
	      printF(5, $sp, scalar(@{$spTbl{$sp}}), $spKVL);
	    }
	  }
	}
      }
    }
  }
}

sub printLineageToFile
{
  my ($fh, $msg) = @_;

  print $fh "\n$msg\n";

  my @phs = keys %phTbl;
  my $phKVL = getKeyValStrLengths(\%phTbl);

  for my $ph (@phs)
  {
    printF2(0, $ph, scalar(@{$phTbl{$ph}}), $phKVL, $fh);
    my @cls = keys %{$children{$ph}};
    my $clKVL = getKeyValStrLengths(\%clTbl, \@cls);
    for my $cl ( sort{scalar(@{$clTbl{$a}}) <=> scalar(@{$clTbl{$b}})} @cls)
    {
      printF2(1, $cl, scalar(@{$clTbl{$cl}}), $clKVL, $fh);
      my @ors = keys %{$children{$cl}};
      my $orKVL = getKeyValStrLengths(\%orTbl, \@ors);
      for my $or ( sort{scalar(@{$orTbl{$a}}) <=> scalar(@{$orTbl{$b}})} @ors)
      {
	printF2(2, $or, scalar(@{$orTbl{$or}}), $orKVL, $fh);
	my @fas = keys %{$children{$or}};
	my $faKVL = getKeyValStrLengths(\%faTbl, \@fas);
	for my $fa ( sort{scalar(@{$faTbl{$a}}) <=> scalar(@{$faTbl{$b}})} @fas)
	{
	  printF2(3, $fa, scalar(@{$faTbl{$fa}}), $faKVL, $fh);
	  my @ges = keys %{$children{$fa}};
	  my $geKVL = getKeyValStrLengths(\%geTbl, \@ges);
	  for my $ge ( sort{scalar(@{$geTbl{$a}}) <=> scalar(@{$geTbl{$b}})} @ges)
	  {
	    printF2(4, $ge, scalar(@{$geTbl{$ge}}), $geKVL, $fh);
	    my @sps = keys %{$children{$ge}};
	    my $spKVL = getKeyValStrLengths(\%spTbl, \@sps);
	    for my $sp ( sort{scalar(@{$spTbl{$a}}) <=> scalar(@{$spTbl{$b}})} @sps)
	    {
	      printF2(5, $sp, scalar(@{$spTbl{$sp}}), $spKVL, $fh);
	    }
	  }
	}
      }
    }
  }
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
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

  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  for (@args)
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "$_$pad" . $rTbl->{$_} . "\n";
  }
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

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    print "\n\nERROR in readTbl() at line " . __LINE__ . ": $file does not exist\n\n\n";
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


# read taxon table excluding outgroup seq ID and any seqID with OUTGROUP.* type label
sub readTxTbl{

  my ($file, $outgroupSeqID) = @_;

  if ( ! -f $file )
  {
    print "\n\nERROR in readTxTbl() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( defined $outgroupSeqID )
  {
    foreach (<IN>)
    {
      chomp;
      my ($id, $t) = split /\s+/,$_;

      if ( $t !~ /OUTGROUP.*/ )
      {
	$tbl{$id} = $t;
      }
    }
  }
  else
  {
    foreach (<IN>)
    {
      chomp;
      my ($id, $t) = split /\s+/,$_;

      if ( ($id ne $outgroupSeqID) && ($t !~ /OUTGROUP.*/) )
      {
	$tbl{$id} = $t;
      }
    }
  }
  close IN;

  return %tbl;
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

## parse fasta file and check if any of the seq IDs matches the first argument
sub findSeqIDinFasta
{
  my ($seqId, $inFile) = @_;

  if ( ! -f $inFile )
  {
    ##print "\n\nERROR in findSeqIDinFasta() at line " . __LINE__ . ": $inFile and does not exist\n\n\n";
    print "\n\nERROR in findSeqIDinFasta(): $inFile and does not exist\n\n\n";
    exit;
  }

  my $found = 0;
  open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
  $/ = ">";
  my $junkFirstOne = <FASTA>;
  while (<FASTA>)
  {
    chomp;
    my ($def,@seqlines) = split /\n/, $_;
    my $seq = join '', @seqlines;
    my ($id) = split /\s+/, $def;
    if ( $id eq $seqId )
    {
      $found = 1;
      last;
    }
  }
  $/ = "\n";
  close FASTA;

  return $found;
}


# parse table with row ids (and optional header)
# output written to
# 1. tbl: <rowId> => <column of the table associated with rowId>
# 2. rowIds: array of row ids, in the order as they appear in the input file
# 3. header array: if header present
#
sub readBigTbl{

  my ($file, $hasHeader) = @_;

  if ( ! -f $file )
  {
    print "\n\nERROR in readBigTbl() at line " . __LINE__ . ": $file does not exist\n\n\n";
    exit;
  }

  my %tbl;
  my @header;
  my @rowIds;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";

  if ( defined $hasHeader )
  {
    my $headerStr = <IN>;
    chomp $headerStr;
    @header = split "\t", $headerStr;
    shift @header; # get rid of <rowId> label
  }

  foreach my $line (<IN>)
  {
    chomp $line;
    my @f = split "\t", $line;
    my $id = shift @f;
    push @rowIds, $id;
    $tbl{ $id } = \@f;
  }
  close IN;

  return (\%tbl, \@rowIds, \@header);
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

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

sub generateREADMEfile
{
  my $file = shift;

  my $text = qq~

The algorithm of automatic curation of full length sequences has the following structure

A given taxonomic rank (say phylum or family) a set of sequences from this
taxonomic rank members is chosed from RDP db, so that there are no more than 100
sequences per species. Outgroup sequences are chosen from a sibling taxonomic
rank.

I use SATe to generate a MSA of the input fasta file

A QC step consists of the following components

  mothur > summary.seqs()
  mothur > screen.seqs()
  trim_align.pl
  mothur > filter.seqs()

  screen.seqs() removes sequences using 95 percentile criteria for the start and
  the end of alignment postions and min length of the sequence. Sequences with
  ambiguity codes are also filtered out.

The remaining steps of the algorithm are

   * align outgroup sequence to the input Multiple Sequence Alignment
   * add aligned outgroup to the alignment
   * do the same for the input taxon file
   * build tree; reroot
   * generate annotation and query file; query file contains sequences with
     <genus>_sp taxonomy; annotation file contains the remaining sequences.
   * run vicut with the generated annotation and query files
   * remove outliers and rename query taxonomy based on vicut results; outliers are singleton clusters.
   * repeat the above two steps unit there are no outliers


The taxonomy_cleanup.log file contains the standard output of the taxonomic cleanup command.

At the end is a list of lost taxons and transformation table from old to new
taxon assignment. NA means that sequences were removed in the QC step. For example

Mycoplasma_felis:            NA(6)
Mycoplasma_imitans:          Mycoplasma_gallisepticum(1)

The first line means that 6 sequences of Mycoplasma_felis were removed in the QC
step. There was only one sequence of Mycoplasma_imitans and it was renamed to
Mycoplasma_gallisepticum.

Please check out *_aln.bad.summary file for taxonomy and alignment parameters of
all sequences removed as a result of QC. The last column contains a reason for
removal.

trim_align.pl trims sequences using 95% criterion for the start and end positions.

filter.seqs() removes all-gap columns from the trimmed MSA.

*_nr_SATe_aln_good95_tf_wOG.fa file contains final MSA including the outgroup
 sequence after the QC step and removal of gap columns.

* _nr_SATe_aln.fa contains the initial alignment

Here is a set of phylogenetic trees that are generated by the pipeline

*_nr_SATe_aln_good95_tf_wOG_rr_sppSeqIDs_before_tx_change.tree
*_nr_SATe_aln_good95_tf_wOG_rr_spp_before_tx_change.tree

The above two trees (generated using FastTree) are trees with taxonomy before
taxonomy was changed using vicut.

_spp contains species lables at the leaf nodes
_sppSeqIDs contains species_seqID lables at the leaf nodes

The following trees are generated after the initial taxonomy change

*_nr_SATe_aln_good95_tf_wOG_rr_spp.tree
*_nr_SATe_aln_good95_tf_wOG_rr_sppSeqIDs.tree
*_nr_SATe_aln_good95_tf_wOG_rr_sppCltr.tree
*_nr_SATe_aln_good95_tf_wOG_rr_sppCondensed.tree

_sppCltr contains <species name>_<vicut cluster ID> at the leaf nodes
_sppCondensed is a condencsed _spp tree (clades collapsed to one leaf).

If there are outliers in the tree, they will be removed in the second prunning
step of the algorithm and the corresonding trees are labeled as

*_nr_SATe_aln_good95_tf_wOG_rr_pruned1_spp.tree
*_nr_SATe_aln_good95_tf_wOG_rr_pruned1_sppCltr.tree
*_nr_SATe_aln_good95_tf_wOG_rr_pruned1_sppCondensed.tree
*_nr_SATe_aln_good95_tf_wOG_rr_pruned1_sppSeqIDs.tree

~;

  open OUT, ">$file" or die "Cannot open $file for writing: $OS_ERROR\n";
  print OUT "$text\n";
  close OUT;
}

exit;
