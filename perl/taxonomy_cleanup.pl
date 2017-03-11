#!/usr/bin/env perl

=head1 NAME

  taxonomy_cleanup.pl

=head1 VERSION 0.95

=head1 CHANGE LOG
  Mar 8, 2017:
  - Cleanup of genus level annotations added to lines 841-879

  Mar 5, 2017:
  - In the second run of vicut the way _sp seq's are used has been changed.
    In the first run if an _sp species was the only species of a genus, then
    the corresponding seq's were put in the annotation file. But now in the second
    vicut run, only the largest cluster of that species is treated the same way
    and all other treated as query sequences.

  - Removed the source of Unclassified taxonomy in the second run of vicut. It
    stemmed from the fact that there were leaves in the vicut input tree that
    were neither in the query nor in the annotation files.


  Mar 3, 2017:
  - Lines 445-493: Check if _sp sequences remain after first vicut run. If so,
  	push to array, print to file, and re-run vicut.
  - Line 464: If taxa labeled as "Unclassified", add to query2 file.

  Feb 28, 2017:
  - change input to group name only
  - add chDir to group name directory
  - alter handle assignments to be files within group directory
  - added in "--johanna" flag to direct to mothur location

  Feb 21, 2017: changing interface/input data
   - aligned fasta file
   - lineage file with the same sample IDs as in the alignment file
   - list of outgroup seqIDs

  Getting rid of taxon file.

  Feb 20, 2017: Using lineage file instead of taxon file, so that the hole
  lineage can be updated not only species assignment.

=head1 DESCRIPTION

  Given Multiple Sequence Alignment, lineage file and an outgroup sequence ID,
  QC alignment sequences using 95 percentile criteria for the start and the end
  of alignment postion and minlength. Also, filter out sequences with ambiguity
  codes.

   * align outgroup sequence to the input Multiple Sequence Alignment
   * add aligned outgroup to the alignment
   * do the same for the input taxon file
   * build tree; reroot
   * generate annotation and query file; query file contains sequences with
     <genus>_sp taxonomy (excluding singleton _sp species - the only species of a genus); annotation file contains the remaining sequences.
   * run vicut with the generated annotation and query files
   * remove outliers and rename query taxonomy based on vicut results

  It is also possible to pass as input a tree file.

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

  taxonomy_cleanup.pl -i Firmicutes_group_6

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
  ##$refFile = "/home/pgajer/projects/16S_MC2/rdp_Archaea_Bacteria_all_good_headers.fa";
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
  warn "WARNING: $trimmedAlgnFile does not exist. Creaint a symbolic link to $algnFile.\n";
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


print "--- Testing if $treeFile is consistent with the trimmed alignment file\n";
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
  print "Number of seq IDs of $algnFile: " . @seqIDs . "\n";
  print "Number of common seq IDs elements: " . @commTL . "\n";
  exit;
}

if ($debug)
{
  print "\n\n\tNumber of seq IDs of $treeFile: " . @treeLeaves . "\n";
  print      "\tNumber of seq IDs of $algnFile: " . @seqIDs . "\n";
  print      "\tNumber of common seq IDs elements: " . @commTL . "\n\n";
}

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


print "\nTaxonomy before cleanup\n";
printLineage();

my $summaryStatsFile = "summary_stats.txt";
open my $SRYOUT, ">$summaryStatsFile" or die "Cannot open $summaryStatsFile for writing: $OS_ERROR\n";
printLineageToFile($SRYOUT);


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
  print "\n\nGenera and their species\n";
  for my $g ( sort { scalar(@{$genus{$b}}) <=> scalar(@{$genus{$a}}) } keys %genus)
  {
    print $g . "\n";
    for (@{$genus{$g}})
    {
      print "\t$_\n";
    }
  }

  print "Frequency of the number of sequences per species\n";
  for ( sort { $sppFreq2{$b} <=> $sppFreq2{$a} || $a <=> $b } keys %sppFreq2 )
  {
    print "$_\t$sppFreq2{$_}\n";
  }
  print "\n";


  print "\n\nNumber of sequences per species\n";
  for ( sort { $sppFreq{$b} <=> $sppFreq{$a} } keys %sppFreq )
  {
    print "$_\t$sppFreq{$_}\n";
  }
  print "\n\n";
}

print $SRYOUT "\n\nNumber of species: " . (scalar( keys %sppFreq )-1) . "\n"; # -1 because OUTGROUP is subtracted
print $SRYOUT  "Number of _sp species: $nSpSpp\n";
print $SRYOUT  "Number of genera: " . scalar( keys %genus ) . "\n";
print $SRYOUT  "Number of genera with only one species and the species being _sp: $nSpGenera\n";

print $SRYOUT  "\n_sp genera\n";
for (@spGenera)
{
  print $SRYOUT  "$_\n";
}

print $SRYOUT "\nNumber of singletons (species with only one reference sequence): $sppFreq2{1}\n\n";


print "--- Generating species ann and query files\n";
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

## Rerunning vicut with updated taxons
## Hoping that sequences whose taxonomy changed to <genus>_sp will have their
## taxonomy updated to a proper species taxonomy, as now they will be put in the
## new query file
my $vicutFinalTx = "$vicutDir/updated.tx";
my %newSpp = readTbl($vicutFinalTx);

my $newTxFile = $vicutFinalTx;

## Putting <genus>_sp seq's to query file.

## In the first run of vicut <geuns>_sp, also referred to as _sp species, were
## treated as sequences with known annotation.

## Here we check if they form more than one cluster and if they do, only
## sequences from the largest cluster of the parent genus are treated as having
## 'known taxonomy'. The rest are put in the query file.

## Where do the "Unclassified" come from ?

## changing !defined $sspSpecies{$t}
## to
## push @query2, "$id\n" if !defined $sspSeqID{$id};

## where $sspSeqID{$id} is defined for seq's from largest clusters of the parent
## genus of a singleton _sp.

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

my %sspSeqID; # seqID => 1 if seqID is in the largest cluster of a singleton
	      # species; How to figure out if it is within the genus of that
	      # species?

# Here is an example of a singleton species Iamia_sp before and after taxonomy_cleanup

# BEFORE taxonomy_cleanup

# o_Acidimicrobiales:     1678
# 	f_Iamiaceae:                       450
# 		g_Iamia: 450
# 			Iamia_sp: 450

# AFTER taxonomy_cleanup

# o_Acidimicrobiales:     1678
# 	f_Iamiaceae:                       450
# 		g_Iamia: 450
# 			Aciditerrimonas_sp_tGT_1:    17
# 			Iamia_sp_tGT_1:             431
# 			Iamia_sp_tGT_2:               1
# 	f_Acidimicrobineae_incertae_sedis: 566
# 		g_Aciditerrimonas: 566
# 			Aciditerrimonas_sp_tGT_1:   503
# 			Aciditerrimonas_sp_tGT_2:    11
# 			Aciditerrimonas_sp_tGT_3:     2
# 			Aciditerrimonas_sp_tGT_4:     1
# 			Aciditerrimonas_sp_tGT_5:     1
# 			Iamia_sp_tGT_1:              48
# 	f_Acidimicrobiaceae:               662
# 		g_Ferrithrix:       6
# 			Aciditerrimonas_sp_tGT_1:   6
# 		g_Acidimicrobium:  14
# 			Aciditerrimonas_sp_tGT_1:  14
# 		g_Ferrimicrobium:  31
# 			Aciditerrimonas_sp_tGT_1:  31
# 		g_Ilumatobacter:  611
# 			Aciditerrimonas_sp_tGT_4:     1
# 			Iamia_sp_tGT_1:             594
# 			Ilumatobacter_fluminis:       1
# 			Ilumatobacter_sp:             3
# 			Ilumatobacter_sp_tGT_13:      9

for my $sp (keys %spFreqTbl)
{
  if (defined $sspSpecies{$sp})
  {
    my @cls = sort { $spFreqTbl{$sp}{$b} <=> $spFreqTbl{$sp}{$a} } keys %{$spFreqTbl{$sp}};
    if ( @cls > 1 && $spFreqTbl{$sp}{$cls[0]} > $spFreqTbl{$sp}{$cls[1]} )
    {
      print "\n\nProcessing $sp\tnCltrs: " . @cls . "\n";
      print "Cluster sizes: ";
      map { print $spFreqTbl{$sp}{$_} . ", "} @cls;

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
      print "Query2: $id\t$t\n" if $debug;
    }
    else
    {
      push @ann2, "$id\t$t\n";
    }
  }
  elsif ( $g eq "Unclassified" )
  {
    push @query2, "$id\n";
    print "Query2b: $id\t$t\n" if $debug;
  }
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
}

print "--- Running genotype_spp.pl\n";
$cmd = "genotype_spp.pl -d $vicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating final lineage and taxonomy files\n";

##my $vicutFinalTx = "$vicutDir/minNodeCut.cltrs";

$newTxFile = "$vicutDir/updated2.tx";
my %newTx = readTbl($newTxFile);

my $newTxFileNoTGTs = "$vicutDir/updated.tx";
my %newTxNoTGTs = readTbl($newTxFileNoTGTs);


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

my $finalLineageFile = $grPrefix . "_final.lineage";
my $finalLineageFile2 = $grPrefix . "_final_no_tGTs.lineage";
open LOUT, ">$finalLineageFile" or die "Cannot open $finalLineageFile for writing: $OS_ERROR\n";
open LOUT2, ">$finalLineageFile2" or die "Cannot open $finalLineageFile2 for writing: $OS_ERROR\n";
for my $id (keys %newTx)
{
  my $lineage = $lineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  print LOUT "$id\t";
  print LOUT2 "$id\t";
  for (@f)
  {
    print LOUT "$_;";
    print LOUT2 "$_;";
  }
  if ( ! exists $ogInd{$id} )
  {
    print LOUT $newTx{$id} . "\n";
    print LOUT2 $newTxNoTGTs{$id} . "\n";
  }
  else
  {
    print LOUT $sp . "\n";
    print LOUT2 $sp . "\n";
  }
}
close LOUT;
close LOUT2;

my $finalTxFile = $grPrefix . "_final.tx";
if ($rmOGfromTx) # if OG seq's are to be removed from the final tx file
{
  print "   Removing outgroup sequences from final taxonomy file, as requested";
  open OUT, ">$finalTxFile" or die "Cannot open $finalTxFile for writing: $OS_ERROR\n";
  for my $id (keys %newTx)
  {
    if ( ! exists $ogInd{$id} )
    {
      print OUT "$id\t" . $newTx{$id} . "\n";
    }
    else
    {
      print "Outgroup sequence $id (" . $newTx{$id} . ") excluded from $finalTxFile\n";
    }
  }
  close OUT;
}
else
{
  my $ap = abs_path( $newTxFile );
  my $cmd = "rm -f $finalTxFile; ln -s $ap $finalTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

%newTx = readTbl($finalTxFile);

## final alignment
my $finalAlgnFile = $grPrefix . "_algn_trimmed_final.fa";

my $finalnSeqs  = keys %newTx;
my $nSeqs       = keys %spp;

if ($nSeqs != $finalnSeqs)
{
  print "\nInitial number of sequences: $nSeqs\n";
  print "Final number of sequences: $finalnSeqs\n\n";

  print "--- Generating final alignment files (removing sequences lost in taxonomy cleanup)\n";
  $cmd = "cut -f1 $finalTxFile > final.seqIDs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "select_seqs.pl --quiet -s final.seqIDs -i $trimmedAlgnFile -o $finalAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
else
{
  my $ap = abs_path( $trimmedAlgnFile );
  $cmd = "rm -f $finalAlgnFile; ln -s $ap $finalAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## comparison between old and new species assignments
print "--- Comparing old and new species assignments\n";
$cmd = "cmp_tx.pl -i $txFile -j $newTxFile -o old_vs_new_spp_cmp";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "\nTaxonomy AFTER cleanup\n";
printLineageAfterSpp();

printLineageAfterSppToFile();


##
## phylo tree with final taxonomy
##

## extracting leave IDs
$treeLeavesFile = "$grPrefix" . "_tree.leaves";
$cmd = "nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## looking at the difference between leaf IDs and newTx keys
@treeLeaves = readArray($treeLeavesFile);
my @survivedIDs = keys %newTx;
my @lostLeaves = diff(\@treeLeaves, \@survivedIDs);

## prunning
my $prunedTreeFile = $grPrefix . "_tx_pruned.tree";
if (@lostLeaves>0)
{
  print "\n\tSpecies cleanup eliminated " . @lostLeaves . " sequences\n";
  print "--- Pruning lost seqIDs from the current phylo tree\n";
  $cmd = "nw_prune $treeFile @lostLeaves > $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}
else
{
  $prunedTreeFile = $treeFile;
}


print "\n--- Generating a tree with final species names at leaves\n";
my $finalSppTreeFile = "$grPrefix" . "_final_spp.tree";
$cmd = "rm -f $finalSppTreeFile; nw_rename $prunedTreeFile $newTxFile | nw_order - > $finalSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Generating tree with <final species name>_<seqID> labels at leaves\n";
my $finalSppSeqIDsFile = "$grPrefix" . "_final_spp.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalTxFile > $finalSppSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalSppSeqIdTreeFile = "$grPrefix" . "_final_sppSeqIDs.tree";
$cmd = "rm -f $finalSppSeqIdTreeFile; nw_rename $prunedTreeFile $finalSppSeqIDsFile | nw_order -  > $finalSppSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final species clades collapsed to a single node \n";
my $finalCondSppTreeFile = "$grPrefix" . "_final_sppCondensed.tree";
$cmd = "rm -f $finalCondSppTreeFile; nw_condense $finalSppTreeFile > $finalCondSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


##
## genus-level cleanup
##
print "--- Generating genus ann and query files\n";

my $geFile = "$grPrefix" . ".genus";
print "--- Writing genus assignments to $geFile file\n";
writeTbl(\%gen, $geFile);


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
my $geVicutDir = "genus_vicut_dir";
$cmd = "vicut -t $prunedTreeFile -a $annFile -o $geVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $genVicutTxFile = "$vicutDir/minNodeCut.cltrs";
if ( ! -f $genVicutTxFile )
{
  warn "\nERROR: $genVicutTxFile does not exist";
  exit;
}

my %id2genCltrTbl;  # seqID => <genus cltr ID>
#my %cltrTbl;     # cltrID => ref to array of IDs in the given cluster
open IN, "$genVicutTxFile" or die "Cannot open $genVicutTxFile for reading: $OS_ERROR\n";
$header = <IN>;
foreach (<IN>)
{
  chomp;
  my ($id, $cl, $tx) = split /\s+/,$_;
  $id2genCltrTbl{$id} = $cl;
  #push @{$cltrTbl{$cl}}, $id;
}
close IN;

my $geprunedTreeFile = $grPrefix . "_genus.tree";
$cmd = "cp $prunedTreeFile $geprunedTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running update_tx.pl for genus level\n";
$cmd = "update_tx.pl $debugStr -a $finalTxFile -d $geVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Running genotype_genus.pl for genus level\n";
$cmd = "genotype_genus.pl -d $geVicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalgeTx = "$geVicutDir/updated2.tx";

print "\n--- Generating a tree with final genus names at leaves\n";
my $finalgeTreeFile = "$grPrefix" . "_final_genus.tree";
$cmd = "rm -f $finalgeTreeFile; nw_rename $geprunedTreeFile $finalgeTx | nw_order - > $finalgeTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating tree with <final genus name>_<seqID> labels at leaves\n";
my $finalgeSeqIDsFile = "$grPrefix" . "_final_genus.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $finalgeTx > $finalgeSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalgeSeqIdTreeFile = "$grPrefix" . "_final_geSeqIDs.tree";
$cmd = "rm -f $finalgeSeqIdTreeFile; nw_rename $geprunedTreeFile $finalgeSeqIDsFile | nw_order -  > $finalgeSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final genera collapsed to a single node \n";
my $finalCondgeTreeFile = "$grPrefix" . "_final_geCondensed.tree";
$cmd = "rm -f $finalCondgeTreeFile; nw_condense $finalgeTreeFile > $finalCondgeTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;




## ---------------------------------------
##    Taxonomy Cleanup Summary Stats
## ---------------------------------------

my $txBFile = "$vicutDir/updated.tx";
my %sppB = read2colTbl($txBFile);

my @txVals = values %spp;
my @uqTx = unique(\@txVals);
# my @uqTxAll = unique(\@txVals);
# my @uqTx;
# for (@uqTxAll)
# {
#   if (!/OUTGROUP.*/) # this is redundant given that outgroup seq's are filtered out in readTxTbl()
#   {
#     push @uqTx, $_;
#   }
# }

## debugging
# printArray(\@uqTx, "uqTx");
# exit;


my @txBVals = values %sppB;
my @uqTxBAll = unique(\@txBVals);

my @uqTxB;
for (@uqTxBAll)
{
  if (!/OUTGROUP/)
  {
    push @uqTxB, $_;
  }
}

my %sppTrB; ## species transition table
for ( keys %spp )
{
  if ( exists $sppB{$_} )
  {
    $sppTrB{ $spp{$_} }{ $sppB{$_} }++;
  }
  else
  {
    $sppTrB{ $spp{$_} }{ "NA" }++;
  }
}

## array of lost taxons
my @diffTxB = diff(\@uqTx, \@uqTxB);

@diffTxB = sort { $sppFreq{$b} <=> $sppFreq{$a} } @diffTxB;


my %sppFreqB; ## table of number of sequences per species
map { $sppFreqB{$_}++ } values %sppB;

my %sppFreqB2; ## frequency table of species sequence frequencies
map { $sppFreqB2{$_}++ } values %sppFreqB;

## number of _sp species
my $nSpSppB = 0;
my %genusB;
for (keys %sppFreqB)
{
  my ($g, $s) = split "_";
  $nSpSppB++ if $s && $s eq "sp";
  push @{$genusB{$g}}, $_;
}

## Number of genera with only one species and the species being _sp
my $nSpGeneraB = 0;
my @spGeneraB;
my %sspSpeciesB;
for my $g (keys %genusB)
{
  if ( scalar(@{$genusB{$g}})==1 )
  {
    my $sp = $genusB{$g}->[0];
    my ($g, $s) = split "_", $sp;
    if ($s && $s eq "sp")
    {
      $nSpGeneraB++;
      push @spGeneraB, $g;
      $sspSpeciesB{$sp} = 1;
    }
  }
}


print $SRYOUT "\n--- After vicut majority vote taxonomy changes (before genotyping and size filtering of taxon table)\n";

print $SRYOUT "\n\nNumber of species: " . (scalar( keys %sppFreqB )-1) . "\n"; # subtracting OUTGROUP

print $SRYOUT "Number of lost species: " . scalar( @diffTxB ) . "\n";
# print $SRYOUT "\nLost species:\n";
# map{ print $SRYOUT $_ . "\n"} @diffTxB;

print $SRYOUT "\nLost species taxonomy changes:\n";
my $maxStrLenB = 0;
map { $maxStrLenB = length($_) if( length($_) > $maxStrLenB )} @diffTxB;

for my $sp ( @diffTxB )
{
  my @newSpp = keys %{$sppTrB{$sp}};
  my @foo = grep{ !/^$sp$/ } @newSpp;

  # print "\n\n$sp: newSpp: @newSpp\n";
  # print "foo: @foo\n";

  my $n = $maxStrLenB - length($sp);
  my $pad = "";
  for (my $i=0; $i<$n; $i++)
  {
    $pad .= " ";
  }

  if ( @foo > 0 )
  {
    ##print "\n$sp: @newSpp\n";
    print $SRYOUT "$sp:$pad";
    for (@newSpp)
    {
      print $SRYOUT  "  " . $_ . "(" . $sppTrB{$sp}->{$_} . ")";
    }
    print $SRYOUT "\n";
  }
}


print $SRYOUT "\nNumber of _sp species: $nSpSppB\n";
print $SRYOUT "Number of genera: " . scalar( keys %genusB ) . "\n";
print $SRYOUT "Number of genera with only one species and the species being _sp: $nSpGeneraB\n";

print $SRYOUT "\n_sp genera\n";
for (@spGeneraB)
{
  print $SRYOUT "\t$_\n";
}

my @lostSpGeneraB = diff(\@spGenera, \@spGeneraB);

print $SRYOUT  "\nLost _sp genera:\n";
for (@lostSpGeneraB)
{
  print $SRYOUT  "\t" . $_ . "\n";
}
print $SRYOUT  "\n";


print $SRYOUT "Number of singletons (species with only one reference sequence): $sppFreqB2{1}\n\n";



## reporting taxonomic changes
## 1. which taxons has been lost
## 2. log file of old and new taxonomic assignments for each sequence
## 3. Enumerate transition types <old tx> => <new tx> and their frequencies


print $SRYOUT "\n--- After genotyping and size filtering of _sp species clusters\n";

## changes that took place after genotyping and size filtering of _sp species clusters

##%newTx = readTbl($newTxFile);
my @newTxVals = values %newTx;
my @uqNewTxAll = unique(\@newTxVals);

my @uqNewTx;
for (@uqNewTxAll)
{
  if (!/OUTGROUP/)
  {
    push @uqNewTx, $_;
  }
}

## array of lost taxons
my @diffTx = diff(\@uqTx, \@uqNewTx);

@diffTx = sort { $sppFreq{$b} <=> $sppFreq{$a} } @diffTx;

print "\n--- Generating list of lost species and with the number of sequences that represented it";
my $diffFile = "lost_spp.txt";
open OUT, ">$diffFile" or die "Cannot open $diffFile for writing: $OS_ERROR\n";
##for ( sort { $sppFreq{$b} <=> $sppFreq{$a} } @diffTx)
for ( @diffTx )
{
  print OUT "$_\t$sppFreq{$_}\n";
  ##print "\t$_\t$sppFreq{$_}\n" if $debug;
}
close OUT;


## 2. log file of old and new taxonomic assignments for each sequence
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

print "\n\tLog file of old and new taxonomic assignments for each sequence written to $grDir/$oldNewTxFile\n";

if ($verbose)
{
  print "\n";
  printArray(\@uqTx, "\n--- Species BEFORE filtering");

  # print "\n";
  # printArray(\@uqNewTx, "\n--- Species AFTER filtering");
  # print "\n";

  print "\n--- Lost species\n";
  printFormatedTbl(\%sppFreq, \@diffTx);
  print "\n";

  print "\n--- Species AFTER filtering";

  my %sppFreq2; ## table of number of sequences per species
  map { $sppFreq2{$_}++ } @newTxVals;
  for ( sort { $sppFreq2{$b} <=> $sppFreq2{$a} } keys %sppFreq2 )
  {
    print "$_\t$sppFreq2{$_}\n";
  }
  print "\n";

  # printArray(\@diffTx, "diffTx");
}

## 3. Enumerate transition types <old tx> => <new tx> and their frequencies
## only species for which there was a change of assignment are reported
##print "\n--- Enumerating taxon transitions <old tx> => <new tx>\n";
my $sppTrTypes = "spp_transitions.txt";
open OUT, ">$sppTrTypes" or die "Cannot open $sppTrTypes for writing: $OS_ERROR\n";
##for my $sp ( keys %sppTr )

my $maxStrLen = 0;
map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @diffTx;

for my $sp ( @diffTx )
{
  ##next if $sp =~ /^c\d+/;

  my @newSpp = keys %{$sppTr{$sp}};
  my @foo = grep{ !/^$sp$/ } @newSpp;

  # print "\n\n$sp: newSpp: @newSpp\n";
  # print "foo: @foo\n";

  my $n = $maxStrLen - length($sp);
  my $pad = "";
  for (my $i=0; $i<$n; $i++)
  {
    $pad .= " ";
  }

  if ( @foo > 0 )
  {
    ##print "\n$sp: @newSpp\n";
    print OUT "$sp:$pad";
    print "$sp:$pad" if $verbose;
    for (@newSpp)
    {
      print OUT "  " . $_ . "(" . $sppTr{$sp}->{$_} . ")";
      print  "  " . $_ . "(" . $sppTr{$sp}->{$_} . ")" if $verbose;
    }
    print OUT "\n";
    print "\n" if $verbose;
  }
}
close OUT;

print "\n\tTaxon transitions <old tx> => <new tx> written to $grDir/$sppTrTypes\n";


my $readmeFile = "README";
generateREADMEfile($readmeFile);

print "\n\tSee $grDir/$readmeFile for info about the algorithm and the content of the output directory.\n";

print $SRYOUT "\n\n--- Final summary\n";

print $SRYOUT  "\n\tNumber of species before spp cleanup: " . scalar( keys %sppFreq ) . "\n";
print $SRYOUT    "\tNumber of species after spp cleanup:  " . @uqNewTx . "\n";
print $SRYOUT    "\tNumber of lost species: " . @diffTx . "\n\n";


## number of _sp after spp cleanup species
my $nSpSppFinal = 0;
my %genusFinal;
for (@uqNewTx)
{
  my ($g, $s) = split "_";
  if (defined $s)
  {
    $nSpSppFinal++ if $s eq "sp";
  }
  else
  {
    warn "WARNING in genusFinal loop: genus/species split of $_ was not performed";
  }

  push @{$genusFinal{$g}}, $_;
}

print $SRYOUT  "\tNumber of _sp species before spp cleanup: $nSpSpp\n";
print $SRYOUT  "\tNumber of _sp species after spp cleanup: $nSpSppFinal\n";
print $SRYOUT  "\tNumber of lost _sp species: " . ($nSpSpp - $nSpSppFinal) . "\n\n";

my %sppFreqFinal; ## table of number of sequences per species
map { $sppFreqFinal{$_}++ } values %newTx;

my %sppFreq2Final; ## frequency table of species sequence frequencies
map { $sppFreq2Final{$_}++ } values %sppFreqFinal;

print $SRYOUT  "\tNumber of singletons species before spp cleanup: $sppFreq2{1}\n";
print $SRYOUT  "\tNumber of singletons species after spp cleanup: $sppFreq2Final{1}\n\n";

## Number of genera with only one species and the species being _sp
my $nSpGeneraFinal = 0;
my @spGeneraFinal;
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
    }
  }
}

print $SRYOUT  "\tNumber of genera before spp cleanup: " . scalar( keys %genus ) . "\n";
print $SRYOUT  "\tNumber of genera after spp cleanup: " . scalar( keys %genusFinal ) . "\n";
my $nLostGenera = (scalar( keys %genus ) - scalar( keys %genusFinal ));
print $SRYOUT  "\tNumber of lost genera: $nLostGenera\n\n";

## Lost genera
if ($nLostGenera)
{
  my @beforeGenera = keys %genus;
  my @afterGenera = keys %genusFinal;
  my @lostGenera = diff(\@beforeGenera, \@afterGenera);

  print $SRYOUT  "\tLost genera:\n";
  for (@lostGenera)
  {
    print $SRYOUT  "\t\t" . $_ . "\n";
  }
  print $SRYOUT  "\n";
}

print $SRYOUT "\n\tFinal genera and their species were written to $grDir/finalGenera.txt\n\n";
my $gFile = "finalGenera.txt";
open OUT, ">$gFile" or die "Cannot open $gFile for writing: $OS_ERROR\n";
for my $g ( sort { scalar(@{$genusFinal{$b}}) <=> scalar(@{$genusFinal{$a}}) } keys %genusFinal)
{
  print OUT $g . "\n";
  for (@{$genusFinal{$g}})
  {
    print OUT"\t$_\n";
  }
}
close OUT;


print $SRYOUT  "\tNumber of genera with only one species and the species being _sp before spp cleanup: $nSpGenera\n";
print $SRYOUT  "\tNumber of genera with only one species and the species being _sp after spp cleanup: $nSpGeneraFinal\n";
my $nLostSpGenera = $nSpGenera - $nSpGeneraFinal;
print $SRYOUT  "\tNumber of lost genera with only _sp species: $nLostSpGenera\n\n";

## Lost _sp genera
if ($nLostSpGenera)
{
  my @beforeSpGenera = @spGenera;
  my @afterSpGenera = @spGeneraFinal;
  my @lostSpGenera = diff(\@beforeSpGenera, \@afterSpGenera);

  print $SRYOUT  "\tLost _sp genera:\n";
  for (@lostSpGenera)
  {
    print $SRYOUT  "\t\t" . $_ . "\n";
  }
  print $SRYOUT  "\n";
}

close $SRYOUT;

print "\n\tSummary stats written to $grDir/$summaryStatsFile\n\n";

# print "\tFinal alignment file: $grDir/$finalAlgnFile\n";
# print "\tNOTE: The final alignment and taxonomy files contain an outgroup sequence.\n";
print "\tFinal taxonomy file: $grDir/$finalTxFile\n\n";
print "\tNOTE: The final taxonomy file contain outgroup sequences.\n";

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
  my $fh = shift;

  my @phs = keys %phTbl;
  my $phKVL = getKeyValStrLengths(\%phTbl);

  print $fh "\nTaxonomy before cleanup\n";
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
