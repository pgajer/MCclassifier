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

=item B<--taxon-size-thld>
  Upper limit for the number of elements within each taxon

=item B<--build-model-data>
  Generate
    - species lineage file (grPrefix_final.spLineage)
    - taxon file (grPrefix_final.tx)
    - ungapped fasta file corresponding to sequences present in the taxon file

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

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Firmicutes_dir

  taxonomy_cleanup.pl --debug -i Firmicutes_group_6

  taxonomy_cleanup.pl --do-not-pop-pdfs --build-model-data --use-long-spp-names -i Firmicutes_group_6_V3V4

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

my $maxStrLen     = 150;  # maximal number of characters in a taxon string; if
			  # larger file names using the taxon name will use an
			  # abbreviated version of the taxon name
my $taxonSizeThld = 20;

GetOptions(
  "input-group|i=s" 	=> \my $grPrefix,
  "rm-OGs-from-tx|r"    => \my $rmOGfromTx,
  "skip-FastTree"       => \my $skipFastTree,
  "use-long-spp-names"  => \my $useLongSppNames,
  "taxon-size-thld"     => \$taxonSizeThld,
  "show-all-trees"      => \my $showAllTrees,
  "do-not-pop-pdfs"     => \my $doNotPopPDFs,
  "build-model-data"    => \my $buildModelData,
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
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

my $mothur = "/Users/pgajer/bin/mothur";
my $readNewickFile = "/Users/pgajer/.Rlocal/read.newick.R";

my $igsStr = "";
if ( defined $igs )
{
  $mothur = "/usr/local/packages/mothur-1.36.1/mothur";
  $igsStr = "--igs";
  $readNewickFile = "??";
}

my $johannaStr = "";
if ( defined $johanna )
{
  $mothur = "/Users/pgajer/bin/mothur";
  $johannaStr = "--johanna";
  $readNewickFile = "/Users/jholm/MCclassifier/perl/read.newick.R";
}

## Export LD_LIBRARY_PATH=/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64

local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/readline/lib:/usr/local/packages/gcc-5.3.0/lib64";

my $useLongSppNamesStr = "";
if ( defined $useLongSppNames )
{
  $useLongSppNamesStr = "--use-long-spp-names";
}

my $showAllTreesStr = "";
if ($showAllTrees)
{
  $showAllTreesStr = "--show-tree";
}

####################################################################
##                               MAIN
####################################################################

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

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  exit;
}

my $section = qq~

##
## Species-level cleanup
##

~;
print "$section";

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

print "--- Testing if OG seq's form a monophylectic clade at the top or bottom of the input tree\n";
test_OG($treeFile, \%ogInd);

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
@ogSeqIDs = comm( \@ogSeqIDs, \@commSeqIDs );

if ( scalar(@ogDiff) != 0 )
{
  warn "\n\tWARNING the following outgroup seq IDs are not in the trimmed alignment file:\n\n";
  printArray(\@ogDiff);
}

if ( scalar(@ogSeqIDs) == 0 )
{
  warn "\n\tERROR: All outgroup seq's were lost";
  print "\n\n";
  exit;
}

if ($debug)
{
  print "\n\tNumber of seq's in the trimmed alignment/lineage files: " . @seqIDs . "\n";
  print "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n\n";
}

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
  print     "\tNumber of seq IDs of $trimmedAlgnFile: " . @seqIDs . "\n";
  print     "\tNumber of common seq IDs elements: " . @commTL . "\n\n";
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
my %geLineage; # genus   => lineage of the genus (skipping species)
my %spParent;

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

  $spParent{$sp} = $ge;

  $spLineage{$sp} = $lineage;

  @f = split ";", $lineage;
  $sp = pop @f;
  $lineage = join ";", @f;
  $geLineage{$ge} = $lineage;

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

if ($debug)
{
  print "\nTaxonomy BEFORE cleanup\n";
  printLineage();
}

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
  $cmd = "vicut $quietStr -t $treeFile -a $annFile -q $queryFile -o $vicutDir";
}
else
{
  $cmd = "vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
}
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## update_spp_tx.pl takes vicut clusters and updates taxonomy using majority vote.
## Moreover, if a species is present in more than 1 cluster, it changes
## taxonomies of sequences from the small clusters to <genus>_sp
print "--- Running update_spp_tx.pl\n";
$cmd = "update_spp_tx.pl $quietStr $debugStr $useLongSppNamesStr -a $txFile -d $vicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## Rerunning vicut with updated taxons giving a chance sequences whose taxonomy
## changed to <genus>_sp to have their taxonomy updated to a proper species
## taxonomy, as now they will be put in the new query file.

my $updatedTxFile = "$vicutDir/updated.tx";
my %newTx = readTbl($updatedTxFile);

my %idCl;    # seqID => cltrID
my %cltrTbl; # cltrID => ref to array of IDs in the given cluster
my $cltrFile = "$vicutDir/minNodeCut.cltrs";
open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR\n";
my $header = <IN>;
foreach (<IN>)
{
  chomp;
  my ($id, $cl, $tx) = split /\s+/,$_;
  $idCl{$id} = $cl;
  push @{$cltrTbl{$cl}}, $id;
}
close IN;

my %spFreqTbl;
for my $id ( keys %newTx )
{
  my $sp = $newTx{$id};
  my $cl = $idCl{$id};
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
for my $id ( keys %newTx )
{
  my $t = $newTx{$id};
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
  $cmd = "vicut $quietStr -t $treeFile -a $annFile2 -q $queryFile2 -o $vicutDir";
  ##$cmd = "vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Running update_spp_tx.pl\n";
  $cmd = "update_spp_tx.pl $quietStr $debugStr $useLongSppNamesStr -a $updatedTxFile -d $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $updatedTxFile = "$vicutDir/updated.tx";
  %newTx = readTbl($updatedTxFile);
}

print "--- Updating lineageTbl after second vicut and detecting seq's that had their taxonomy changed to OG\n";
my @extraOG; # array of seqIDs that had their taxonomy changed to OG - they will be removed
for my $id (keys %lineageTbl)
{
  next if exists $ogInd{$id};

  if ( exists $newTx{$id} )
  {
    my $newSp = $newTx{$id};

    if ( $newSp =~ /OG$/ )
    {
      print "\n\tWARNING: non-OG seq $id was changed to OG species => $newSp; Scheduling removal of the sequence.\n";
      print "\tlineageTbl{$id}: " . $lineageTbl{$id} . "\n\n";
      push @extraOG, $id;
      delete $lineageTbl{$id};
      delete $newTx{$id};
      next;
    }

    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;

    if ($newSp ne $sp)
    {
      my ($g, $s) = split "_", $newSp;
      if ( exists $geLineage{$g} )
      {
	$lineageTbl{$id} = $geLineage{$g} . ";$newSp";
      }
      else
      {
	if (exists $spParent{$newSp})
	{
	  $g = $spParent{$newSp};
	  $lineageTbl{$id} = $geLineage{$g} . ";$newSp";
	}
	else
	{
	  warn "\n\n\tERROR: $g not found in geLineage";
	  print "\tand spParent{$newSp} not found\n";
	  print "\tnewSp: $newSp\n";
	  print "\tlineageTbl{$id}: " . $lineageTbl{$id} . "\n";
	  print "\n\n";
	  exit;
	}
      }
    }
  } # end of if ( exists $newTx{$id} )
}


if ( @extraOG>0 )
{
  delete @lineageTbl{@extraOG};
  delete @newTx{@extraOG};

  # removing elements of extraOG from the updated.tx file
  $updatedTxFile = "$vicutDir/updated.tx";
  print "--- Removing elements of extraOG from $updatedTxFile\n";
  my $origNewTxFileNoTGTs = "$vicutDir/updated_orig.tx";
  my $extraOGfile = "$vicutDir/extraOG.seqIDs";
  writeArray(\@extraOG, $extraOGfile);

  $cmd = "mv $updatedTxFile $origNewTxFileNoTGTs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "select_tx.pl $quietStr -e $extraOGfile -i $origNewTxFileNoTGTs -o $updatedTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # pruning alignment
  print "--- Pruning extraOG seq's from $trimmedAlgnFile\n";
  my $prunedAlgnFile = $grPrefix . "_algn_trimmed_pruned.fa";
  $cmd = "select_seqs.pl $quietStr -e $extraOGfile -i $trimmedAlgnFile -o $prunedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $trimmedAlgnFile = $prunedAlgnFile;

  # pruning tree
  print "--- Pruning extraOG seq's from $treeFile\n";
  my $prunedTreeFile = $grPrefix . "_pruned.tree";
  $cmd = "rm -f $prunedTreeFile; nw_prune $treeFile @extraOG > $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $treeFile = $prunedTreeFile;
}

my $sppLineageFile = $grPrefix . "_spp.lineage";
open OUT, ">$sppLineageFile" or die "Cannot open $sppLineageFile for writing: $OS_ERROR\n";
my $wCount = 1;
for my $id (keys %newTx)
{
  if ( exists $lineageTbl{$id} )
  {
    print OUT "$id\t" . $lineageTbl{$id} . "\n";
  }
  elsif ( exists $ogLineageTbl{$id} )
  {
    print OUT "$id\t" . $ogLineageTbl{$id} . "\n";
  }
  else
  {
    warn "$wCount: WARNING: $id does not exist in lineageTbl and ogLineageTbl";
    $wCount++;
  }
}
close OUT;

## Creating symbolic link to the most recent trimmed alignment file for use with
## truncate_algn.pl
my $finalAlgnFile = $grPrefix . "_algn_trimmed_final.fa";
my $ap = abs_path( $trimmedAlgnFile );
$cmd = "rm -f $finalAlgnFile; ln -s $ap $finalAlgnFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


## Removing outgroup sequences from the updated.tx file and the phylogenetic
## tree.

## Deleting outgroup sequences from %newTx !!!!
delete @newTx{@ogSeqIDs};
delete @lineageTbl{@ogSeqIDs};

my %newTxNoTGTs = %newTx; ## readTbl($updatedTxFile);

# removing OG seq's from the updated.tx file
#$updatedTxFile = "$vicutDir/updated.tx";
print "--- Removing OG seq's from $updatedTxFile\n";
my $origNewTxFileNoTGTs = "$vicutDir/updated_orig2.tx";
my $ogFile = "$vicutDir/og.seqIDs";
writeArray(\@ogSeqIDs, $ogFile);

$cmd = "mv $updatedTxFile $origNewTxFileNoTGTs";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "select_tx.pl $quietStr -e $ogFile -i $origNewTxFileNoTGTs -o $updatedTxFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# pruning tree from OG seq's
print "--- Pruning OG seq's from $treeFile\n";
my $prunedTreeFile = $grPrefix . "_no_OG_seqs.tree";
$cmd = "rm -f $prunedTreeFile; nw_prune $treeFile @ogSeqIDs > $prunedTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$treeFile = $prunedTreeFile;



## Testing consistency between the phylo tree and sequences after taxonomy
## cleanup. In particular, the set of leave seqIDs has to be equal to the set of
## seqIDs of the new taxonomy.

## extracting leave IDs
$cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

%newTx = readTbl($updatedTxFile);

## looking at the difference between leaf IDs and newTx keys
@treeLeaves = readArray($treeLeavesFile);
my @survivedIDs = keys %newTx;
my @lostLeaves = diff(\@treeLeaves, \@survivedIDs);

## prunning tree from lost IDs
if (@lostLeaves>0)
{
  $prunedTreeFile = $grPrefix . "_no_OG_seqs_pruned.tree";
  print "\n\tSpecies cleanup eliminated " . @lostLeaves . " sequences\n";
  print "--- Pruning lost seqIDs from the current phylo tree\n";
  $cmd = "rm -f $prunedTreeFile; nw_prune $treeFile @lostLeaves > $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $treeFile = $prunedTreeFile;
}

print "--- Removing _sp singletons clusters\n";
%idCl = ();
my %clSize;
my %clSeqs;
$cltrFile = "$vicutDir/minNodeCut.cltrs";
open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR\n";
$header = <IN>;
foreach (<IN>)
{
  chomp;
  my ($id, $cl, $tx) = split /\s+/,$_;
  if (exists $newTx{$id})
  {
    $idCl{$id} = $cl;
    $clSize{$cl}++;
    push @{$clSeqs{$cl}}, $id;
  }
}
close IN;


%spFreqTbl = ();
for my $id ( keys %newTx )
{
  my $sp = $newTx{$id};
  my $cl = $idCl{$id};
  $spFreqTbl{$sp}{$cl}++;
}

my @spSingletons;
for my $sp (keys %spFreqTbl)
{
  my @f = split "_", $sp;
  my $g = shift @f;
  my $s = shift @f;

  if ($s eq "sp")
  {
    my @cls = keys %{$spFreqTbl{$sp}};
    for my $cl (@cls)
    {
      if ($clSize{$cl}==1)
      {
	push @spSingletons, @{$clSeqs{$cl}};
      }
    }
  }
}

if ( @spSingletons )
{
  delete @newTx{@spSingletons};
  delete @lineageTbl{@spSingletons};

  # removing elements of extraOG from the updated.tx file
  my $updatedTxFile = "$vicutDir/updated.tx";
  print "--- Removing elements of spSingletons from $updatedTxFile\n";
  my $origNewTxFileNoTGTs = "$vicutDir/updated_orig2.tx";
  my $spSingletonsFile = "$vicutDir/spSingletons.seqIDs";
  writeArray(\@spSingletons, $spSingletonsFile);

  $cmd = "mv $updatedTxFile $origNewTxFileNoTGTs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "select_tx.pl $quietStr -e $spSingletonsFile -i $origNewTxFileNoTGTs -o $updatedTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # pruning alignment
  print "--- Pruning spSingletons seq's from $trimmedAlgnFile\n";
  my $prunedAlgnFile = $grPrefix . "_algn_trimmed_pruned2.fa";
  $cmd = "select_seqs.pl $quietStr -e $spSingletonsFile -i $trimmedAlgnFile -o $prunedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $trimmedAlgnFile = $prunedAlgnFile;

  # pruning tree
  print "--- Pruning spSingletons seq's from $treeFile\n";
  my $prunedTreeFile = $grPrefix . "_pruned.tree";
  $cmd = "rm -f $prunedTreeFile; nw_prune $treeFile @spSingletons > $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $treeFile = $prunedTreeFile;

  ## running vicut again
  ## NOTE: I am abandoning the strategy of protecting singleton _sp species (the only species of a genus that is an _sp species)
  my @query3;
  my @ann3;
  for my $id ( keys %newTx )
  {
    my $t = $newTx{$id};
    my @f = split "_", $t;
    my $g = shift @f;
    my $suffix = shift @f;
    my $suffix2 = shift @f;

    if ( defined $suffix && $suffix eq "sp" )
    {
      push @query3, "$id\n";
    }
    else
    {
      push @ann3, "$id\t$t\n";
    }
  }

  my $queryFile3 = "spp_query3.seqIDs";
  my $annFile3   = "spp_ann3.tx";
  ## If _sp sequences, print to file, and set vicutDir to 2nd run directory
  $vicutDir  = "spp_vicut_dir3";
  open QOUT, ">$queryFile3" or die "Cannot open $queryFile3 for writing: $OS_ERROR\n";
  print QOUT @query3;
  close QOUT;

  open ANNOUT, ">$annFile3" or die "Cannot open $annFile3 for writing: $OS_ERROR\n";
  print ANNOUT @ann3;
  close ANNOUT;

  print "--- Running vicut on species data the 3rd time\n";
  if (@query3)
  {
    $cmd = "vicut $quietStr -t $treeFile -a $annFile3 -q $queryFile3 -o $vicutDir";
  }
  else
  {
    $cmd = "vicut $quietStr -t $treeFile -a $annFile3 -o $vicutDir";
  }
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Running update_spp_tx.pl\n";
  $cmd = "update_spp_tx.pl $quietStr $debugStr $useLongSppNamesStr -a $updatedTxFile -d $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $updatedTxFile = "$vicutDir/updated.tx";
  %newTx = readTbl($updatedTxFile);

  print "--- Updating lineageTbl after 3rd run of vicut\n";
  for my $id (keys %lineageTbl)
  {
    if ( exists $newTx{$id} )
    {
      my $lineage = $lineageTbl{$id};
      my @f = split ";", $lineage;
      my $sp = pop @f;
      my $newSp = $newTx{$id};
      if ($newSp ne $sp)
      {
	my ($g, $s) = split "_", $newSp;
	if ( exists $geLineage{$g} )
	{
	  $lineageTbl{$id} = $geLineage{$g} . ";$newSp";
	}
	else
	{
	  if (exists $spParent{$newSp})
	  {
	    $g = $spParent{$newSp};
	    $lineageTbl{$id} = $geLineage{$g} . ";$newSp";
	  }
	  else
	  {
	    warn "\n\n\tERROR: $g not found in geLineage";
	    print "\tand spParent{$newSp} not found\n";
	    print "\tnewSp: $newSp\n";
	    print "\tlineageTbl{$id}: " . $lineageTbl{$id} . "\n";
	    print "\n\n";
	    exit;
	  }
	}
      }
    } # end of if ( exists $newTx{$id} )
  }

}

## Checking the number of species, to make sure its > 1
## if not we end here (after updating lineage table)

my %sppFreqFinal; ## table of number of sequences per species
map { $sppFreqFinal{$_}++ } values %newTx;
my $nSpecies = keys %sppFreqFinal;

if ( $nSpecies == 1 )
{
  print "\n\n\tWARNING:  This group of sequences consists of only one species !!!\n";
  print "\tAbandoning taxonomic cleanup at the higher taxonomic ranks\n\n";

  print "--- Updating lineage file at the species level\n";

  my @ids = keys %newTx;
  my $id = $ids[0];
  my $newSp = $newTx{$id};
  my $newSpNoExt = $newSp;
  $newSpNoExt =~ s/_\d+$//;

  my $lineage = $lineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $ge = pop @f;
  my @t1 = (@f, $ge, $newSpNoExt);
  my @t2 = (@f, $ge, $ge, $newSp); # making genus also a sub-genus for the sake of consistency with other phyla
  my $lineageNoExt = join ';', @t1;
  my $lineageExt = join ';', @t2;

  my $finalLineageFile = $grPrefix . "_final.lineage";
  open OUT1, ">$finalLineageFile" or die "Cannot open $finalLineageFile for writing: $OS_ERROR\n";
  my $finalLineageFile2 = $grPrefix . "_final_no_tGTs.lineage";
  open OUT2, ">$finalLineageFile2" or die "Cannot open $finalLineageFile2 for writing: $OS_ERROR\n";
  for my $id (keys %newTx)
  {
    my $lineage = $lineageTbl{$id};
    print OUT1 "$id\t$lineageExt\n";
    print OUT2 "$id\t$lineageNoExt\n";
  }
  for my $id ( keys %ogLineageTbl )
  {
    my $lineage = $ogLineageTbl{$id};
    print OUT2 "$id\t$lineage\n";
  }
  close OUT2;
  close OUT1;

  print "\n\tUpdated lineage tables written to:\n\t\t$finalLineageFile\n\t\t$finalLineageFile2\n\n";

  exit;
}


## Generating species curation diagnostic plots
## condensed species (using nw_condense2 showing the number of sequences in
## each species) with vicut cluster number

print "--- Generating a tree with species names at leaves\n";
my $sppTreeFile = "$grPrefix" . "_preGenotyping_spp.tree";
$cmd = "rm -f $sppTreeFile; nw_rename $treeFile $updatedTxFile | nw_order -c n  - > $sppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

if (0)
{
  print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
  my $sppSeqIDsFile = "$grPrefix" . "_preGenotyping_spp.seqIDs";
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $updatedTxFile > $sppSeqIDsFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $sppSeqIdTreeFile = "$grPrefix" . "_preGenotyping_sppSeqIDs.tree";
  $cmd = "rm -f $sppSeqIdTreeFile; nw_rename $treeFile $sppSeqIDsFile | nw_order -c n  -  > $sppSeqIdTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

# print "--- Generating a condensed tree with species clades collapsed to a single node \n";
# my $condSppTreeFile = "$grPrefix" . "_preGenotyping_spp_condensed.tree";
# $cmd = "rm -f $condSppTreeFile; nw_condense2 $sppTreeFile > $condSppTreeFile";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# print "--- Extracting leaves of the species condensed tree\n";
# my $csppTreeLeavesFile = "$grPrefix" . "_preGenotyping_spp_condensed_tree.leaves";
# my $cmd = "nw_labels -I $condSppTreeFile > $csppTreeLeavesFile";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# my @csppTreeLeaves = readArray($csppTreeLeavesFile);
# printArray(\@csppTreeLeaves, "\ncsppTreeLeaves\n") if $debug;

$cltrFile = "$vicutDir/minNodeCut.cltrs";
%idCl = ();
%clSize = ();
open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR\n";
$header = <IN>;
foreach (<IN>)
{
  chomp;
  my ($id, $cl, $tx) = split /\s+/,$_;
  if (exists $newTx{$id})
  {
    #$sppCltr{$newTx{$id}} = $cl;
    $idCl{$id} = $cl;
    $clSize{$cl}++;
  }
}
close IN;

my %clIdx; # assigning each cluster of index from 1 to <number of cluster> with
# 1 assigned to the largest one, 2 to second largest etc
my $idx = 1;
for my $cl (sort {$clSize{$b} <=> $clSize{$a}} keys %clSize)
{
  $clIdx{$cl} = $idx;
  $idx++;
}

my %idSppCltr;     # species => sppSizeCltr
my %idSppCltrIdx;  # sppSizeCltr => idx
for my $id (keys %newTx)
{
  my $cl = $idCl{$id};
  $idSppCltr{$id} = $newTx{$id} . "_n" . $clSize{$cl} . "_cl_" . $clIdx{$cl};
  $idSppCltrIdx{$idSppCltr{$id}} = $clIdx{$cl};
}

# print "\nidSppCltrIdx table:\n";
# my @a = sort { $idSppCltrIdx{$a} <=> $idSppCltrIdx{$b} } keys %idSppCltrIdx;
# printFormatedTbl(\%idSppCltrIdx, \@a);
# print "\n\n";

print "--- Creating species_nSize_clID tree\n";
my $spSizeCltrFile = "$grPrefix" . "_preGenotyping.spp_size_cltr";
open OUT, ">$spSizeCltrFile" or die "Cannot open $spSizeCltrFile for writing: $OS_ERROR\n";
for my $id (keys %idSppCltr)
{
  print OUT "$id\t" . $idSppCltr{$id} . "\n";
}
close OUT;

print "--- Generating a tree with species names sizes and cluster index at leaves\n";
my $sppSizeCltrTreeFile = "$grPrefix" . "_preGenotyping_sppSizeCltr.tree";
$cmd = "rm -f $sppSizeCltrTreeFile; nw_rename $treeFile $spSizeCltrFile | nw_order -c n  - > $sppSizeCltrTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $condSppSizeCltrTreeFile = "$grPrefix" . "_preGenotyping_sppSizeCltr_condensed.tree";
$cmd = "rm -f $condSppSizeCltrTreeFile; nw_condense $sppSizeCltrTreeFile > $condSppSizeCltrTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $pdfTreeFile = abs_path( "$grPrefix" . "_preGenotyping_sppSizeCltr_condensed_tree.pdf" );
my $condSppSizeCltrTreeFileAbsPath = abs_path( $condSppSizeCltrTreeFile );
# my $idSppCltrIdxFile = abs_path( "$grPrefix" . "_preGenotyping_sppSizeCltr.idx" );
# writeTbl(\%idSppCltrIdx, $idSppCltrIdxFile);

plot_tree_bw($condSppSizeCltrTreeFileAbsPath, $pdfTreeFile);

if ( $showAllTrees && $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


print "--- Running genotype_spp.pl\n";
## Here species that appear more than twice are being indexed so each species
## forms a monophyletic clade on the reference tree
$cmd = "genotype_spp.pl $debugStr -d $vicutDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $updatedTxFile2 = "$vicutDir/updated2.tx";
%newTx = readTbl($updatedTxFile2);

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

print "--- Updating lineageTbl after genotype_spp.pl\n";
for my $id (keys %lineageTbl)
{
  if ( exists $newTx{$id} )
  {
    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;

    my $newSp = $newTx{$id};
    my ($g, $s) = split "_", $newSp;
    if ($s)
    {
      $spParent{$newSp} = $g;
    }
    elsif (!exists $spParent{$newSp})
    {
      warn "\n\n\tERROR: count not find a genus corresponding to $newSp";
      print "\n\n";
      exit;
    }

    if ($newSp ne $sp)
    {
      if ( exists $geLineage{$g} )
      {
	$lineageTbl{$id} = $geLineage{$g} . ";$newSp";
      }
      else
      {
	warn "\n\n\tERROR: $g not found in geLineage";
	print "\tnewSp: $newSp\n";
	print "\tlineageTbl{$id}: " . $lineageTbl{$id} . "\n";
	print "\n\n";
	exit;
      }
    }
  } # end of if ( exists $newTx{$id} )
  else
  {
    delete $lineageTbl{$id};
  }
}

print "--- Testing again consistency between the phylo tree and sequences after taxonomy cleanup\n";

## extracting leave IDs
$cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## looking at the difference between leaf IDs and newTx keys
@treeLeaves = readArray($treeLeavesFile);
@survivedIDs = keys %newTx;
@lostLeaves = diff(\@treeLeaves, \@survivedIDs);

## prunning tree
if (@lostLeaves>0)
{
  $prunedTreeFile = $grPrefix . "_no_OG_seqs_pruned2.tree";
  print "\n\tSpecies cleanup eliminated " . @lostLeaves . " sequences\n\n" if $debug;
  print "--- Pruning lost seqIDs from the current phylo tree\n";
  $cmd = "rm -f $prunedTreeFile; nw_prune $treeFile @lostLeaves > $prunedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $treeFile = $prunedTreeFile;
}

print "--- Creating a symbolic link to the most recent version of the phylogenetic tree\n";
my $finalTreeFile = $grPrefix . "_final.tree";
$ap = abs_path( $treeFile );
$cmd = "rm -f $finalTreeFile; ln -s $ap $finalTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Testing consistency between lineageTbl and newTx keys\n";
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

## Generating some summary tables
undef %sppFreqFinal; ## table of number of sequences per species
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

# ## Creating a symbolic link for final.tx to point to vicutDir/updated2.tx
# my $finalTxFile = $grPrefix . "_final.tx";
# $ap = abs_path( $updatedTxFile2 );
# $cmd = "rm -f $finalTxFile; ln -s $ap $finalTxFile";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

## comparison between old and new species assignments
print "--- Comparing old and new species assignments\n";
$cmd = "cmp_tx.pl $quietStr -i $txFile -j $updatedTxFile2 -o old_vs_new_spp_cmp";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


##
## phylogenetic tree with final taxonomy
##

print "--- Generating a tree with final species names at leaves\n";
my $finalSppTreeFile = "$grPrefix" . "_final_spp.tree";
$cmd = "rm -f $finalSppTreeFile; nw_rename $treeFile $updatedTxFile2 | nw_order -c n  - > $finalSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


print "--- Generating tree with <final species name>_<seqID> labels at leaves\n";
my $finalSppSeqIDsFile = "$grPrefix" . "_final_spp.seqIDs";
$cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $updatedTxFile2 > $finalSppSeqIDsFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $finalSppSeqIdTreeFile = "$grPrefix" . "_final_sppSeqIDs.tree";
$cmd = "rm -f $finalSppSeqIdTreeFile; nw_rename $treeFile $finalSppSeqIDsFile | nw_order -c n  -  > $finalSppSeqIdTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final species clades collapsed to a single node \n";
my $finalCondSppTreeFile = abs_path( $grPrefix . "_final_spp_condensed.tree" );
$cmd = "rm -f $finalCondSppTreeFile; nw_condense $finalSppTreeFile > $finalCondSppTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


$section = qq~

##
## Genus-level cleanup
##

~;
print "$section";


#
# Using phyloPart followed by vicut to generate genus taxonomy
#
# Given a condenset tree at some taxonomic rank (say species), a percentile
# threshold and a taxonomic parent table, the script performs a partition of a
# phylogenetic tree with the given percentile distance threshold using phyloPart
# and then does vicut clustering using all except 0 clusters of phyloPart as
# annotation and elements of cluster 0 as query leaves and uses the taxonomic
# parent table to assign taxonomically relevant names to the clusters.

print "--- Running cluster_taxons.pl on condensed species tree\n";

## Synching %spParent with $finalCondSppTreeFile leaves
## right now spParent may have more species as the keys

$cmd = "rm -f $treeLeavesFile; nw_labels -I $finalCondSppTreeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

@treeLeaves = readArray($treeLeavesFile);

my @spParentKeys = keys %spParent;
my @excessSpp = diff(\@spParentKeys, \@treeLeaves);
if (@excessSpp)
{
  if ($debug)
  {
    print "\n\nWARNING: Detected spParent keys that are not tree leaves\n";
    print "Here is the list of deleted keys\n";
    printArray(\@excessSpp);
    print "\n\tDeleting excess species from spParent table\n\n";
  }
  delete @spParent{@excessSpp};
}

@spParentKeys = keys %spParent;

my $spParentFile = $grPrefix . ".spParent";
writeTbl(\%spParent, $spParentFile);

if ( @spParentKeys != @treeLeaves )
{
  warn "\n\n\tERROR: number of spParent species and $finalCondSppTreeFile leaves are not the same";

  print "\n\tspParentKeys: " . @spParentKeys . "\n";
  print "\ttreeLeaves: " . @treeLeaves . "\n";

  print "\n\tspParent tbl written to $spParentFile\n";

  my @excessLeaves = diff(\@treeLeaves, \@spParentKeys);
  print "\n\tHere is the list of tree leaves (species) that are missing in spParent table\n";
  printArray(\@excessLeaves);

  exit;
}


if ($debug)
{
  print "\nspParent table:\n";
  printFormatedTbl(\%spParent);
  print "\n\n";
}

my $sppGenusFile = $grPrefix . "_spp.genusTx";
$cmd = "cluster_taxons.pl $quietStr $igsStr $johannaStr $debugStr $showAllTreesStr -i $finalCondSppTreeFile -p 0.1 -f $spParentFile -t species -o $sppGenusFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my %sppGenus = readTbl($sppGenusFile);

print "--- Updating lineage table at the genus level\n";
my $finalGenusTxFile = $grPrefix . "_final_genus.tx";
open OUT, ">$finalGenusTxFile" or die "Cannot open $finalGenusTxFile for writing: $OS_ERROR\n";
my %geParent;
# my $outFile = "tmp.txt";
# open TEMPOUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
for my $id (keys %lineageTbl)
{
  my $lineage = $lineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;

  if ( exists $sppGenus{$sp} )
  {
    my $ge = pop @f;
    my $newGenus = $sppGenus{$sp};
    print OUT "$id\t$newGenus\n";

    if ( $newGenus ne $ge )
    {
      my @t = (@f, $newGenus, $sp);
      $lineageTbl{$id} = join ";", @t;

      # print TEMPOUT "\nge: $ge\tnewGe: $newGenus\n";
      # print TEMPOUT "lineage BEFORE: $lineage\n";
      # print TEMPOUT "lineage AFTER: " . $lineageTbl{$id} . "\n";
    }

    my $fa = pop @f;
    $geParent{$newGenus} = $fa;
  }
  else
  {
    warn "\n\n\tERROR: $id with sp: $sp not detected in sppGenus table";
    print "\n\n";
    exit;
    #delete $lineageTbl{$id};
  }
}
close OUT;
#close TEMPOUT;

my %geChildren;
for my $sp (keys %sppGenus)
{
  my $ge = $sppGenus{$sp};
  $geChildren{$ge}{$sp}++;
}

if ($debug)
{
  print "\n\nNumber of phylo-partition-vicut based genera: " . scalar(keys %geChildren) . "\n";
  print "\nSpecies frequencies in phylo-partition-vicut based genera:\n";
  my @a = sort { scalar(keys %{$geChildren{$b}}) <=> scalar(keys %{$geChildren{$a}}) } keys %geChildren;
  printFormatedTableValuedTbl(\%geChildren, \@a);
  print "\n\n"
}

print "--- Generating a tree with final genus names at leaves\n";
my $finalGenusTreeFile = "$grPrefix" . "_final_genus.tree";
$cmd = "rm -f $finalGenusTreeFile; nw_rename $treeFile $finalGenusTxFile | nw_order -c n  - > $finalGenusTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating a condensed tree with final genera collapsed to a single node \n";
my $finalCondGenusTreeFile = abs_path( "$grPrefix" . "_final_genus_condensed.tree" );
$cmd = "rm -f $finalCondGenusTreeFile; nw_condense $finalGenusTreeFile > $finalCondGenusTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


my $pdfCondGenusTreeFile = abs_path( "$grPrefix" . "_final_genus_condensed_tree.pdf" );
plot_tree_bw($finalCondGenusTreeFile, $pdfCondGenusTreeFile);

if ( $showAllTrees &&  $OSNAME eq "darwin")
{
  $cmd = "open $pdfCondGenusTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## generating pdf figure of species condensed tree colored by genus

# Assigning to each genus an index
my @genera = values %sppGenus;
@genera = unique(\@genera);
my $generaCount = 1;
my %genusIdx = map{ $_ => $generaCount++ } @genera;

if ($debug)
{
  print "\n\ngenusIdx:\n";
  printFormatedTbl(\%genusIdx);
  print "\n\n";
}

my $genusIdxFile = abs_path( "$grPrefix" . "_genus.idx" );
open OUT, ">$genusIdxFile" or die "Cannot open $genusIdxFile for writing: $OS_ERROR\n";
for my $sp (keys %sppGenus)
{
  print OUT "$sp\t" . $genusIdx{$sppGenus{$sp}} . "\n";
}
close OUT;

my $pdfCsppWithGenusColorsTreeFile = abs_path( $grPrefix . "_final_spp_condensed_tree_with_genus_colors.pdf" );

my $title = $grPrefix . " - genera";
plot_tree($finalCondSppTreeFile, $genusIdxFile, $pdfCsppWithGenusColorsTreeFile, $title);

if ( !$doNotPopPDFs && $OSNAME eq "darwin")
{
  $cmd = "open $pdfCsppWithGenusColorsTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


print "--- Splitting phylo-vicut-based genera if their sizes are above taxonSizeThld\n";

my %subSppGenus = %sppGenus;
my $detectedLargeGenus = 0;
my @a = sort { scalar(keys %{$geChildren{$b}}) <=> scalar(keys %{$geChildren{$a}}) } keys %geChildren;
for my $ge (@a)
{
  if ( scalar(keys %{$geChildren{$ge}}) > $taxonSizeThld ) # try to split this cluster into smaller pieces
  {
    $detectedLargeGenus = 1;
    print "\n\tSplitting $ge\n" if $debug;
    my @spp = keys %{$geChildren{$ge}};
    print "Species of $ge:\n" if $debug;
    printArray(\@spp) if $debug;

    my $geStr = $ge;
    if ( length($geStr) > $maxStrLen )
    {
      $geStr = substr( $ge, 0, $maxStrLen);
    }

    # prune the final species condensed tree to contain the given genus only
    my $prunedTreeFile = $grPrefix . "_pruned_$geStr.tree";
    print "--- Restricting $finalCondSppTreeFile to the species of $ge\n";
    $cmd = "rm -f $prunedTreeFile; nw_clade $finalCondSppTreeFile @spp > $prunedTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    # testing if the resulting tree has the same number of leaves as @spp array
    my $prunedTreeLeavesFile = $grPrefix . "_pruned_$geStr.leaves";
    $cmd = "rm -f $prunedTreeLeavesFile; nw_labels -I $prunedTreeFile > $prunedTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    my @cladeLeaves = readArray($prunedTreeLeavesFile);

    my @commSpp = comm(\@cladeLeaves, \@spp);
    if (@commSpp != @cladeLeaves || @commSpp != @spp)
    {
      warn "ERROR: leaves of $prunedTreeLeavesFile and genus $ge species do not match";
      print "Number of leaves of $prunedTreeLeavesFile: " . @cladeLeaves . "\n";
      print "Number of species in $ge: " . @spp . "\n";
      print "Number of common species: " . @commSpp . "\n";
      exit;
    }

    # make this a rooted tree
    $cmd = "root_tree.pl -i $prunedTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    # Creating parent table
    # Assigning to each species its genus.

    print "\nCreating spParent2 table\n" if $debug;
    $spParentFile = $grPrefix . "_$geStr.spParent";
    my %spParent2;
    for my $sp (@spp)
    {
      #print "\nsp: $sp\n" if $debug;
      my ($g, $s) = split "_", $sp;
      if (!defined $g)
      {
	print "Genus undef for $sp\n" if $debug;
	exit;
      }
      else
      {
	$spParent2{$sp} = $g;
	#print "g: $g\n" if $debug;
      }
    }

    print "Writing spParent2 to $spParentFile\n" if $debug;
    writeTbl(\%spParent2, $spParentFile);

    $sppGenusFile = $grPrefix . "_$geStr" . "_spp.genusTx";
    $cmd = "cluster_taxons.pl $quietStr $igsStr $johannaStr $debugStr $showAllTreesStr -i $prunedTreeFile -p 0.1 -f $spParentFile -t species -o $sppGenusFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    my %sppGenus2 = readTbl($sppGenusFile);

    ## If the genus has a numeric suffix, idx, we will have to change genus
    ## assignment after splitting this genus to something like
    ## <genus>_<idx>s<newIdx>.

    my $coreGenusName = $ge;
    $coreGenusName =~ s/_(\d+)//;
    # print "coreGenusName: $coreGenusName\n" if $debug;
    # if ($1)
    # {
    #   print "Index suffix: $1\n" if $debug;
    # }
    # else
    # {
    #   print "Index suffix undef\n" if $debug;
    # }

    if ($1)
    {
      # Changing cluster IDs to <genus>_<idx>s<newIdx>
      for my $sp (keys %sppGenus2)
      {
	my $g = $sppGenus2{$sp};
	my @f = split "_", $g;
	my $i = pop @f;
	my $b = join "_", @f;
	$sppGenus2{$sp} = $b . "_$1" . "s$i";
      }
    }

    my %geChildren2;
    for my $sp (keys %sppGenus2)
    {
      my $ge = $sppGenus2{$sp};
      $geChildren2{$ge}{$sp}++;
    }

    if ($debug)
    {
      print "\nNumber of phylo-partition-vicut based sub-genera of $ge: " . scalar(keys %geChildren2) . "\n";
      print "\nSpecies frequencies in phylo-partition-vicut based sub-genera of $ge:\n";
      my @a = sort { scalar(keys %{$geChildren2{$b}}) <=> scalar(keys %{$geChildren2{$a}}) } keys %geChildren2;
      printFormatedTableValuedTbl(\%geChildren2, \@a);
      print "\n\n"
    }

    @subSppGenus{@spp} = @sppGenus2{@spp};

    print "--- Updating lineage table at the sub-genus of $ge level\n" if $debug;
    for my $id (keys %lineageTbl)
    {
      my $lineage = $lineageTbl{$id};
      my @f = split ";", $lineage;
      my $sp = pop @f;

      if ( exists $sppGenus2{$sp} )
      {
	##my $ge = pop @f; # keeping genus
	my $newSubGenus = $sppGenus2{$sp};
	my @t = (@f, $newSubGenus, $sp);
	$lineageTbl{$id} = join ";", @t;

	my $ge = pop @f;
	#$subgeParent{$newSubGenus} = $ge;
      }
    }

  }
  else
  {
    my @spp = keys %{$geChildren{$ge}};
    my %sppInd = map{$_ =>1} @spp;

    print "--- Updating lineage table at the sub-genus level for small genus $ge\n" if $debug;
    for my $id (keys %lineageTbl)
    {
      my $lineage = $lineageTbl{$id};
      my @f = split ";", $lineage;
      my $sp = pop @f;

      if ( exists $sppInd{$sp} )
      {
	my $ge = pop @f; # keeping genus
	my @t = (@f, $ge, $ge, $sp);
	$lineageTbl{$id} = join ";", @t;
	#$subgeParent{$ge} = $ge;
      }
    }
  }

}

if ($detectedLargeGenus)
{
  my %subgeChildren;
  for my $sp (keys %subSppGenus)
  {
    my $ge = $subSppGenus{$sp};
    $subgeChildren{$ge}{$sp}++;
  }

  if ($debug)
  {
    print "\n\nAFTER splitting of genera the number of phylo-partition-vicut based sub-genera: " . scalar(keys %subgeChildren) . "\n";
    print "\nSpecies frequencies in the new phylo-partition-vicut based genera:\n";
    my @a = sort { scalar(keys %{$subgeChildren{$b}}) <=> scalar(keys %{$subgeChildren{$a}}) } keys %subgeChildren;
    printFormatedTableValuedTbl(\%subgeChildren, \@a);
    print "\n\n"
  }

  my $finalSubGenusTxFile = $grPrefix . "_final_sub_genus.tx";
  open OUT, ">$finalSubGenusTxFile" or die "Cannot open $finalSubGenusTxFile for writing: $OS_ERROR\n";
  for my $id (keys %lineageTbl)
  {
    my $lineage = $lineageTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $subge = pop @f;
    print OUT "$id\t$subge\n";
  }
  close OUT;

  print "--- Generating a tree with final sub genus names at leaves\n";
  my $finalSubGenusTreeFile = "$grPrefix" . "_final_sub_genus.tree";
  $cmd = "rm -f $finalSubGenusTreeFile; nw_rename $treeFile $finalSubGenusTxFile | nw_order -c n  - > $finalSubGenusTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Generating a condensed tree with final genera collapsed to a single node \n";
  my $finalCondSubGenusTreeFile = abs_path( "$grPrefix" . "_final_sub_genus_condensed.tree" );
  $cmd = "rm -f $finalCondSubGenusTreeFile; nw_condense $finalSubGenusTreeFile > $finalCondSubGenusTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $pdfCondSubGenusTreeFile = abs_path( "$grPrefix" . "_final_sub_genus_condensed_tree.pdf" );
  plot_tree_bw($finalCondSubGenusTreeFile, $pdfCondSubGenusTreeFile);

  if ( $showAllTrees && $OSNAME eq "darwin")
  {
    $cmd = "open $pdfCondSubGenusTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }


  ## generating pdf figure of species condensed tree colored by sub-genus

  # Assigning to each sub-genus an index
  my @subgenera = values %subSppGenus;
  @subgenera = unique(\@subgenera);
  my $subgeneraCount = 1;
  my %subgenusIdx = map{ $_ => $subgeneraCount++ } @subgenera;

  if ($debug)
  {
    print "\n\nsubgenusIdx:\n";
    printFormatedTbl(\%subgenusIdx);
    print "\n\n";
  }

  my $subgenusIdxFile = abs_path( "$grPrefix" . "_subgenus.idx" );
  open OUT, ">$subgenusIdxFile" or die "Cannot open $subgenusIdxFile for writing: $OS_ERROR\n";
  for my $sp (keys %subSppGenus)
  {
    print OUT "$sp\t" . $subgenusIdx{$subSppGenus{$sp}} . "\n";
  }
  close OUT;

  my $pdfCsppTreeFile = abs_path( $grPrefix . "_final_spp_condensed_tree_with_sub_genus_colors.pdf" );

  $title = $grPrefix . " - sub-genera";
  plot_tree($finalCondSppTreeFile, $subgenusIdxFile, $pdfCsppTreeFile, $title);

  if ( !$doNotPopPDFs && $OSNAME eq "darwin")
  {
    $cmd = "open $pdfCsppTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }
}
else
{
  print "\n\n\tNo large genera detected\n\n";
}


$section = qq~

##
## family-level cleanup
##

~;
print "$section";

print "--- Running cluster_taxons.pl on condensed genus tree\n";

my $geParentFile = $grPrefix . ".geParent";
writeTbl(\%geParent, $geParentFile);

if ($debug)
{
  print "\ngeParent table:\n";
  printFormatedTbl(\%geParent);
  print "\n\n";
}

my $genusFamilyFile = $grPrefix . "_genus.familyTx";
$cmd = "cluster_taxons.pl $quietStr $igsStr $johannaStr $debugStr $showAllTreesStr -i $finalCondGenusTreeFile -p 0.1 -f $geParentFile -t genus -o $genusFamilyFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my %genusFamily = readTbl($genusFamilyFile);

if ($debug)
{
  print "\n\ngenus => family phyloPart-vicut table:\n";
  printFormatedTbl(\%genusFamily);
  print "\n\n";
}

print "--- Updating lineage table at the family level\n";
my $finalFamilyTxFile = $grPrefix . "_final_family.tx";
open OUT, ">$finalFamilyTxFile" or die "Cannot open $finalFamilyTxFile for writing: $OS_ERROR\n";
my %faParent;
for my $id (keys %lineageTbl)
{
  my $lineage = $lineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $sge = pop @f;
  my $ge = pop @f;

  if ( exists $genusFamily{$ge} )
  {
    my $fa = pop @f;
    my $newFamily = $genusFamily{$ge};
    print OUT "$id\t$newFamily\n";
    my @t = (@f, $newFamily, $ge, $sge, $sp);
    $lineageTbl{$id} = join ";", @t;

    my $or = pop @f;
    $faParent{$newFamily} = $or;
  }
  else
  {
    warn "\n\n\tERROR: $id with ge: $ge not detected in genusFamily table";
    print "\n\n";
    exit;
  }
}
close OUT;

my %faChildren;
for my $ge (keys %genusFamily)
{
  my $fa = $genusFamily{$ge};
  $faChildren{$fa}{$ge}++;
}

if ($debug)
{
  print "\n\nNumber of phylo-partition-vicut based families: " . scalar(keys %faChildren) . "\n";
  print "\nGenus frequencies in phylo-partition-vicut based families:\n";
  my @a = sort { scalar(keys %{$faChildren{$b}}) <=> scalar(keys %{$faChildren{$a}}) } keys %faChildren;
  printFormatedTableValuedTbl(\%faChildren, \@a);
  print "\n\n"
}


if ( scalar(keys %faChildren) > 1 )
{
  print "--- Generating a tree with final family names at leaves\n";
  my $finalFamilyTreeFile = "$grPrefix" . "_final_family.tree";
  $cmd = "rm -f $finalFamilyTreeFile; nw_rename $treeFile $finalFamilyTxFile | nw_order -c n  - > $finalFamilyTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "--- Generating a condensed tree with final families collapsed to a single node \n";
  my $finalCondFamilyTreeFile = abs_path( "$grPrefix" . "_final_family_condensed.tree" );
  $cmd = "rm -f $finalCondFamilyTreeFile; nw_condense $finalFamilyTreeFile > $finalCondFamilyTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $pdfCondFamilyTreeFile = abs_path( "$grPrefix" . "_final_family_condensed_tree.pdf" );
  plot_tree_bw($finalCondFamilyTreeFile, $pdfCondFamilyTreeFile);

  if ( $showAllTrees &&  $OSNAME eq "darwin")
  {
    $cmd = "open $pdfCondFamilyTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }


  ## generating pdf figure of species condensed tree colored by families

  # Assigning to each family an index
  my @families = values %genusFamily;
  @families = unique(\@families);
  my $familiesCount = 1;
  my %familyIdx = map{ $_ => $familiesCount++ } @families;

  if ($debug)
  {
    print "\n\nfamilyIdx:\n";
    printFormatedTbl(\%familyIdx);
    print "\n\n";
  }

  my $familyIdxFile = abs_path( "$grPrefix" . "_family.idx" );
  open OUT, ">$familyIdxFile" or die "Cannot open $familyIdxFile for writing: $OS_ERROR\n";
  for my $sp (keys %sppGenus)
  {
    my $ge = $sppGenus{$sp};
    print OUT "$sp\t" . $familyIdx{$genusFamily{$ge}} . "\n";
  }
  close OUT;

  my $pdfFamilyColorsCsppTreeFile = abs_path( $grPrefix . "_final_species_condensed_tree_with_family_colors.pdf" );

  $title = $grPrefix . " - families";
  plot_tree($finalCondSppTreeFile, $familyIdxFile, $pdfFamilyColorsCsppTreeFile, $title);

  if ( !$doNotPopPDFs && $OSNAME eq "darwin")
  {
    $cmd = "open $pdfFamilyColorsCsppTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }

  $section = qq~

##
## Order-level cleanup
##

~;
  print "$section";

  print "--- Running cluster_taxons.pl on condensed family tree\n";

  my $faParentFile = $grPrefix . ".faParent";
  writeTbl(\%faParent, $faParentFile);

  my $familyOrderFile = $grPrefix . "_family.orderTx";
  $cmd = "cluster_taxons.pl $quietStr $igsStr $johannaStr $debugStr $showAllTreesStr -i $finalCondFamilyTreeFile -p 0.1 -f $faParentFile -t family -o $familyOrderFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my %familyOrder = readTbl($familyOrderFile);

  if ($debug)
  {
    print "\n\nfamily => order table:\n";
    printFormatedTbl(\%familyOrder);
    print "\n\n";
  }

  print "--- Updating lineage table at the order level\n";
  my %orParent;

  # testing if the resulting clustering does not consists of singlentons
  # if it does, keep lineage as it is and stop the taxonomy update
  my $nFamilies = keys %faParent;
  my @orders = values %familyOrder;
  @orders = unique(\@orders);

  my $finalOrderTxFile = $grPrefix . "_final_order.tx";
  if ( @orders < $nFamilies )
  {
    open OUT, ">$finalOrderTxFile" or die "Cannot open $finalOrderTxFile for writing: $OS_ERROR\n";
    for my $id (keys %lineageTbl)
    {
      my $lineage = $lineageTbl{$id};
      my @f = split ";", $lineage;
      my $sp = pop @f;
      my $sge = pop @f;
      my $ge = pop @f;
      my $fa = pop @f;

      if ( exists $familyOrder{$fa} )
      {
	my $or = pop @f;
	my $newOrder = $familyOrder{$fa};
	print OUT "$id\t$newOrder\n";
	my @t = (@f, $newOrder, $fa, $ge, $sge, $sp);
	$lineageTbl{$id} = join ";", @t;

	my $cl = pop @f;
	$orParent{$newOrder} = $cl;
      }
      else
      {
	warn "\n\n\tERROR: $id with fa: $fa not detected in familyOrder table";
	print "\n\n";
	exit;
      }
    }
    close OUT;
  }

  my %orChildren;
  for my $fa (keys %familyOrder)
  {
    my $or = $familyOrder{$fa};
    $orChildren{$or}{$fa}++;
  }

  if ($debug)
  {
    print "\n\nNumber of phylo-partition-vicut based orders: " . scalar(keys %orChildren) . "\n";
    print "\nFamily frequencies in phylo-partition-vicut based orders:\n";
    my @a = sort { scalar(keys %{$orChildren{$b}}) <=> scalar(keys %{$orChildren{$a}}) } keys %orChildren;
    printFormatedTableValuedTbl(\%orChildren, \@a);
    print "\n\n"
  }


  if ( (scalar(keys %orChildren) > 1) && (@orders < $nFamilies) )
  {
    print "--- Generating a tree with final order names at leaves\n";
    my $finalOrderTreeFile = "$grPrefix" . "_final_order.tree";
    $cmd = "rm -f $finalOrderTreeFile; nw_rename $treeFile $finalOrderTxFile | nw_order -c n  - > $finalOrderTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    print "--- Generating a condensed tree with final families collapsed to a single node \n";
    my $finalCondOrderTreeFile = abs_path( "$grPrefix" . "_final_order_condensed.tree" );
    $cmd = "rm -f $finalCondOrderTreeFile; nw_condense $finalOrderTreeFile > $finalCondOrderTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    my $pdfCondOrderTreeFile = abs_path( "$grPrefix" . "_final_order_condensed_tree.pdf" );
    plot_tree_bw($finalCondOrderTreeFile, $pdfCondOrderTreeFile);

    if ( $showAllTrees &&  $OSNAME eq "darwin")
    {
      $cmd = "open $pdfCondOrderTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }


    ## generating pdf figure of species condensed tree colored by orders

    # Assigning to each order an index
    my $ordersCount = 1;
    my %orderIdx = map{ $_ => $ordersCount++ } @orders;

    if ($debug)
    {
      print "\n\norderIdx:\n";
      printFormatedTbl(\%orderIdx);
      print "\n\n";
    }

    my $orderIdxFile = abs_path( "$grPrefix" . "_order.idx" );
    open OUT, ">$orderIdxFile" or die "Cannot open $orderIdxFile for writing: $OS_ERROR\n";
    for my $sp (keys %sppGenus)
    {
      my $ge = $sppGenus{$sp};
      print OUT "$sp\t" . $orderIdx{$familyOrder{$genusFamily{$ge}}} . "\n";
    }
    close OUT;

    my $pdfOrderColorsCsppTreeFile = abs_path( $grPrefix . "_final_species_condensed_tree_with_order_colors.pdf" );

    $title = $grPrefix . " - orders";
    plot_tree($finalCondSppTreeFile, $orderIdxFile, $pdfOrderColorsCsppTreeFile, $title);

    if ( !$doNotPopPDFs && $OSNAME eq "darwin")
    {
      $cmd = "open $pdfOrderColorsCsppTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    $section = qq~

##
## Class-level cleanup
##

~;
    print "$section";

    print "--- Running cluster_taxons.pl on condensed order tree\n";

    my $orParentFile = $grPrefix . ".orParent";
    writeTbl(\%orParent, $orParentFile);

    my $orderClassFile = $grPrefix . "_order.classTx";
    $cmd = "cluster_taxons.pl $quietStr $igsStr $johannaStr $debugStr $showAllTreesStr -i $finalCondOrderTreeFile -p 0.1 -f $orParentFile -t order -o $orderClassFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    my %orderClass = readTbl($orderClassFile);

    if ($debug)
    {
      print "\n\norder => class table:\n";
      printFormatedTbl(\%orderClass);
      print "\n\n";
    }

    print "--- Updating lineage table at the class level\n";
    my %clParent;

    # testing if the resulting clustering does not consists of singlentons
    # if it does, keep lineage as it is and stop the taxonomy update
    my $nOrders = keys %orParent;
    my @classes = values %orderClass;
    @classes = unique(\@classes);

    my $finalClassTxFile = $grPrefix . "_final_class.tx";
    if ( @classes < $nOrders )
    {
      open OUT, ">$finalClassTxFile" or die "Cannot open $finalClassTxFile for writing: $OS_ERROR\n";
      for my $id (keys %lineageTbl)
      {
	my $lineage = $lineageTbl{$id};
	my @f = split ";", $lineage;
	my $sp = pop @f;
	my $sge = pop @f;
	my $ge = pop @f;
	my $fa = pop @f;
	my $or = pop @f;

	if ( exists $orderClass{$or} )
	{
	  my $cl = pop @f;
	  my $newClass = $orderClass{$or};
	  print OUT "$id\t$newClass\n";
	  my @t = (@f, $newClass, $or, $fa, $ge, $sge, $sp);
	  $lineageTbl{$id} = join ";", @t;

	  my $ph = pop @f;
	  $clParent{$newClass} = $ph;
	}
	else
	{
	  warn "\n\n\tERROR: $id with or: $or not detected in orderClass table";
	  print "\n\n";
	  exit;
	}
      }
      close OUT;
    }

    my %clChildren;
    for my $or (keys %orderClass)
    {
      my $cl = $orderClass{$or};
      $clChildren{$cl}{$or}++;
    }

    if ($debug)
    {
      print "\n\nNumber of phylo-partition-vicut based classes: " . scalar(keys %clChildren) . "\n";
      print "\nOrder frequencies in phylo-partition-vicut based classes:\n";
      my @a = sort { scalar(keys %{$clChildren{$b}}) <=> scalar(keys %{$clChildren{$a}}) } keys %clChildren;
      printFormatedTableValuedTbl(\%clChildren, \@a);
      print "\n\n"
    }


    if ( scalar(keys %clChildren) > 1 )
    {
      print "--- Generating a tree with final class names at leaves\n";
      my $finalClassTreeFile = "$grPrefix" . "_final_class.tree";
      $cmd = "rm -f $finalClassTreeFile; nw_rename $treeFile $finalClassTxFile | nw_order -c n  - > $finalClassTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

      print "--- Generating a condensed tree with final families collapsed to a single node \n";
      my $finalCondClassTreeFile = abs_path( "$grPrefix" . "_final_class_condensed.tree" );
      $cmd = "rm -f $finalCondClassTreeFile; nw_condense $finalClassTreeFile > $finalCondClassTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

      my $pdfCondClassTreeFile = abs_path( "$grPrefix" . "_final_class_condensed_tree.pdf" );
      plot_tree_bw($finalCondClassTreeFile, $pdfCondClassTreeFile);

      if ( $showAllTrees &&  $OSNAME eq "darwin")
      {
	$cmd = "open $pdfCondClassTreeFile";
	print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
      }


      ## generating pdf figure of species condensed tree colored by classes

      # Assigning to each class an index
      my $classesCount = 1;
      my %classIdx = map{ $_ => $classesCount++ } @classes;

      if ($debug)
      {
	print "\n\nclassIdx:\n";
	printFormatedTbl(\%classIdx);
	print "\n\n";
      }

      my $classIdxFile = abs_path( "$grPrefix" . "_class.idx" );
      open OUT, ">$classIdxFile" or die "Cannot open $classIdxFile for writing: $OS_ERROR\n";
      for my $sp (keys %sppGenus)
      {
	my $ge = $sppGenus{$sp};
	print OUT "$sp\t" . $classIdx{$orderClass{$familyOrder{$genusFamily{$ge}}}} . "\n";
      }
      close OUT;

      my $pdfClassColorsCsppTreeFile = abs_path( $grPrefix . "_final_species_condensed_tree_with_class_colors.pdf" );

      $title = $grPrefix . " - classes";
      plot_tree($finalCondSppTreeFile, $classIdxFile, $pdfClassColorsCsppTreeFile, $title);

      if ( !$doNotPopPDFs && $OSNAME eq "darwin")
      {
	$cmd = "open $pdfClassColorsCsppTreeFile";
	print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
      }
    }
  }
} # end if ( scalar(keys %orChildren) > 1 )


#   $section = qq~

# ##
# ## Generating final lineage, spLineage, fasta and taxon files
# ##

# ~;
#   print $section;

my $initNumSpecies  = scalar( keys %spTbl );
my $initNumGenera   = scalar( keys %geTbl );
my $initNumFamilies = scalar( keys %faTbl );
my $initNumOrders   = scalar( keys %orTbl );
my $initNumClasses  = scalar( keys %clTbl );
my $initNumPhyla    = scalar( keys %phTbl );

undef %spTbl;
undef %geTbl;
my %subGeTbl;
undef %faTbl;
undef %orTbl;
undef %clTbl;
undef %phTbl;
undef %children;
undef %parent;
undef %spLineage;

for my $id ( keys %lineageTbl )
{
  my $lineage = $lineageTbl{$id};
  my @f = split ";", $lineage;
  my $sp = pop @f;
  my $subGe = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;
  my $or = pop @f;
  my $cl = pop @f;
  my $ph = pop @f;

  ##$sp = "s_$sp";
  $sp .= "_OG" if ( exists $ogInd{$id} );
  $subGe = "sg_$subGe";
  $ge = "g_$ge";
  $fa = "f_$fa";
  $or = "o_$or";
  $cl = "c_$cl";
  $ph = "p_$ph";

  $spLineage{$sp} = "$subGe\t$ge\t$fa\t$or\t$cl\t$ph\td_Bacteria";

  $parent{$sp} = $subGe;
  $parent{$subGe} = $ge;
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
  $children{$ge}{$subGe}++;
  $children{$subGe}{$sp}++;

  push @{$spTbl{$sp}}, $id;
  push @{$subGeTbl{$subGe}}, $id;
  push @{$geTbl{$ge}}, $id;
  push @{$faTbl{$fa}}, $id;
  push @{$orTbl{$or}}, $id;
  push @{$clTbl{$cl}}, $id;
  push @{$phTbl{$ph}}, $id;
}

if ($debug)
{
  print "\nTaxonomy AFTER cleanup\n";
  printLineage2();
}
printLineageToFile2($SRYOUT, "\n\n====== Taxonomy AFTER cleanup ======\n");

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
  my $sge = pop @f; # getting rid of sub-genus
  push @f, $newTxNoTGTs{$id};
  my $l = join ";", @f;
  print OUT "$id\t$l\n";
}
for my $id ( keys %ogLineageTbl )
{
  my $lineage = $ogLineageTbl{$id};
  print OUT "$id\t$lineage\n";
}
close OUT;

if ($buildModelData)
{
  # Generate
  #   - species lineage file (grPrefix_final.spLineage)
  #   - taxon file (grPrefix_final.tx)
  #   - ungapped fasta file corresponding to sequences present in the taxon file

  $section = qq~

##
## Building MC models
##

~;
  print $section;

  print "--- Creating final taxonomy file\n";
  my $txFile = $grPrefix . "_final.tx";
  $cmd = "rm -f $txFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  open OUT, ">$txFile" or die "Cannot open $txFile for writing: $OS_ERROR\n";
  for my $id (keys %newTx)
  {
    print OUT "$id\t" . $newTx{$id} . "\n";
  }
  close OUT;

  print "--- Creating species lineage file\n";
  my $spLineageFile = $grPrefix . "_final.spLineage";
  open OUT, ">$spLineageFile" or die "Cannot open $spLineageFile for writing: $OS_ERROR\n";
  for my $sp (keys %spLineage)
  {
    print OUT "$sp\t" . $spLineage{$sp} . "\n";
  }
  close OUT;

  print "--- Creating reference gap-free fasta file\n";
  my $tmpFile = $grPrefix . "_tmp.fa";
  $cmd = "rmGaps $quietStr -i $trimmedAlgnFile -o $tmpFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $faFile = $grPrefix . "_final.fa";
  $cmd = "select_seqs.pl $quietStr -s $txFile -i $tmpFile -o $faFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # testing if faFile and txFile have the same seq IDs
  my @faSeqIDs  = get_seqIDs_from_fa($faFile);
  my @newTxKeys = keys %newTx;

  my @commIDs = comm(\@newTxKeys, \@faSeqIDs);
  if (@commIDs != @newTxKeys || @commIDs != @faSeqIDs)
  {
    warn "\n\nERROR: seq IDs of new taxonomy table and the fasta file do not match";
    print "\n\tNumber of elements in the new taxonomy table: " . @newTxKeys . "\n";
    print "\tNumber of elements in the fasta file: " . @faSeqIDs . "\n";
    print "\tNumber of common elements: " . @commIDs . "\n";

    writeArray(\@newTxKeys, "newTxKeys.txt");
    writeArray(\@faSeqIDs, "faSeqIDs.txt");

    print "\n\tNew taxon keys and fasta IDs written to newTxKeys.txt and faSeqIDs.txt, respectively\n\n";

    if (@newTxKeys > @faSeqIDs)
    {
      my @d = diff(\@newTxKeys, \@faSeqIDs);
      print "\nElements in new taxonomy that are not in new fasta:\n";
      for (@d)
      {
	print "\t$_\t" . $newTx{$_} . "\n";
      }
      print "\n\n";
    }

    if (@faSeqIDs > @newTxKeys)
    {
      my @d = diff(\@faSeqIDs, \@newTxKeys);
      print "\nElements in new fasta that are not in the new taxonomy:\n";
      for (@d)
      {
	print "\t$_\t" . $lineageTbl{$_} . "\n";
      }
      print "\n\n";
    }

    exit;
  }

  # rm -rf Firmicutes_group_6_V3V4_MC_models_dir; buildModelTree -l Firmicutes_group_6_V3V4_final.spLineage -i Firmicutes_group_6_V3V4_final.fa -t Firmicutes_group_6_V3V4_final.tx -o Firmicutes_group_6_V3V4_MC_models_dir
  print "--- Building model tree and creating taxon's reference fasta files\n";
  my $mcDir = $grPrefix . "_MC_models_dir";
  $cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $spLineageFile -i $faFile -t $txFile -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # buildMC -t Firmicutes_group_6_V3V4_MC_models_dir/spp_paths.txt -k 8 -d Firmicutes_group_6_V3V4_MC_models_dir
  print "--- Building MC models\n";
  $cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # clError -d Firmicutes_group_6_V3V4_MC_models_dir -o Firmicutes_group_6_V3V4_MC_models_clError_dir
  print "--- Generating random sequences from the MC models\n";
  my $errorDir = $grPrefix . "_MC_models_clError_dir";
  $cmd = "rm -rf $errorDir; clError -d $mcDir -o $errorDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # generate ncProbThlds.txt
  print "--- Generating ncProbThlds.txt\n";
  build_clError(abs_path($mcDir), abs_path($errorDir));

  # classify  -d Firmicutes_group_6_V3V4_MC_models_dir -i Firmicutes_group_6_V3V4_final.fa -o Firmicutes_group_6_V3V4_MC_models_dir
  print "--- Running classify on $faFile\n";
  $cmd = "classify -d $mcDir -i $faFile -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  #$ classify -r vaginal_319_806_v2_dir/refTx.tree -d vaginal_319_806_v2_MCdir -i vaginal_319_806_v2.fa -o mcDir_319_806_v2
  #$ cmp_tx.pl -i vaginal_319_806_v2.tx -j mcDir_319_806_v2_no_err_thld//MC.order7.results.txt -o mcDir_319_806_v2_no_err_thld/
  print "--- Comparing ref seq's taxonomy with the classification results\n";
  $cmd = "cmp_tx.pl --verbose $quietStr -i $txFile -j $mcDir/MC_order7_results.txt -o $mcDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


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



## Number of species per sub-genus
my @subGenera = keys %subGeTbl;
my %subGenusSize; # <sub-genus> => <number of species of that genus>
for my $g (@subGenera)
{
  $subGenusSize{$g} = keys %{$children{$g}};
}
@subGenera = sort {$subGenusSize{$b} <=> $subGenusSize{$a}} keys %subGenusSize;
printFormatedTblToFile(\%subGenusSize, \@subGenera, $SRYOUT, "\n==== Number of species per sub-genus ====");
if ($debug)
{
  print "\n==== Number of species per sub-genus ====\n";
  printFormatedTbl(\%subGenusSize, \@subGenera);
}

## Number of sub-genera per genus
@genera = keys %geTbl;
my %genusSize; # <genus> => <number of species of that genus>
for my $g (@genera)
{
  $genusSize{$g} = keys %{$children{$g}};
}
@genera = sort {$genusSize{$b} <=> $genusSize{$a}} keys %genusSize;
printFormatedTblToFile(\%genusSize, \@genera, $SRYOUT, "\n==== Number of sub-genera per genus ====");
if ($debug)
{
  print "\n==== Number of sub-genera per genus ====\n";
  printFormatedTbl(\%genusSize, \@genera);
}

## Number of genera per order
my @families = keys %faTbl;
my %familySize; # <family> => <number of genera of that family>
for my $f (@families)
{
  $familySize{$f} = keys %{$children{$f}};
}
@families = sort {$familySize{$b} <=> $familySize{$a}} keys %familySize;
printFormatedTblToFile(\%familySize, \@families, $SRYOUT, "\n==== Number of genera per family ====");
if ($debug)
{
  print "\n==== Number of genera per family ====\n";
  printFormatedTbl(\%familySize, \@families);
}

## Number of families per order
my @orders = keys %orTbl;
my %orderSize; # <order> => <number of families of that order>
for my $o (@orders)
{
  $orderSize{$o} = keys %{$children{$o}};
}
@orders = sort {$orderSize{$b} <=> $orderSize{$a}} keys %orderSize;
printFormatedTblToFile(\%orderSize, \@orders, $SRYOUT, "\n==== Number of families per order ====");
if ($debug)
{
  print "\n==== Number of families per order ====\n";
  printFormatedTbl(\%orderSize, \@orders);
}

## Number of orders per class
my @classes = keys %clTbl;
my %classSize; # <class> => <number of orders of that class>
for my $c (@classes)
{
  $classSize{$c} = keys %{$children{$c}};
}
@classes = sort {$classSize{$b} <=> $classSize{$a}} keys %classSize;
printFormatedTblToFile(\%classSize, \@classes, $SRYOUT, "\n==== Number of orders per class ====");
if ($debug)
{
  print "\n==== Number of orders per class ====\n";
  printFormatedTbl(\%classSize, \@classes);
}

print $SRYOUT "\n\n--- Final summary\n";

print $SRYOUT  "\n\tNumber of species (with OG seq's) BEFORE taxonomic cleanup: " . scalar( keys %sppFreq ) . "\n";
print $SRYOUT    "\tNumber of species (with OG seq's) AFTER taxonomic cleanup and tentative ribotype splits:  " . scalar( keys %sppFreqFinal ) . "\n\n";

print $SRYOUT "\tNumber of _sp species BEFORE taxonomic cleanup: $nSpSpp\n";
print $SRYOUT "\tNumber of _sp species AFTER taxonomic cleanup: $nSpSppFinal\n\n";

print $SRYOUT  "\tNumber of singletons species BEFORE taxonomic cleanup: $sppFreq2{1}\n";
if (exists $sppFreqFinal2{1})
{
  print $SRYOUT  "\tNumber of singletons species AFTER taxonomic cleanup: $sppFreqFinal2{1}\n\n";
}
else
{
  print $SRYOUT  "\tNo singletons species AFTER taxonomic cleanup\n\n";
}

print $SRYOUT  "\tNumber of sub-genera AFTER taxonomic cleanup: " . scalar( keys %subGeTbl ) . "\n\n";

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


if ($debug)
{
  print  "\n\n--- Final summary\n";

  print   "\n\n\tNumber of species (with OG seq's) BEFORE taxonomic cleanup: " . scalar( keys %sppFreq ) . "\n";
  print     "\tNumber of species (with OG seq's) AFTER taxonomic cleanup and tentative ribotype splits:  " . scalar( keys %sppFreqFinal ) . "\n\n";

  print  "\tNumber of _sp species BEFORE taxonomic cleanup: $nSpSpp\n";
  print  "\tNumber of _sp species AFTER taxonomic cleanup: $nSpSppFinal\n\n";

  print   "\tNumber of singletons species BEFORE taxonomic cleanup: $sppFreq2{1}\n";
  if (exists $sppFreqFinal2{1})
  {
    print   "\tNumber of singletons species AFTER taxonomic cleanup: $sppFreqFinal2{1}\n\n";
  }
  else
  {
    print   "\tNo singletons species AFTER taxonomic cleanup\n\n";
  }

  print  "\tNumber of sub-genera AFTER taxonomic cleanup: " . scalar( keys %subGeTbl ) . "\n\n";

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
}

print "\n\n\tSee $grDir/$readmeFile for info about the taxonomy curation algorithm and the content of the output directory.\n";
print     "\tSummary stats written to $grDir/$summaryStatsFile\n\n";




####################################################################
##                               SUBS
####################################################################

sub plot_tree
{
  my ($treeFile, $clFile, $pdfFile, $title) = @_;

  my $showBoostrapVals = "F";

  if (!defined $title)
  {
    $title = $grPrefix;
  }

  my $Rscript = qq~

clTbl <- read.table(\"$clFile\", header=F)
str(clTbl)

cltr <- clTbl[,2]
names(cltr) <- clTbl[,1]

source(\"$readNewickFile\")
require(phytools)

tr1 <- read.newick(file=\"$treeFile\")
tr1 <- collapse.singles(tr1)

tip.cltr <- cltr[tr1\$tip.label]

colIdx <- 1
tip.colors <- c()
tip.colors[1] <- colIdx
for ( i in 2:length(tip.cltr) )
{
    if ( tip.cltr[i] != tip.cltr[i-1] )
    {
        colIdx <- colIdx + 1
        if ( colIdx==9 )
        {
            colIdx <- 1
        }
    }
    tip.colors[i] <- colIdx
    if ( colIdx==7 )
    {
        tip.colors[i] <- "brown" # using brown instead of yellow
    }
}

(nLeaves <- length(tr1\$tip.label))

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr1,type=\"phylogram\", no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8, tip.color=tip.colors, main=\"$title\")
par(op)
dev.off()
~;

  runRscript( $Rscript );
}


## plot tree with without any colors
sub plot_tree_bw
{
  my ($treeFile, $pdfFile, $title) = @_;

  my $showBoostrapVals = "F";

  if (!defined $title)
  {
    $title = $grPrefix;
  }

  my $Rscript = qq~

source(\"$readNewickFile\")
require(phytools)

tr1 <- read.newick(file=\"$treeFile\")
tr1 <- collapse.singles(tr1)

(nLeaves <- length(tr1\$tip.label))

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr1,type=\"phylogram\", no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8, main=\"$title\")
par(op)
dev.off()
~;

  runRscript( $Rscript );
}

  # execute an R-script
sub runRscript{

  my ($Rscript, $noErrorCheck) = @_;

  my $outFile = "rTmp.R";
  open OUT, ">$outFile",  or die "cannot write to $outFile: $!\n";
  print OUT "$Rscript";
  close OUT;

  my $cmd = "R CMD BATCH $outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

  if (!$noErrorCheck)
  {
    my $outR = $outFile . "out";
    open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
    my $exitStatus = 1;
    foreach my $line (<IN>)
    {
      if ( $line =~ /Error/ )
      {
	print "R script crashed at\n$line";
	print "check $outR for details\n";
	$exitStatus = 0;
	exit;
      }
    }
    close IN;
  }
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

## print lineage after adding subgenus taxon
sub printLineage2
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
	    my @subGes = keys %{$children{$ge}};
	    my $subgeKVL = getKeyValStrLengths(\%subGeTbl, \@subGes);
	    for my $sge ( sort{scalar(@{$subGeTbl{$a}}) <=> scalar(@{$subGeTbl{$b}})} @subGes)
	    {
	      printF(5, $sge, scalar(@{$subGeTbl{$sge}}), $subgeKVL);
	      my @sps = keys %{$children{$sge}};
	      my $spKVL = getKeyValStrLengths(\%spTbl, \@sps);
	      for my $sp ( sort{scalar(@{$spTbl{$a}}) <=> scalar(@{$spTbl{$b}})} @sps)
	      {
		printF(6, $sp, scalar(@{$spTbl{$sp}}), $spKVL);
	      }
	    }
	  }
	}
      }
    }
  }
}

sub printLineageToFile2
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
	    my @subGes = keys %{$children{$ge}};
	    my $subgeKVL = getKeyValStrLengths(\%subGeTbl, \@subGes);
	    for my $sge ( sort{scalar(@{$subGeTbl{$a}}) <=> scalar(@{$subGeTbl{$b}})} @subGes)
	    {
	      printF2(5, $sge, scalar(@{$subGeTbl{$sge}}), $subgeKVL, $fh);
	      my @sps = keys %{$children{$sge}};
	      my $spKVL = getKeyValStrLengths(\%spTbl, \@sps);
	      for my $sp ( sort{scalar(@{$spTbl{$a}}) <=> scalar(@{$spTbl{$b}})} @sps)
	      {
		printF2(6, $sp, scalar(@{$spTbl{$sp}}), $spKVL, $fh);
	      }
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
  print "\n";
}

# print elements of a hash table whose values are reference to a hash table so
# that sizes of the value hash tables are aligned
sub printFormatedTableValuedTbl{

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
    my $size = keys %{$rTbl->{$_}};
    print "$_$pad$size\n";
  }
  print "\n";
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTblToFile{

  my ($rTbl, $rSub, $fh, $msg) = @_; # the second argument is a subarray of the keys of the table

  print $fh "\n$msg\n";

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
    print $fh "$_$pad" . $rTbl->{$_} . "\n";
  }
  print $fh "\n";
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


# test if OG seq's form one or more clusters in the tree
sub test_OG
{
  my ($treeFile, $rogInd) = @_;

  my %ogInd = %{$rogInd};

  my $debug_test_OG = 0;

  print "\t--- Extracting leaves\n" if $debug_test_OG;
  my $treeLeavesFile = "$grPrefix" . "_sppSeqIDs.leaves";
  my $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "\t--- Reading leaves\n" if $debug_test_OG;
  my @leaves = readArray($treeLeavesFile);

  print "\t--- Checking the number of clusters formed by OG seqs\n" if $debug_test_OG;
  my @ogIdx;
  for my $i (0..$#leaves)
  {
    if ( exists $ogInd{$leaves[$i]} )
    {
      push @ogIdx, $i;
    }
  }

  printArray(\@ogIdx, "\nPositions of OG seq's") if ($debug_test_OG);

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

      if (0 && $debug_test_OG)
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

  printArrayByRow(\@og, "\nOutgroup elements") if ($debug_test_OG);

  if ( scalar(@start) != scalar(@end) )
  {
    warn "$grPrefix\nERROR: start and end arrays have different lengths!";
    print "length(start): " . @start . "\n";
    print "length(end): " . @end . "\n";
    exit;
  }

  my @rangeSize;
  for my $i (0..$#start)
  {
    push @rangeSize, ($end[$i] - $start[$i]+1);
  }

  if ($debug_test_OG)
  {
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";
  }

  if (@rangeSize>1)
  {
    warn "\n\n\tERROR: Detected multiple OG clusters";
    print "\n\n";

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

	#print "\t--- Extracting the clade of OG sequences\n";
	my $ogCladeTreeFile = "$grPrefix" . "_clade.tree";
	$cmd = "rm -f $ogCladeTreeFile; nw_clade $treeFile @og > $ogCladeTreeFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
	system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

	#print "\t--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$grPrefix" . "_clade.leaves";
	$cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
	system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

	#print "\t--- Reading the leaves\n" if $debug_test_OG;
	my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

	print "i: $i\t$start[$i]\t$end[$i]\t$rangeSize[$i]\t" . @ogCladeLeaves . "\n";
	if ( @ogCladeLeaves < $minCladeSize )
	{
	  $minCladeSize = @ogCladeLeaves;
	  $minCladeSizeIdx = $i;
	}
      }
    }

    # $imax = $minCladeSizeIdx;
    # print "\nUpdated imax: $imax\n";
    exit;
  }
  elsif ( !( $start[0] == 0 || $end[0] == $#leaves) )
  {
    warn "$\n\n\tERROR: In the pruned tree outgroups sequences are not at the top or bottom of the tree!";

    print "\n\nNumber of leaves: " . @leaves . "\n";
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

    print "\t--- Extracting the clade of OG sequences\n";
    my $ogCladeTreeFile = "$grPrefix" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug_test_OG;
    my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

    my $maxCladeSize = 100;
    if ( @ogCladeLeaves < $maxCladeSize )
    {
      printArrayByRow(\@ogCladeLeaves, "OG Clade Leaves");
    }
    else
    {
      print "\n\tLeaves of the OG clade written to $ogCladeTreeLeavesFile\n"
    }

    print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
    print   "\tNumber of OG sequences: " . @og . "\n\n";

    exit;
  }
  else
  {
    print "\t--- Extracting the clade of OG sequences\n" if $debug_test_OG;
    my $ogCladeTreeFile = "$grPrefix" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug_test_OG;
    my @ogCladeLeaves = readArray($ogCladeTreeLeavesFile);

    if ( @ogCladeLeaves != @og )
    {
      warn "\n\n\tERROR: The outgroup sequences do not form a monophyletic clade!";

      my $maxCladeSize = 100;
      if ( @ogCladeLeaves < $maxCladeSize )
      {
	printArrayByRow(\@ogCladeLeaves, "OG Clade Leaves");
      }
      else
      {
	print "\n\tLeaves of the OG clade written to $ogCladeTreeLeavesFile\n"
      }

      print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
      print   "\tNumber of OG sequences: " . @og . "\n\n";
      exit;
    }
  }

  print "\n\tOG seq's form a monophylectic clade at the top or bottom of the input tree\n\n" if $debug_test_OG;
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

sub build_clError
{
  my ($mcDir, $errorDir) = @_;

  my $Rscript = qq~

files <- list.files(path=\"$errorDir\", pattern=\"*.txt\", full.names=TRUE)
length(files)

files <- setdiff(files, c(\"spp_paths.txt\", \"modelIds.txt\"))
length(files)

(ii <- grep(\"error\", files))
length(ii)
if ( length(ii) )
    files <- files[-ii]
length(files)

(i <- grep(\"ncProbThlds.txt\", files))
length(i)
if ( length(i) )
files <- files[-i]
length(files)

thldTbl <- c()
errTbl <- c()
outDir <- \"$mcDir\"
(figFile <- paste(outDir,\"/errFig.pdf\",sep=\"\"))
pdf(figFile, width=11, height=11)
op <- par(mar=c(4,4,4.5,0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
for ( file in files )
{
    print(file)

    ## now, reading log posterior probabilities files is more complicated as for
    ## higher taxons the number of these log posteriors is greater than 1000, but
    ## for siblings its a 1000, so we have a 'table' with different number of
    ## columns.

    nFields <- count.fields(file, sep = \'\\t\')
    maxNumFields <- max(nFields)

    tbl <- read.table(file, sep=\"\\t\", col.names = paste0(\"V\",seq_len(maxNumFields)), fill = TRUE)
    dim(tbl)
    ids <- tbl[,1]
    refID <- ids[1]
    idx <- ids==refID
    sum(idx)
    x.ref <- as.numeric(as.matrix(tbl[idx,2:nFields[1]]))
    x.sib <- as.numeric(as.matrix(tbl[!idx,2:ncol(tbl)]))

    if ( length(x.ref)==0 )
    {
        stop(\"length(x.ref) is 0\")
    }
    else if ( length(x.sib)==0 )
    {
        stop(\"length(x.sib) is 0\")
    }

    id <- ids[1]

    d.ref <- density(x.ref)
    d.sib <- density(x.sib)

    xmin.ref <- min(x.ref)
    p0 <- min( c( (max(x.sib) + xmin.ref) / 2, xmin.ref ) )

    d.ref.fun <- approxfun(d.ref\$x, d.ref\$y, yleft=0, yright=0)
    d.sib.fun <- approxfun(d.sib\$x, d.sib\$y, yleft=0, yright=0)

    ff <- function(x) d.ref.fun(x)  - d.sib.fun(x)

    xmin <- min(c(d.ref\$x, d.sib\$x))
    xmax <- max(c(d.ref\$x, d.sib\$x))
    x <- seq(from=xmin, to=xmax, length=1000)
    y <- ff(x)
    r <- uniroot(ff, c(x[which.min(y)], x[which.max(y)]))
    ##r <- uniroot(ff, c( xmin.ref, d.ref\$x[which.max(d.ref\$y)]))
    p0 <- r\$root

    plot(x,y,main=id, type='l')
    abline(h=0, col='gray80'); abline(v=p0,col='gray80')
    abline(v=x[which.min(y)], col='red')
    abline(v=x[which.max(y)], col='red')

    if ( ff(x[which.max(y)]) > ff(x[which.min(y)]) )
    {
        dx <- diff(d.sib\$x)[1]
        idx <- d.sib\$x<0
        F1.sib <- sum(d.sib\$y[idx])*dx - cumsum(d.sib\$y[idx])*dx
        F1.sib.fun <- approxfun(d.sib\$x[idx], F1.sib, yleft=1, yright=0)
        x <- seq(p0, 0, length=100)
        tx <- cbind(x, F1.sib.fun(x))   # I don't think / F1.sib.fun(x[1])
                                        # normalization is correct. I should be
                                        # computing probabilities of hitting a
                                        # given value or more
    } else {
      tx <- cbind(p0, 0)
    }

    thldTbl[ids[1]] <- p0

    ## probability of an error = integral of d.sib from p0 to +inf
    errTbl[ids[1]] <- 0
    if ( p0 < max(d.sib\$x) )
    {
        e1 <- integrate(d.sib.fun, p0, max(d.sib\$x))[[1]]
        e2 <- integrate(d.ref.fun, min(d.ref\$x), p0)[[1]]
        errTbl[ids[1]] <- max(c(e1, e2))

        plot(d.ref, xlim=c(min(x.sib),0), xlab=\"log10[ p(x|M) ]\", ylim=c(0, max(c(d.ref\$y, d.sib\$y))), las=1,
             main=paste(ids[1], sprintf(\"\\nthld=%.2f err1=%.2f err2=%.2f\",p0, e1, e2)))
        lines(d.sib, col=2)
        if ( length(ids) > 5 )
        {
            legend(\"topleft\",legend=c(ids[1], \"Siblings\"), fill=c(1,2), inset=0.05, title=\"MC model (M)\", cex=1)
        } else {
            legend(\"topleft\",legend=ids, fill=c(1,rep(2,length(ids)-1)), inset=0.05, title=\"MC model (M)\", cex=0.7)
        }
        abline(v=p0, col='gray80')
    } else {
        plot(d.ref, xlim=c(min(x.sib),0), xlab=\"log10[ p(x|M) ]\", ylim=c(0, max(c(d.ref\$y, d.sib\$y))), las=1,
             main=paste(ids[1], sprintf(\"\\nthld=%.2f\",p0))) # main=\"\")
        lines(d.sib, col=2)
        if ( length(ids) > 5 )
        {
            legend(\"topleft\",legend=c(ids[1], \"Siblings\"), fill=c(1,2), inset=0.05, title=\"MC model (M)\", cex=1)
        } else {
            legend(\"topleft\",legend=ids, fill=c(1,rep(2,length(ids)-1)), inset=0.05, title=\"MC model (M)\", cex=0.7)
        }
        abline(v=p0, col='gray80')
    }

    outfile <- paste(outDir,\"/\",ids[1],\"_error.txt\", sep=\"\")
    write.table(tx, file=outfile, sep=\"\t\", row.names=F, col.names=F)
}

myHist(errTbl)
myHist(log(errTbl))
par(op)
dev.off()

thldTbl2 <- cbind(names(thldTbl),thldTbl)
colnames(thldTbl2) <- c(\"Taxon\",\"Threshold\")
dim(thldTbl2)

(outFile <- paste(outDir,\"/ncProbThlds.txt\",sep=\"\"))
write.table(thldTbl2,file=outFile, sep=\"\\t\",row.names=F, col.names=T, quote=F)

  ~;

  runRscript( $Rscript, "noErrorCheck" );
}

exit;
