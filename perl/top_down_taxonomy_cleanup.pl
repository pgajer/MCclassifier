#!/usr/bin/env perl

=head1 NAME

  top_down_taxonomy_cleanup.pl

=head1 DESCRIPTION

  Given a rooted phylogenetic tree and a lineage data of a set of sequences, the
  script iteratively curates the taxonomy of the sequences starting at the first
  high taxonomic rank where differentiation between sequences takes place and
  then successively moving down the taxonomy.

  The curation is done utilizing the vicut clustering algorithm that identifies
  internal nodes on the phylogentic tree whose induced clustering structure is
  most consistent with the taxonomy of the leaf sequences.

=head1 SYNOPSIS

  top_down_taxonomy_cleanup.pl -i <input group name>

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

=item B<--rm-ref-outliers>
  Removing sequences identified as outliers in ref_sib_pp_models.R

=item B<--show-trees>
  Show relevant for taxonomy cleanup condensed trees.

=item B<--show-lineage>
  Print lineage before any taxonomic modifications

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--quiet>
  Do not print progress messages.

=item B<--debug>
  Prints system commands

=item B<--debug-spp>
  Debug some isolated species taxonomy cleanup steps.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Firmicutes_dir

  top_down_taxonomy_cleanup.pl --show-trees --debug -i Firmicutes_group_6

  top_down_taxonomy_cleanup.pl --do-not-pop-pdfs --build-model-data --use-long-spp-names -i Firmicutes_group_6_V3V4

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Cwd qw(abs_path);
use List::Util qw( sum );
use Data::Dumper;
use File::Temp qw/ tempfile /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $maxStrLen     = 150;  # maximal number of characters in a taxon string; if
			  # larger file names using the taxon name will use an
			  # abbreviated version of the taxon name
my $taxonSizeThld = 20;
my $offsetCoef    = 0.9;
my $txSizeThld    = 10;

GetOptions(
  "input-group|i=s" 	=> \my $grPrefix,
  "rm-OGs-from-tx|r"    => \my $rmOGfromTx,
  "skip-FastTree"       => \my $skipFastTree,
  "use-long-spp-names"  => \my $useLongSppNames,
  "taxon-size-thld"     => \$taxonSizeThld,
  ##"show-all-trees"      => \my $showAllTrees,
  "do-not-pop-pdfs"     => \my $doNotPopPDFs,
  "show-trees"          => \my $showTrees,
  "show-lineage"        => \my $showLineage,
  "build-model-data"    => \my $buildModelData,
  "rm-ref-outliers"     => \my $rmRefOutliers,
  "pp-embedding"        => \my $ppEmbedding,
  "offset-coef|c=f"     => \$offsetCoef,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "igs"                 => \my $igs,
  "johanna"             => \my $johanna,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "debug-spp"           => \my $debugSpp,
  "debug2"              => \my $debug2,## file name debug
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$grPrefix )
{
  warn "\n\n\tERROR: Missing input group name";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
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

# my $showAllTreesStr = "";
# if ($showAllTrees)
# {
#   $showAllTreesStr = "--show-tree";
# }

my @suffixes;

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
  warn "\n\n\tERROR: $grDir does not exist";
  print "\n\n";
  exit 1;
}


####################################################################
##                               MAIN
####################################################################


my $section = qq~

##
## Species-level cleanup
##

~;
print "$section";

my $tmpDir = $grDir . "/temp";
my $cmd = "mkdir -p $tmpDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed: $?" if !$dryRun;

chdir $grDir;
print "--- Changed dir to $grDir\n";

my $lineageFile     = $grPrefix . ".lineage";
my $algnFile	    = $grPrefix . "_algn.fa";
my $trimmedAlgnFile = $grPrefix . "_algn_trimmed.fa";
my $outgroupFile    = $grPrefix . "_outgroup.seqIDs";
my $treeFile	    = $grPrefix . "_ginsi_rr.tree";
my $txFile          = $grPrefix . ".tx";

if ( ! -e $lineageFile )
{
  warn "\n\n\tERROR: $lineageFile does not exist";
  print "\n\n";
  exit 1;
}
elsif ( ! -e $algnFile )
{
  warn "\n\n\tERROR: $algnFile does not exist";
  print "\n\n";
  exit 1;
}
elsif ( ! -e $trimmedAlgnFile )
{
  warn "WARNING: $trimmedAlgnFile does not exist. Creating a symbolic link to $algnFile.\n";
  my $ap = abs_path( $algnFile );
  my $cmd = "rm -f $trimmedAlgnFile; ln -s $ap $trimmedAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}
elsif ( ! -e $treeFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n";
  exit 1;
}
elsif ( ! -e $outgroupFile )
{
  warn "\n\n\tERROR: $outgroupFile does not exist";
  print "\n\n";
  exit 1;
}

my $treesDir = "trees_dir";
if ( $showTrees )
{
   my $cmd = "mkdir -p $treesDir";
   print "\tcmd=$cmd\n" if $dryRun || $debug;
   system($cmd) == 0 or die "system($cmd) failed: $?" if !$dryRun;
}

## Gathering outgroup data
print "--- Parsing $outgroupFile\n";
my @ogSeqIDs = readArray($outgroupFile);
my %ogInd = map{$_ =>1} @ogSeqIDs; # outgroup elements indicator table

print "--- Testing if OG seq's form a monophylectic clade at the top or bottom of the input tree\n";
if ( test_OG($treeFile, \%ogInd) != 0 )
{
  warn "\n\n\tERROR: There is an issue with input tree outgroup seq's";
  print "\n\n";
  exit 1;
} ;

print "--- Extracting seq IDs from trimmed alignment fasta file\n";
my @seqIDs = get_seqIDs_from_fa($trimmedAlgnFile);

print "--- Parsing lineage table\n";
my %lineageTbl = readLineageTbl($lineageFile);

print "--- Checking parent consistency of the lineage table\n";
if ( check_parent_consistency(\%lineageTbl) )
{
  warn "";
  print "\n\n";
  exit 1;
}

## testing if lineage and fa files has the same seq IDs
print "--- Checking if seqIDs of $trimmedAlgnFile and $lineageFile are the same\n";
my @lSeqIDs = keys %lineageTbl;

if ( ! setequal(\@seqIDs, \@lSeqIDs) )
{
  warn "WARNING: seq IDs of trimmed alignment fasta file and lineage file do not match";
  print "Number of elements in the trimmed alignment file: " . @seqIDs . "\n";
  print "Number of elements in the lineage file: " . @lSeqIDs . "\n\n";
}

print "--- Testing if outgroup sequences are part of seqIDs\n";
my @ogDiff = diff( \@ogSeqIDs, \@seqIDs );
@ogSeqIDs = comm( \@ogSeqIDs, \@seqIDs );

if ( scalar(@ogDiff) != 0 )
{
  warn "\n\tWARNING the following outgroup seq IDs are not in the trimmed alignment file:\n\n";
  printArray(\@ogDiff);
}

if ( scalar(@ogSeqIDs) == 0 )
{
  warn "\n\n\tERROR: All outgroup seq's were lost";
  print "\n\n";
  exit 1;
}

if ($debug)
{
  print "\n\tNumber of seq's in the trimmed alignment/lineage files: " . @seqIDs . "\n";
  print   "\tNumber of outgroup seq's: " . @ogSeqIDs . "\n\n";
}

print "--- Testing if $treeFile is consistent with $trimmedAlgnFile\n";
## extracting leaves' IDs
my $treeLeavesFile = $grPrefix . "_tree.leaves";
my $cmd = "nw_labels -I $treeFile > $treeLeavesFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

## looking at the difference between leaf IDs and newTx keys
my @treeLeaves = readArray($treeLeavesFile);

if ( !setequal( \@treeLeaves, \@seqIDs ) )
{
  my @commTL = comm(\@treeLeaves, \@seqIDs);

  warn "\n\n\tERROR: seq IDs of $treeFile and $algnFile do not match";
  print "\tNumber of seq IDs of $treeFile: " . @treeLeaves . "\n";
  print "\tNumber of seq IDs of $trimmedAlgnFile: " . @seqIDs . "\n";
  print "\tNumber of common seq IDs elements: " . @commTL . "\n\n";

  exit 1;
}

if ($debug)
{
  my @commTL = comm(\@treeLeaves, \@seqIDs);

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
my %geParent;
my %faParent;
my %orParent;
my %clParent;

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
  $geParent{$ge} = $fa;
  $faParent{$fa} = $or;
  $orParent{$or} = $cl;
  $clParent{$cl} = $ph;

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

# if ($debug)
# {
#   print "\nTaxonomy BEFORE cleanup\n";
#   printLineage();
# }

# my $summaryStatsFile = "summary_stats.txt";
# open my $SRYOUT, ">$summaryStatsFile" or die "Cannot open $summaryStatsFile for writing: $OS_ERROR";
# printLineageToFile($SRYOUT, "\n\n====== Taxonomy BEFORE cleanup ======\n");

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


print  "\n\n--- Taxonomy summary\n";

print   "\tNumber of phyla:      " . scalar( keys %phTbl ) . "\n";
print   "\tNumber of classes:    " . scalar( keys %clTbl ) . "\n";
print   "\tNumber of orders:     " . scalar( keys %orTbl ) . "\n";
print   "\tNumber of families:   " . scalar( keys %faTbl ) . "\n";
print   "\tNumber of genera:     " . scalar( keys %geTbl ) . "\n";
print   "\tNumber of species:    " . scalar( keys %spTbl ) . "\n\n";

if ( $showLineage )
{
  printLineage();
}


## Taxonomic curation at the family level
if ( keys %faTbl > 1 ) # this makes sense only for at least 2 families
{
  my @seqIDs = keys %fam;

  ## Generating condensed genus tree
  for ( @ogSeqIDs )
  {
    $fam{$_} = "OG";
  }

  my @seqIDsWithOGs = keys %fam;

  my $faTxFile = $grPrefix . "_family.tx";
  writeTbl( \%fam, $faTxFile );

  ##
  ## Running vicut on family taxonomy
  ##

  ## Here an idea is that vicut clusters seq’s trying to be as consistent with
  # the given taxonomic annotation (here family taxonomy) as it can be. If there
  # are a small number of sequences that have taxonomic annotation inconsistent
  # with phylogeny at the given level, they will be clustered with sequences of
  # their correct taxonomy and so we can use this to re-annotate them
  # appropriately.

  my $vicutDir = "family_vicut_dir";

  print "--- Running vicut on the family taxonomy\n";
  $cmd = "vicut $quietStr -t $treeFile -a $faTxFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "--- Running update_tx.pl\n";
  my $vicutFaCltrFile = "$vicutDir/minNodeCut.cltrs";
  ## NOTE that the pre-vicut taxonomy table $faTxFile is passed to update_tx.pl
  ## (this is how it suppose to be :)
  $cmd = "update_tx.pl $quietStr $debugStr $useLongSppNamesStr -a $faTxFile -d $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $vicutFaTxFile =  "$vicutDir/TDupdated.tx";

  ## Its good to assess how much things improved due to vicut taxonomy update.
  print "--- Generating a vicut family condensed tree\n";
  my $vicutCondFaTreeFile = $treesDir . "/vicut_fa_cond.tree";
  condense_tree( $treeFile, $vicutFaTxFile, $vicutCondFaTreeFile );

  if ( $showTrees )
  {
    ## get family condensed tree
    print "--- Generating a family condensed tree\n";
    my $condFaTreeFile = $treesDir . "/family_cond2.tree";
    condense2_tree( $treeFile, $faTxFile, $condFaTreeFile );

    ## using nw_condense2
    my $vicutCondFaTreeFile2 = $treesDir . "/vicut_fa_cond2.tree";
    condense2_tree( $treeFile, $vicutFaTxFile, $vicutCondFaTreeFile2 );

    if ( $OSNAME eq "darwin")
    {
      ## show family condensed tree
      my $pdfCondFaTreeFile = abs_path( $treesDir . "/family_cond2_tree.pdf" );
      my $title = $grPrefix . " - family cond tree";
      plot_tree_bw( $condFaTreeFile, $pdfCondFaTreeFile, $title );
      $cmd = "open $pdfCondFaTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

      ## show vicut updated family condensed tree
      my $pdfVicutCondFaTreeFile = abs_path( $treesDir . "/vicut_fa_cond2_tree.pdf" );
      $title = $grPrefix . " - vicut family cond tree";
      plot_tree_bw( $vicutCondFaTreeFile2, $pdfVicutCondFaTreeFile, $title );
      $cmd = "open $pdfVicutCondFaTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    }
  }

  ## Its quite possible that vicut did not solve the problem of the same family
  ## appearing at different places of a condensed family tree.

  ## This can be solved by running vicut again this time making the small clusters of families appearing more than once, query sequences.

  ## Looking at the frequency of condensed family tree leaf names to identify
  ## families appearing more than once in the condensed family tree after the
  ## first run of vicut.

  print "--- Extracting leaf IDs of the _vicut_fa_cond.tree\n";
  # my $leaveslFile = $grPrefix . "_vicut_fa_cond_tree.leaves";
  my @vicutFaCondLeaves = get_leaves( $vicutCondFaTreeFile ); #, $leaveslFile );

  ## family => the number of times it appears in the condensed family tree (at the leaf level).
  my %vicutFaCondLeafFreq;
  for my $fa ( @vicutFaCondLeaves )
  {
    $vicutFaCondLeafFreq{$fa}++;
  }

  my @vicutFaCondCountGr1Leaves = grep { $vicutFaCondLeafFreq{$_} > 1 } keys %vicutFaCondLeafFreq;

  if ( $debug )
  {
    print "\nFamily frequencies (if greater than 1)\n";
    printFormatedTbl( \%vicutFaCondLeafFreq, \@vicutFaCondCountGr1Leaves );
    print "\n\n";

  }

  if ( @vicutFaCondCountGr1Leaves ) # if there are any families that occure more than once in the condensed family tree
  {
    ## for each family that occurs more than once keep the taxonomy of seq's from
    ## the larger cluster and those of size >= minTxSize and put seq IDs of
    ## sequences from other clusters to a query file

    my $minTxSize = 20;

    my %txCltrTbl = parse_cltr_tbl( $vicutFaTxFile, $vicutFaCltrFile ); # txCltrTbl{tx}{cl} ref to seq IDs of tx taxonomy and cl cluster.

    my $vicutDir  = "family_vicut2_dir";
    my $queryFile = "query2.seqIDs";
    my $annFile   = "ann2.tx";

    open QOUT, ">$queryFile" or die "Cannot open $queryFile for writing: $OS_ERROR";
    open AOUT, ">$annFile"   or die "Cannot open $annFile for writing: $OS_ERROR";
    my @query;
    my @ann;
    for my $fa ( keys %txCltrTbl )
    {
      if ( $fa eq "OG" )
      {
	for my $id ( @ogSeqIDs )
	{
	  print AOUT "$id\t$fa\n";
	  push @ann, $id;
	}
	next;
      }

      my %cltTbl = %{$txCltrTbl{$fa}};
      my @cltrs = sort { @{$cltTbl{$b}} <=> @{$cltTbl{$a}} } keys %cltTbl;

      if ( $debug )
      {
	#$Data::Dumper::Pair = " : ";     # specify hash key/value separator
	#print Dumper( $txCltrTbl{$fa} );

	print "\nfa: $fa\ncltTbl:\n";
	#printFormatedTableValuedTbl( $txCltrTbl{$fa}, \@cltrs );

	for ( @cltrs )
	{
	  print "$_: " . @{$cltTbl{$_}} . "\n";
	}
	#print "\n";
      }

      ## If the family appears in only one cluster, put all seq IDs of that
      ## family into an annotation table. If the family appears in more than one
      ## cluster, put the seq's of the largest one to the annotation table and
      ## the seq's of # the remaining clusters to the query file.

      my $cl = shift @cltrs; # Note, that the clsuters are sorted by size (see above).
      my @ids = @{$cltTbl{$cl}};
      print "Largest cl: $cl\nn(ids): " . @ids . "\n" if $debug;
      for my $id (@ids)
      {
	print AOUT "$id\t$fa\n";
	push @ann, $id;
      }

      for my $cl ( @cltrs ) # if @cltrs is non-empty, we are dealing with a
	                    # family that apprears in at least two clusters put
      {
	my @ids = @{$cltTbl{$cl}};
	if ( @ids < $minTxSize )
	{
	  print "Adding $cl to Query\n" if $debug;
	  for my $id (@ids)
	  {
	    print QOUT "$id\n";
	    push @query, $id;
	  }
	}
	else
	{
	  print "Adding $cl to Annotation\n" if $debug;
	  for my $id (@ids)
	  {
	    print AOUT "$id\t$fa\n";
	    push @ann, $id;
	  }
	}
      }
    } # end of   for my $fa ( keys %txCltrTbl )
    close QOUT;
    close AOUT;

    if ( $debug )
    {
      my @all = (@ann, @query);

      print "\n\nsize(ann):   " . @ann . "\n";
      print "size(query): " . @query . "\n";
      print "size(all):   " . @all . "\n";

      my @d = diff( \@seqIDs, \@all );
      print "d: @d\n\n";
    }

    ## vicut run 2
    print "--- Vicut run (2) on the family taxonomy\n";
    $vicutDir = "family_vicut2_dir";
    if ( @query )
    {
      $cmd = "vicut $quietStr -t $treeFile -a $annFile -q $queryFile -o $vicutDir";
    }
    else
    {
      $cmd = "vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
    }

    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    print "--- Running family update_tx.pl (2)\n";
    ## NOTE: passing pre 2 run of vicut taxonomy $vicutFaTxFile
    $cmd = "update_tx.pl $quietStr $debugStr $useLongSppNamesStr -a $vicutFaTxFile -d $vicutDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    my $vicutFaTx3File = "$vicutDir/TDupdated2.tx";
    my %faTx = read_tx3_tbl($vicutFaTx3File);

    $vicutFaTxFile = "$vicutDir/TDupdatedIdx.tx";
    writeTbl(\%faTx, $vicutFaTxFile);


    if ( $showTrees )
    {
      print "--- Generating a vicut(2) family condensed tree\n";
      my $vicutCondFaTreeFile = $treesDir . "/vicut2_fa_cond2.tree";
      condense2_tree( $treeFile, $vicutFaTxFile, $vicutCondFaTreeFile );

      ## show vicut updated family condensed tree
      my $pdfVicutCondFaTreeFile = abs_path( $treesDir . "/vicut2_fa_cond2_tree.pdf" );
      my $title = $grPrefix . " - vicut (2) family cond tree";
      plot_tree_bw( $vicutCondFaTreeFile, $pdfVicutCondFaTreeFile, $title );
      $cmd = "open $pdfVicutCondFaTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    }
  }

  exit;
}


## Taxonomic curation at the genus level


## Generating condensed genus tree
for (@ogSeqIDs)
{
  $gen{$_} = "OG";
}
my $geTxFile = $grPrefix . "_genus.tx";
writeTbl(\%gen, $geTxFile);

print "--- Generating a tree with genus names at leaves\n";
my $geTreeFile = $grPrefix . "_genus.tree";
$cmd = "rm -f $geTreeFile; nw_rename $treeFile $geTxFile | nw_order -c n  - > $geTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Generating a genus condensed tree\n";
my $condGeTreeFile = $grPrefix . "_genus_cond.tree";
$cmd = "rm -f $condGeTreeFile; nw_condense2 $geTreeFile > $condGeTreeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

if ( $showTrees &&  $OSNAME eq "darwin")
{
  my $pdfCondGeTreeFile = abs_path( "$grPrefix" . "_genus_cond_tree.pdf" );
  plot_tree_bw($condGeTreeFile, $pdfCondGeTreeFile);

  $cmd = "open $pdfCondGeTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}




####################################################################
##                               SUBS
####################################################################

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

sub plot_color_tree
{
  my ($treeFile, $pdfFile, $title) = @_;

  my $showBoostrapVals = "T";

  if (!defined $title)
  {
    $title = "";
  }

  my $Rscript = qq~

source(\"$readNewickFile\")
require(phytools)
library(ade4)

treeStr <- readChar(\"$treeFile\", file.info(\"$treeFile\")\$size)
tr1 <- newick2phylog(treeStr)

parent <- c()
for ( node in names(tr1\$parts))
{
    chld <- tr1\$parts[[node]]
    for ( ch in chld )
    {
        parent[ch] <- node
    }
}

tr <- read.newick(file=\"$treeFile\")
tr <- collapse.singles(tr)

pars <- parent[tr\$tip.label]

colIdx <- 1
tip.colors <- c()
tip.colors[1] <- colIdx
for ( i in 2:length(pars) )
{
    if ( pars[i] != pars[i-1] )
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

(nLeaves <- length(tr\$tip.label))

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr, tip.color=tip.colors, type=\"phylogram\", no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8, main=\"$title\")
par(op)
dev.off()
~;

  runRscript( $Rscript );
}

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

  my ($fh, $inFile) = tempfile("rTmpXXXX", SUFFIX => '.R', OPEN => 1, DIR => $tmpDir);
  print $fh "$Rscript";
  close $fh;

  my $outFile = $inFile . "out";
  my $cmd = "R CMD BATCH $inFile $outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

  if (!$noErrorCheck)
  {
    open IN, "$outFile" or die "Cannot open $outFile for reading: $OS_ERROR";
    my $exitStatus = 1;
    foreach my $line (<IN>)
    {
      if ( $line =~ /Error/ )
      {
	print "R script crashed at\n$line";
	print "check $outFile for details\n";
	$exitStatus = 0;
	exit 1;
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

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR";
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


# print elements of a hash table whose values are reference to a hash table so
sub printTableValuedTbl{

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

  for (@args)
  {
    print "$_\n";
    for my $e ( keys %{$rTbl->{$_}} )
    {
      print "\t$e\n";
    }
  }
  print "\n";
}


# print elements of a hash table whose values are reference to an array so
# that indents of the value array elements are the same
sub printFormatedArrayValuedTbl{

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

  for (@args)
  {
    print "$_\n";
    for my $e (@{$rTbl->{$_}})
    {
      print "\t$e\n";
    }
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
    warn "\n\n\tERROR in readArray(): $file does not exist";
    print "\n\n";
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

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
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



#
# read 3 column table
#
# format
#
# S002234001	f_Staphylococcaceae	1
# S002351864	f_Staphylococcaceae	1
# S002965910	f_Staphylococcaceae	1
# S002098342	f_Staphylococcaceae	1
#
sub read_tx3_tbl
{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $tx, $idx) = split /\s+/,$_;
    if ( $idx > 0 )
    {
      $tbl{$id} = $tx . "_$idx";
    }
    else
    {
      $tbl{$id} = $tx;
    }
  }
  close IN;

  return %tbl;
}


# read taxon table excluding outgroup seq ID and any seqID with OUTGROUP.* type label
sub readTxTbl{

  my ($file, $outgroupSeqID) = @_;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTxTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
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
    open OUT, ">$file" or die "Cannot open file $file to write: $!";
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
    warn "\n\n\tERROR in findSeqIDinFasta(): $inFile and does not exist";
    print "\n\n";
    exit 1;
  }

  my $found = 0;
  open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR";
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
    warn "\n\n\tERROR in readBigTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  my @header;
  my @rowIds;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";

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

# read lineage table
sub readLineageTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readLineageTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
    ## test for '/' characters
    if ($t =~ /\//)
    {
      warn "\n\n\tERROR: Discovered '/' for id: $id\t$t";
      print "\n\n";
      exit 1;
    }
  }
  close IN;

  return %tbl;
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
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


The top_down_taxonomy_cleanup.log file contains the standard output of the taxonomic cleanup command.

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

  open OUT, ">$file" or die "Cannot open $file for writing: $OS_ERROR";
  print OUT "$text\n";
  close OUT;
}


# test if OG seq's form one or more clusters in the tree
sub test_OG
{
  my ($treeFile, $rogInd) = @_;

  my %ogInd = %{$rogInd};

  my $debug_test_OG = 0;

  my $ret = 0;

  print "\t--- Extracting leaves\n" if $debug_test_OG;
  my $treeLeavesFile = "$grPrefix" . "_sppSeqIDs.leaves";
  my $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

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
    warn "$grPrefix\n\n\tERROR: start and end arrays have different lengths!";
    print "length(start): " . @start . "\n";
    print "length(end): " . @end . "\n\n";
    $ret = 1;
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
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$grPrefix" . "_clade.leaves";
	$cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

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
    $ret = 1;
  }
  elsif ( !( $start[0] == 0 || $end[0] == $#leaves) )
  {
    warn "\n\n\tERROR: In the pruned tree outgroups sequences are not at the top or bottom of the tree!";

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
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

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

    $ret = 1;
  }
  else
  {
    print "\t--- Extracting the clade of OG sequences\n" if $debug_test_OG;
    my $ogCladeTreeFile = "$grPrefix" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

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
      $ret = 1;
    }
  }

  print "\n\tOG seq's form a monophylectic clade at the top or bottom of the input tree\n\n" if $debug_test_OG;

  return $ret;
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
  my ($mcDir, $errorDir, $refSibDir) = @_;

  my $Rscript = qq~

myHist <- function(x, main=\"\", br=100, ... )
{
    hist(x, br=br, col=2, las=1, main=main, ...)
}

## integrate d from x0 to x1
myIntegral <- function(d, x0, x1){
    dx <- d\$x[2]-d\$x[1]
    i0 <- which.min(abs(d\$x - x0))
    i1 <- which.min(abs(d\$x - x1))
    sum(d\$y[i0:i1]) * dx
}

file <- paste(\"$refSibDir\",\"/ref.postProbs\",sep=\"\")
ref.nFields <- count.fields(file, sep = \"\\t\")
ref.maxNumFields <- max(ref.nFields)
ref.pp <- read.table(file, sep=\"\\t\", col.names = paste0("V",seq_len(ref.maxNumFields)), fill = TRUE, stringsAsFactors=FALSE, row.names=1, header=F)
ref.nFields <- ref.nFields - 1
names(ref.nFields) <- rownames(ref.pp)

file <- paste(\"$refSibDir\",\"/sib.postProbs\",sep=\"\")
sib.nFields <- count.fields(file, sep = \"\\t\")
sib.maxNumFields <- max(sib.nFields)
sib.pp <- read.table(file, sep=\"\\t\", col.names = paste0("V",seq_len(sib.maxNumFields)), fill = TRUE, stringsAsFactors=FALSE, row.names=1, header=F)
sib.nFields <- sib.nFields - 1
names(sib.nFields) <- rownames(sib.pp)

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
fpError <- c()
fnError <- c()

outDir <- \"$mcDir\"
(figFile <- paste(outDir,\"/error_thld_figs.pdf\",sep=\"\"))
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

    tbl <- read.table(file, sep=\"\\t\", col.names = paste0(\"V\",seq_len(maxNumFields)), fill = TRUE, stringsAsFactors=FALSE, row.names=1)
    dim(tbl)
    ids <- rownames(tbl)
    refID <- ids[1]

    nFields <- nFields - 1 # the first column is row name

    i <- which(ids==refID)
    x.ref <- as.numeric(as.matrix(tbl[i,1:nFields[1]]))

    ii <- setdiff(seq(ids), i)
    x.sib <- c()
    for ( i in ii )
    {
        x.sib <- c(x.sib, as.numeric(as.matrix( tbl[i,1:nFields[i]] )))
    }

    if ( length(x.ref)==0 ) {
        stop(\"length(x.ref) is 0\")
    } else if ( length(x.sib)==0 ) {
        stop(\"length(x.sib) is 0\")
    }

    d.ref <- density(x.ref)
    d.sib <- density(x.sib)

    xmin.ref <- min(x.ref)
    xmax.sib <- max(x.sib)

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

    # plot(x,y,main=refID, type='l')
    # abline(h=0, col='gray80'); abline(v=p0,col='gray80')
    # abline(v=x[which.min(y)], col='red')
    # abline(v=x[which.max(y)], col='red')

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

    if ( xmax.sib < xmin.ref )
    {
      ##p0 <- min( c( (xmax.sib + xmin.ref) / 2, xmin.ref ) )
      p0 <- min( c( (xmax.sib + xmin.ref) / 2, xmax.sib ) )
    }

    thldTbl[refID] <- p0

    ## probability of a FP error = integral of d.sib from p0 to +inf
    fpError[refID] <- 0
    fnError[refID] <- 0
    if ( p0 < max(d.sib\$x) ) {
      fpError[refID] <- myIntegral(d.sib, p0, max(d.sib\$x)) # integrate(d.sib.fun, lower=p0, upper=max(d.sib\$x))[[1]]
      fnError[refID] <- myIntegral(d.ref, min(d.ref\$x), p0) # integrate(d.ref.fun, lower=min(d.ref\$x), upper=p0, subdivisions=2000)[[1]]
    }

    errTbl[refID] <- fpError[refID] + fnError[refID]

    plot(d.ref, xlim=c(min(x.sib),0), xlab=\"log10[ p(x|M) ]\", ylim=c(0, max(c(d.ref\$y, d.sib\$y))), las=1,
             main=paste(refID, " (n.ref=", ref.nFields[refID][[1]], ")", sprintf(\"\\nthld=%.2f  fpError=%.2f  fnError=%.2f\",p0, fpError[refID], fnError[refID])))
    lines(d.sib, col=2)
    if ( length(ids) > 5 ) {
      legend(\"topleft\",legend=c(refID, \"Siblings\"), fill=c(1,2), inset=0.05, title=\"MC model (M)\", cex=1)
    } else {
      legend(\"topleft\",legend=ids, fill=c(1,rep(2,length(ids)-1)), inset=0.05, title=\"MC model (M)\", cex=0.7)
    }
    abline(v=p0, col='gray80')
    hist(as.numeric((ref.pp[refID, 1:ref.nFields[refID]])), add=T, col=1, br=100)
    hist(as.numeric((sib.pp[refID, 1:sib.nFields[refID]])), add=T, col=2, br=100)

    outfile <- paste(outDir,\"/\",refID,\"_error.txt\", sep=\"\")
    write.table(tx, file=outfile, sep=\"\t\", row.names=F, col.names=F)
}
par(op)
dev.off()

thldTbl2 <- cbind(names(thldTbl),thldTbl)
colnames(thldTbl2) <- c(\"Taxon\",\"Threshold\")
dim(thldTbl2)

(outFile <- paste(outDir,\"/ncProbThlds.txt\",sep=\"\"))
write.table(thldTbl2,file=outFile, sep=\"\\t\",row.names=F, col.names=T, quote=F)

(outFile <- paste(outDir,\"/errorTbl.txt\",sep=\"\"))
write.table(cbind(errTbl),file=outFile, sep=\"\\t\",row.names=F, col.names=T, quote=F)

(outFile <- paste(outDir,\"/fpErrorTbl.txt\",sep=\"\"))
write.table(cbind(fpError),file=outFile, sep=\"\\t\",row.names=F, col.names=T, quote=F)

(outFile <- paste(outDir,\"/fnErrorTbl.txt\",sep=\"\"))
write.table(cbind(fnError),file=outFile, sep=\"\\t\",row.names=F, col.names=T, quote=F)

(figFile <- paste(outDir,\"/error_tbl_hists.pdf\",sep=\"\"))
pdf(figFile, width=11, height=11)
op <- par(mar=c(4,4,4.5,0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
myHist(errTbl)
myHist(log(errTbl))
par(op)
dev.off()

  ~;

  runRscript( $Rscript, "noErrorCheck" );
}

## check if each node of the lineage structure has only one parent
sub check_parent_consistency
{
  my $r = shift;
  my %lineageTbl = %{$r};

  my  %prt;
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

    $subGe = "sg_$subGe";
    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $prt{$sp}{$subGe}++;
    $prt{$subGe}{$ge}++;
    $prt{$ge}{$fa}++;
    $prt{$fa}{$or}++;
    $prt{$or}{$cl}++;
    $prt{$cl}{$ph}++;
  }

  my $ret = 0;
  for my $tx (keys %prt)
  {
    my $nPrts = keys %{$prt{$tx}};
    if ( $nPrts > 1 )
    {
      $ret = 1;
      print "\n\n\tERROR: $tx has more than one parent";
      print "\n\t$tx parents\n";
      for (keys %{$prt{$tx}})
      {
	print "\t\t$_\n";
      }
      print "\n\n";
    }
  }

  return $ret;
}

## are two arrays equal set-theoretically
sub setequal
{
  my ($rA, $rB) = @_;

  my @a = @{$rA};
  my @b = @{$rB};
  my @c = comm(\@a, \@b);

  my $ret = 1;

  if (@c != @a || @c != @b)
  {
    warn "\n\n\tERROR: Elements of the two arrays do not match";
    print "\n\tNumber of elements in the first array: " . @a . "\n";
    print "\tNumber of elements in the second array: " . @b . "\n";
    print "\tNumber of common elements: " . @c . "\n";

    print "\na: ";
    map {print "$_ "} @a;
    print "\n\n";

    print "b: ";
    map {print "$_ "} @b;
    print "\n\n";

    # writeArray(\@a, "a.txt");
    # writeArray(\@b, "b.txt");
    #print "\n\tNew taxon keys and fasta IDs written to a.txt and b.txt, respectively\n\n";

    if (@a > @b)
    {
      my @d = diff(\@a, \@b);
      print "\nElements in a, but not b:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    if (@b > @a)
    {
      my @d = diff(\@b, \@a);
      print "\nElements in b that are not in a:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    $ret = 0;
  }

  return $ret;
}

sub build_cond_spp_tree
{
  my ($treeFile, $txFile, $condSppFile) = @_;

  my %txTbl = readTbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  my $treeLeavesFile = $grPrefix . "_tmp_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @treeLeaves = readArray($treeLeavesFile);

  if ( !setequal( \@txs, \@treeLeaves ) )
  {
    warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
    print "\n\n";
    exit 1;
  }

  my $sppTreeFile = "$grPrefix" . "_tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; nw_rename $treeFile $txFile | nw_order -c n  - > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $cmd = "rm -f $condSppFile; nw_condense $sppTreeFile > $condSppFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\n\ntree:        $treeFile\n";
  print     "tx:          $txFile\n";
  print     "sppTree:     $sppTreeFile\n";
  print     "condSppTree: $condSppFile\n\n";

}

## For each species the script looks at

## 1. the clade of that species in the spp.tree
## 2. the clade of that species in the spp_cond.tree

## A purity of a clade is the entropy of the proportions of different species
## there. Purity 0 means, the clade consists of only one species.

## The sum of species clade purities is reported for the spp and spp_cond trees.

## Clades are also reported
sub get_tree_spp_purity
{
  my ($treeFile, $txFile, $condSppTreeFile) = @_;

  my %txTbl = readTbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  my $treeLeavesFile = $grPrefix . "_tmp_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @treeLeaves = readArray($treeLeavesFile);

  if ( !setequal( \@txs, \@treeLeaves ) )
  {
    warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
    print "\n\n";
    exit 1;
  }

  my $sppTreeFile = "$grPrefix" . "_tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; nw_rename $treeFile $txFile | nw_order -c n  - > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  # $cmd = "rm -f $treeLeavesFile; nw_labels -I $sppTreeFile > $treeLeavesFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  # my @sppTreeLeaves = readArray($treeLeavesFile);

  $cmd = "rm -f $condSppTreeFile; nw_condense $sppTreeFile > $condSppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $cmd = "rm -f $treeLeavesFile; nw_labels -I $condSppTreeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @condSppTreeLeaves = readArray($treeLeavesFile);

  my %spFreq; # species frequency
  map { $spFreq{$_}++ } @condSppTreeLeaves;

  # selecting species appearing more than once in the condensed species tree
  my @hiFreqSpp;
  for (keys %spFreq)
  {
    if ( $spFreq{$_} > 1 )
    {
      push @hiFreqSpp, $_;
    }
  }

  @hiFreqSpp = sort { $spFreq{$b} <=> $spFreq{$a} } @hiFreqSpp;

  print "\n\ntree:        $treeFile\n";
  print     "tx:          $txFile\n";
  print     "sppTree:     $sppTreeFile\n";
  print     "condSppTree: $condSppTreeFile\n\n";

  print "\nSpecies freq's within the species condensed tree\n";
  printFormatedTbl(\%spFreq, \@hiFreqSpp);
  print "\n\n";

  ## looking at the clades of of hiFreqSpp

  ## we need to know seq IDs of seq's of a given species
  #
  my %spToSeqIDs;
  for my $id (keys %txTbl)
  {
    push @{$spToSeqIDs{ $txTbl{$id} }}, $id;
  }

  print "\n\nhiFreqSpp\n";
  for my $sp (@hiFreqSpp)
  {
    my @ids = @{$spToSeqIDs{ $sp }};

    print "$sp (n=" . @ids . ")\t";

    ## Generating species' clade tree
    my $cladeFile = $sp . "_clade.tree";
    my $cmd = "rm -f $cladeFile; nw_clade $treeFile @ids | nw_order -c n  - > $cladeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    ## Getting leaves
    my $lFile = "leaves.txt";
    $cmd = "rm -f $lFile; nw_labels -I $cladeFile > $lFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    my @leaves = readArray($lFile);

    ## entropy of the species
    my @lSpp = @txTbl{@leaves}; # leaf species
    my %lSpFreqTbl; # frequency table of leaf species
    map { $lSpFreqTbl{$_}++ } @lSpp;

    my @lSpFreq = values %lSpFreqTbl;
    # my $total = sum( @lSpFreq );
    # my @lSpProp = map{ $_ / $total } @lSpFreq;
    my $H = evenness( \@lSpFreq );

    print "n(clade spp)=" . (keys %lSpFreqTbl) . "\tEvenness: " . sprintf("%.4f", $H) . "\n";
    #  print "\nSpecies freq's within the species' clade\n";
    my @clSpp = sort { $lSpFreqTbl{$b} <=> $lSpFreqTbl{$a} } keys  %lSpFreqTbl;
    printFormatedTbl( \%lSpFreqTbl, \@clSpp );
    print "\n\n";

    ## leaf => sp table
    my $l2spFile = $sp . "_leaf_to_spp.txt";
    open OUT, ">$l2spFile" or die "Cannot open $l2spFile for writing: $OS_ERROR";
    for (@leaves)
    {
      print OUT "$_\t" . $txTbl{$_} . "\n";
    }
    close OUT;
    ## species tree
    my $sFile = $sp . "_clade_spp.tree";
    $cmd = "rm -f $sFile; nw_rename $cladeFile $l2spFile | nw_order -c n  - > $sFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    ## condensed species tree
    my $csFile = $sp . "_clade_cond_spp.tree";
    $cmd = "rm -f $csFile; nw_condense $sFile > $csFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    if ( $showTrees && $OSNAME eq "darwin")
    {
      ## plot the condensed species tree
      my $pdfFile = $treesDir . $sp . "_clade_cond_spp_tree.pdf";
      my $title   = $sp . " clade cond spp tree";

      # $cmd = "root_tree.pl -i $csFile";
      # print "\tcmd=$cmd\n" if $dryRun || $debug;
      # system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

      ##plot_color_tree($csFile, $pdfFile, $title);
      plot_tree_bw($csFile, $pdfFile, $title);

      $cmd = "open $pdfFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    }


  }
}

# Shannon entropy
sub entropy
{
    my $x = shift;
    my $n = sum (@{$x});
    #my @p = map{ $_ / $n } @{$x};

    my $H = 0; # output value
    foreach my $v ( @{$x} )
    {
      if ( $v )
      {
	my $p = $v / $n;
	$H -= $p * log( $p );
      }
    }
    return $H;
}


sub evenness
{
    my $x = shift;
    my $n = sum (@{$x});
    #my @p = map{ $_ / $n } @{$x};

    my $H = 0; # output value
    foreach my $v ( @{$x} )
    {
      if ( $v )
      {
	my $p = $v / $n;
	$H -= $p * log( $p );
      }
    }
    return $H / log( @{$x} );
}

# Parse vicut's Cltrs table
sub parse_cltr_tbl
{
  my ( $txFile, $cltrFile ) = @_;

  if ( $debug )
  {
    print "\nParsing\n";
    print "$txFile\n";
    print "$cltrFile\n\n";
  }

  my %txTbl = readTbl( $txFile );

  my %txCltrTbl; # txCltrTbl{tx}{cl} ref to seq IDs of tx taxonomy and cl cluster.
  open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR";
  my $header = <IN>;
  my $i = 1;
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    $cl = "c$cl";
    ##print "$i: $id  $cl  $tx\n";
    $i++;
    #push @{$txCltrTbl{$txTbl{$id}}{$cl}}, $id;
    push @{$txCltrTbl{$tx}{$cl}}, $id;
  }
  print "i: $i\n\n" if $debug;
  close IN;

  return %txCltrTbl;
}

## Create a condensed tree given a tree, translation table and a file name of the
## condensed tree
sub condense_tree
{
  my ($treeFile, $txFile, $condTreeFile) = @_;

  my %txTbl = readTbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  my $treeLeavesFile = "tmp_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @treeLeaves = readArray ( $treeLeavesFile );
  unlink ( $treeLeavesFile );

  if ( !setequal( \@txs, \@treeLeaves ) )
  {
    warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
    print "\n\n";
    exit 1;
  }

  my $sppTreeFile = "tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; nw_rename $treeFile $txFile | nw_order -c n  - > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $cmd = "rm -f $condTreeFile; nw_condense $sppTreeFile > $condTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ( $debug )
  {
    print "\n\ntree:        $treeFile\n";
    print     "tx:          $txFile\n";
    print     "sppTree:     $sppTreeFile\n";
    print     "condTree: $condTreeFile\n\n";
  }
}

## Version of condense_tree using nw_condense2
## condensed tree
sub condense2_tree
{
  my ($treeFile, $txFile, $condTreeFile) = @_;

  my %txTbl = readTbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  my $treeLeavesFile = "tmp_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @treeLeaves = readArray ( $treeLeavesFile );
  unlink ( $treeLeavesFile );

  if ( !setequal( \@txs, \@treeLeaves ) )
  {
    warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
    print "\n\n";
    exit 1;
  }

  my $sppTreeFile = "tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; nw_rename $treeFile $txFile | nw_order -c n  - > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $cmd = "rm -f $condTreeFile; nw_condense2 $sppTreeFile > $condTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ( $debug )
  {
    print "\n\ntree:        $treeFile\n";
    print     "tx:          $txFile\n";
    print     "sppTree:     $sppTreeFile\n";
    print     "condTree: $condTreeFile\n\n";
  }
}

sub get_leaves
{
  my ($treeFile, $leavesFile) = @_;

  if ( !defined $leavesFile )
  {
    $leavesFile = "tmp.leaves";
  }
  my $cmd = "rm -f $leavesFile; nw_labels -I $treeFile > $leavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @a = readArray($leavesFile);

  return @a;
}

exit 0;
