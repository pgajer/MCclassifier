#!/usr/bin/env perl

=head1 NAME

  vicut_classify.pl

=head1 DESCRIPTION

  Given an annotation table and a set of query sequences and a tree with leaves
  being the union of query and annotation sequences, run vicut on the data and
  use the majority vote to assign annotation to query sequences.

  For each vicut cluster, generate a grandparent tree with seqID_annot and '_q'
  suffix for the quenry sequences that had their annotation label assigned in
  the above process. Unclassified query sequences would have seqID_query labels.

=head1 SYNOPSIS

  vicut_classify.pl -t <tree file> -a <ann file> -q <query file> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--li-file, -l>
  Master lineage file (seqID => lineage)

=item B<--ann-file, -a>
  Annotation file (seqID => labels)

=item B<--tree-file, -t>
  Consensus tree file

=item B<--label, -b>
  Label to use for tree files.

=item B<--output-dir, -o>
  Output dir with the following files.

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

  cd /Users/pgajer/projects/PECAN/data/Banfield_contax/FL

  vicut_classify.pl -b phylum -t cxhb_FL_nr.tree -a cx_hb_phylum.tx -o cxhb_FL_nr_phylum_dir

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd qw(abs_path);
use File::Temp qw/ tempfile /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "tree-file|t=s"   => \my $treeFile,
  "ann-file|a=s"    => \my $annFile,
  "query-file|q=s"  => \my $queryFile,
  "output-dir|o=s"  => \my $outDir,
  "label|b=s"       => \my $label,
  "verbose|v"       => \my $verbose,
  "quiet"           => \my $quiet,
  "debug"           => \my $debug,
  "dry-run"         => \my $dryRun,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$outDir )
{
  print "ERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$annFile )
{
  print "ERROR: Missing annotation file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$treeFile )
{
  print "ERROR: Missing tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}


if ( ! -e $annFile || ! -s $annFile )
{
  warn "ERROR: Annonation file $annFile does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

if ( ! -e $treeFile || ! -s $treeFile )
{
  warn "ERROR: Tree file $treeFile does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

# creating output directory
my $cmd = "mkdir -p $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $tmpDir = $outDir . "/temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my $treeDir = $outDir . "/trees_dir";

if ( ! -e $treeDir )
{
  my $cmd = "mkdir -p $treeDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my $quietStr = "";
if ( $quiet )
{
  $quietStr = "--quiet";
}

my $debugStr = "";
if ( $debug )
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ( $verbose )
{
  $verboseStr = "--verbose";
}

if ( !$label )
{
  $label = "ann";
}


####################################################################
##                               MAIN
####################################################################

print "--- Computing frequencies of annotation labels\n";
my %ann = read_tbl( $annFile );
my %annFreq;
map { $annFreq{$_}++ } values %ann;
my @uqAnn = keys %annFreq;

print "--- Running vicut\n";
my $vicutDir = $outDir . "/vicut_dir";
my @query = run_vicut( $treeFile, $annFile, $vicutDir );

# my $cmd = "vicut -t $treeFile -a $annFile -q $queryFile -o $vicutDir";
# print "\tcmd=$cmd\n" if $dryRun || $debug;
# system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Parsing results of vicut run\n";
my %clTxTbl = parse_cltr_tbl(); # clTxTbl{cl}{tx} ref to seq IDs of tx taxonomy within the cl cluster.

print "--- Majority vote classification of query elements\n";

my $clFile = $outDir . "/$label" . "_classification_results.txt";
open OUT, ">$clFile" or die "Cannot open $clFile for writing: $OS_ERROR\n";

my $chgFile = $outDir . "/$label" . "_taxonomy_changes.txt";
open OUT2, ">$chgFile" or die "Cannot open $chgFile for writing: $OS_ERROR\n";

my %txChg;    # table recording changes of annotation of annotated elements
my %naChg;    # $naChg{$cl} seq IDs of $cl; only clusters containing NA elements are in this table; it will be used to generate grandparent trees
my %newAnn;   # new taxonomy
my %txClTbl;  # taxon to cluster to seqIDs table
for my $cl ( keys %clTxTbl )
{
  my @txs = sort { @{$clTxTbl{$cl}{$b}} <=> @{$clTxTbl{$cl}{$a}} } keys %{$clTxTbl{$cl}};
  my $moAnn = shift @txs; # most abundant annotation
  if ( exists $clTxTbl{$cl}{"NA"} )
  {
    if ( $moAnn eq "NA" && @txs > 0 )
    {
      $moAnn = shift @txs;
    }
    elsif ( $moAnn eq "NA" && @txs == 0 )
    {
      $moAnn = $cl;
    }
  }

  for my $tx ( keys %{$clTxTbl{$cl}} )
  {
    my @ids = @{$clTxTbl{$cl}{$tx}};
    if ( $tx eq "NA" )
    {
      push @{$naChg{$cl}}, @ids;
    }

    if ( $tx ne "NA" && $tx ne $moAnn )
    {
      $txChg{$tx} = $moAnn;
      for my $id ( @ids )
      {
        print OUT2 "$id\t$tx\t$moAnn\n";
      }
    }

    push @{$txClTbl{$moAnn}{$cl}}, @ids;

    for my $id ( @ids )
    {
      print OUT "$id\t$moAnn\n";
      $newAnn{$id} = $moAnn;
    }
  }
}
close OUT;
close OUT2;

print "--- Moving genera present in more than one cluster to query in the next run of vicut\n";
for my $tx ( keys %txClTbl )
{
  my @cls = sort { @{$txClTbl{$tx}{$b}} <=> @{$txClTbl{$tx}{$a}} } keys %{$txClTbl{$tx}};
  if ( @cls > 1 )
  {
    my $winner = shift @cls;
    for my $cl ( @cls )
    {
      my @ids = @{$txClTbl{$tx}{$cl}};
      delete @ann{@ids};
    }
  }
}


$annFile = $outDir . "/$label" . "_ann.txt";
write_tbl( \%ann, $annFile );


print "--- Running 2nd vicut\n";
$vicutDir = $outDir . "/vicut_dir2";
@query = run_vicut( $treeFile, $annFile, $vicutDir );

print "--- Parsing results of vicut run\n";
%clTxTbl = parse_cltr_tbl(); # clTxTbl{cl}{tx} ref to seq IDs of tx taxonomy within the cl cluster.

print "--- Majority vote classification of query elements\n";

$clFile = $outDir . "/$label" . "_classification_results2.txt";
open OUT, ">$clFile" or die "Cannot open $clFile for writing: $OS_ERROR\n";

$chgFile = $outDir . "/$label" . "_taxonomy_changes2.txt";
open OUT2, ">$chgFile" or die "Cannot open $chgFile for writing: $OS_ERROR\n";

undef %txChg;    # table recording changes of annotation of annotated elements
undef %naChg;    # $naChg{$cl} seq IDs of $cl; only clusters containing NA elements are in this table; it will be used to generate grandparent trees
undef %newAnn;   # new taxonomy
for my $cl ( keys %clTxTbl )
{
  my @txs = sort { @{$clTxTbl{$cl}{$b}} <=> @{$clTxTbl{$cl}{$a}} } keys %{$clTxTbl{$cl}};
  my $moAnn = shift @txs; # most abundant annotation
  if ( exists $clTxTbl{$cl}{"NA"} )
  {
    if ( $moAnn eq "NA" && @txs > 0 )
    {
      $moAnn = shift @txs;
    }
    elsif ( $moAnn eq "NA" && @txs == 0 )
    {
      $moAnn = $cl;
    }
  }

  for my $tx ( keys %{$clTxTbl{$cl}} )
  {
    my @ids = @{$clTxTbl{$cl}{$tx}};
    if ( $tx eq "NA" )
    {
      push @{$naChg{$cl}}, @ids;
    }

    if ( $tx ne "NA" && $tx ne $moAnn )
    {
      $txChg{$tx} = $moAnn;
      for my $id ( @ids )
      {
        print OUT2 "$id\t$tx\t$moAnn\n";
      }
    }

    for my $id ( @ids )
    {
      print OUT "$id\t$moAnn\n";
      $newAnn{$id} = $moAnn;
    }
  }
}
close OUT;
close OUT2;


print "--- Generating tree with vicut classification derived labels\n";
my $annTreeFile = $treeDir . "/$label" . ".tree";
build_ann_tree( $treeFile, $clFile, $annTreeFile );

print "--- Generating condensed tree\n";
my $condAnnTreeFile = $treeDir . "/condensed_$label" . ".tree";
condense_tree_only( $annTreeFile, $condAnnTreeFile );

my $condAnnTreeFile2 = $treeDir . "/condensed2_$label" . ".tree";
condense2_tree_only( $annTreeFile, $condAnnTreeFile2 );

print "--- Generating pdf of the condensed tree\n";
my $pdfTreeFile = $treeDir . "/condensed2_$label" . "_tree.pdf";
plot_tree( abs_path( $condAnnTreeFile2 ), "", $pdfTreeFile );

if ( $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

print "--- Computing frequencies of new annotation labels\n";
my %newAnnFreq;
map { $newAnnFreq{$_}++ } values %newAnn;
my @uqNewAnn = keys %newAnnFreq;

# print "\n\n newAnnFreq\n";
# print_formated_tbl( \%newAnnFreq );

print "--- Generating tree with seqID + classification derived labels\n";
## Adding _q to query sequences
my %qInd = map { $_ => 1 } @query;
for my $id ( keys %newAnn )
{
  if ( !exists $newAnn{$id} )
  {
    warn "\n\n\tERROR: $id not defined in newAnn";
    print "\n\n";
    exit;
  }
  if ( exists $qInd{$id} )
  {
    $newAnn{$id} = $id . "_" . $newAnn{$id} . "_q";
  }
  else
  {
    $newAnn{$id} = $id . "_" . $newAnn{$id};
  }
}

my $labsFile = $treeDir . "/seqID_$label" . "_q.txt";
write_tbl( \%newAnn, $labsFile );

my $annTreeFile2 = $treeDir . "/seqID_$label" . ".tree";
build_ann_tree( $treeFile, $labsFile, $annTreeFile2 );

print "--- Generating pdf figure of the tree with seqID + classification derived labels\n";
$pdfTreeFile = $treeDir . "/seqID_$label" . "_tree.pdf";
plot_tree( abs_path( $annTreeFile2 ), "", $pdfTreeFile );

if ( 0 && $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

##print "--- Generating grandparent trees for clusters with query sequences\n";
# for my $cl ( keys %naChg )
# {
#     my @ids = keys @{$naChg{$cl}};
#     my $clTreeFile = $treeDir . "/$cl" . "_ann.tree";
#     my $title = "$cl";
#     my $pdfTreeFile = plot_clade_tree( $treeFile, \@ids, $title, $clTreeFile )
# }
print "--- Computing multiplicities of condensed annotation tree leaves\n";
my @cLeaves = get_leaves( $condAnnTreeFile );
my %leafFreq; ## table of number of sequences per species
map { $leafFreq{$_}++ } @cLeaves;
my @uqLeaves = sort { $leafFreq{$b} <=> $leafFreq{$a} } keys %leafFreq;
my @uqLeaves1 = grep { $leafFreq{$_} > 1 } @uqLeaves;

if ( keys %txChg > 0 )
{
  print "\n\nWARNING: The following annotations have changed\n";
  print "oldTx\tnewTx\n";
  print_formated_tbl( \%txChg );

  print "\n\nDetailed list of <seq IDs> <orig tx> <new tx> columns was written to $chgFile\n"
}
else
{
  unlink( $chgFile );
}

my @missingAnn = diff( \@uqAnn, \@uqNewAnn );
if ( @missingAnn )
{
  print "\n\nWARNING: The following annotations are not present in the new annotation table\n";
  print_array( \@missingAnn );
}

if ( @uqLeaves1 )
{
  print "\n\nWARNING: Annontations appearing more than once on the condensed tree\n";
  print_formated_tbl( \%leafFreq, \@uqLeaves1 );
}

print "\n\nOutput written to $clFile\n\n";

####################################################################
##                               SUBS
####################################################################


sub plot_grandparent_tree
{
  my ( $treeFile, $rids, $title, $clTreeFile ) = @_;

  my ($idsFH, $idsFile) = tempfile( "tmp.XXXX", SUFFIX => 'seqID', OPEN => 1, DIR => $tmpDir );
  map { print $idsFH "$_\n" } @{$rids};
  close $idsFH;


  # make sure @ids has commas below
  # relabel the leaf labels using seqID_ann [_q] labels
  # maybe color with red _q's and _query's

  my ($fh, $cladeTreeFile) = tempfile( "tmp.XXXX", SUFFIX => 'tree', OPEN => 1, DIR => $tmpDir );
  close $fh;

  my $Rscript = qq~

    tr <- read.tree(\"$treeFile\")
    ids <- read.table(\"$idsFile\")[,1]
    p <- first.common.ancestor( tr, ids )
    p <- tree.parent.node( tr, p ) # get the grandparent of p
    ids2.idx <- node.leaves( tr, p )
    ids2 <- tr\$tip.label[ids2.idx]
    cmd <- sprintf(\"nw_clade %s \", treeFile)
    for ( id in ids2 )
  {
  cmd <- paste(cmd, id)
  }
  cmd <- paste(cmd, \" > \", \"$cladeTreeFile\")
    system(cmd)

    cmd <- sprintf("nw_rename %s %s > %s", cladeTreeFile, \"$annFile\", pFile2)
    system(cmd)

    tip.colors <- rep( 1, length(tr\$tip.label) )
    names(tip.colors) <- clade.tr\$tip.label
    tip.colors[rids] <- 4 # blue

    figH <- 8
    figW <- 12
    if ( nLeaves >= 50 )
  {
  figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
  }

  tr\$node.label <- rep("", length(tr\$node.label))

    pdf(pdfFile, width=figW, height=figH)
    op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
    plot(tr, type=\"phylogram\", use.edge.length=FALSE, no.margin=FALSE, show.node.label=T, cex=0.8, tip.color=tip.colors, main=title)
    nodelabels(cl, p)
    par(op)
    dev.off()

    ~;

  unlink $idsFile;
}

# print elements of a hash table
sub print_tbl
{
  my ($rTbl, $r) = @_;
  if ( $r )
  {
    map {print "$_\t" . $rTbl->{$_} . "\n"} @{$r};
  }
  else
  {
    map {print "$_\t" . $rTbl->{$_} . "\n"} keys %{$rTbl};
  }
}

# Parse vicut's Cltrs table
sub parse_cltr_tbl
{
  my $cltrFile = $vicutDir . "/minNodeCut.cltrs";

  if ( ! -e $cltrFile )
  {
    warn "\n\n\tERROR: vicut cltr file $cltrFile does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl; # tbl{cl}{tx} ref to seq IDs of tx taxonomy within the cl cluster.
  open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR";
  my $header = <IN>;
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    $cl = "c$cl";
    ##print "$i: $id  $cl  $tx\n";
    $tx =~ s/Unclassified/NA/;
    push @{$tbl{$cl}{$tx}}, $id;
  }
  close IN;

  return %tbl;
}

sub build_ann_seqID_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $annTreeFile; nw_rename $treeFile $annFile2 | nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


sub build_ann_cltr_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  # $annFile file format

  # readId	clstr	annot
  # Chitinivibrio	0	Fibrobacteres
  # Ignavibacterium	1	Chlorobi
  # Caldisericum	2	Caldiserica

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$3\"__\"\$2}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $annTreeFile; nw_rename $treeFile $annFile2 | nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}



sub build_ann_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  ## get leaves and test if they have the same number of elements as the table or
  ## the table has those plus more

  $cmd = "rm -f $annTreeFile; nw_rename $treeFile $annFile | nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## plot tree with clade colors
sub plot_tree
{
  my ($treeFile, $title, $pdfFile, $type) = @_;

  # type: a character string specifying the type of phylogeny to be
  #         drawn; it must be one of "phylogram" (the default),
  #         "cladogram", "fan", "unrooted", "radial" or any unambiguous
  #         abbreviation of these.

  if ( !$type )
  {
    $type = "phylogram";
  }
  my $showBoostrapVals = "T";

  if (!defined $title)
  {
    $title = "";
  }

  my $readNewickFile = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";

  my $Rscript = qq~

    source(\"$readNewickFile\")
    require(phytools)

    tr <- read.newick(file=\"$treeFile\")
    tr <- collapse.singles(tr)

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
    plot(tr, type=\"$type\", no.margin=FALSE, use.edge.length=FALSE, show.node.label=F, cex=0.8, main=\"$title\")
    par(op)
    dev.off()
    ~;

  run_R_script( $Rscript );
}

# execute an R-script
sub run_R_script
{
  my $Rscript = shift;

  my ($fh, $inFile) = tempfile("rTmpXXXX", SUFFIX => '.R', OPEN => 1, DIR => $tmpDir);
  print $fh "$Rscript";
  close $fh;

  my $outFile = $inFile . "out";
  my $cmd = "R CMD BATCH --no-save --no-restore-data $inFile $outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

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

## Create a condensed tree given a tree with certain labels at the leaves
sub condense_tree_only
{
  my ($treeFile, $condTreeFile) = @_;

  $cmd = "rm -f $condTreeFile; nw_condense $treeFile > $condTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

## Create a condensed tree given a tree with certain labels at the leaves
sub condense2_tree_only
{
  my ($treeFile, $condTreeFile) = @_;

  $cmd = "rm -f $condTreeFile; nw_condense2 $treeFile > $condTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

## Create a condensed tree given a tree, translation table and a file name of the
## condensed tree
sub condense_tree
{
  my ($treeFile, $txFile, $condTreeFile) = @_;

  my %txTbl = read_tbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  # my $treeLeavesFile = "tmp_tree.leaves";
  # $cmd = "rm -f $treeLeavesFile; $nw_labels -I $treeFile > $treeLeavesFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  # my @treeLeaves = read_array ( $treeLeavesFile );
  # unlink ( $treeLeavesFile );

  # if ( !setequal( \@txs, \@treeLeaves ) )
  # {
  #   warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
  #   print "\n\n";
  #   exit 1;
  # }

  my $sppTreeFile = "tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; nw_rename $treeFile $txFile > $sppTreeFile";
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

sub run_vicut
{
  my ($treeFile, $annFile, $vicutDir) = @_;

  my %annTbl = read_tbl( $annFile );
  my @ann = keys %annTbl;
  my @leaves = get_leaves( $treeFile );

  my @query = diff( \@leaves, \@ann );
  if ( @query )
  {
    my ($fh, $qFile) = tempfile("query.XXXX", SUFFIX => 'txt', OPEN => 1, DIR => $tmpDir);
    for ( @query )
    {
      print $fh "$_\n";
    }
    close $fh;

    my $cmd = "vicut $quietStr -t $treeFile -a $annFile -q $qFile -o $vicutDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }
  else
  {
    my $cmd = "vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
  }

  return @query;
}

sub run_vicut_no_query
{
  my ($treeFile, $annFile, $vicutDir) = @_;

  my $cmd = "vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub cleanup_tmp_files
{
  my $cmd = "rm -f $tmpDir/*";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my ($file, $skipHeader) = @_;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( $skipHeader )
  {
    my $header = <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub get_leaves
{
  my $treeFile = shift;

  my ($fh, $leavesFile) = tempfile("leaves.XXXX", SUFFIX => '', OPEN => 1, DIR => $tmpDir);
  close $fh;
  my $cmd = "nw_labels -I $treeFile > $leavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @a = read_array($leavesFile);

  return @a;
}

# common part of two arrays
sub comm
{
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

# read table with one column
sub read_array
{
  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -e $file )
  {
    warn "\n\n\nERROR: $file does not exist";
    print "\n\n";
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

# print elements of a hash table so that arguments are aligned
sub print_formated_tbl{

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
    print "WARNING: tbl value not defined for $_\n" if !exists $rTbl->{$_};
    print "$_$pad" . $rTbl->{$_} . "\n";
  }
  #print "\n";
}

# write array to a file (one column format)
sub write_array
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# write hash table to a file
sub write_tbl
{
  my ($rTbl, $outFile, $header) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  print OUT $header if $header;
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# write hash table to a file
sub write_sorted_tbl
{
  my ($rTbl, $r, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  #map {print OUT $_ . "\t" . $tbl{$_} . "\n"} @{$r};
  for ( @{$r} )
  {
    if ( exists $tbl{$_} )
    {
      print OUT $_ . "\t" . $tbl{$_} . "\n"
    }
  }
  close OUT;
}

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

exit 0;
