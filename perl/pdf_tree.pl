#!/usr/bin/env perl

=head1 NAME

  pdf_tree.pl

=head1 DESCRIPTION

  Create a pdf image of a tree

=head1 SYNOPSIS

  pdf_tree.pl -i <input tree file> -o <output pdf file> -t <title> [Options]

=head1 OPTIONS

=over

=item B<--tree-file, -i>
  Input tree file.

=item B<--pdf-file, -o>
  Output pdf file.

=item B<--title, -t>
  Title of the tree

=item B<--show-tree>
  Opens tree's pdf file

=item B<--tree-type>
  A character string specifying the type of phylogeny to be drawn; it must be
  one of "phylogram" (the default), "cladogram", "fan", "unrooted", "radial"
  or any unambiguous abbreviation of these.

=item B<--labs-file, -l>
  A tab delimited file with two columns <leaf ID> => <label>

=item B<--cltr-file, -c>
  A tab delimited file with two columns <leaf ID> => <cluster ID>

=item B<--condense>
  Only used with the --labs-file flag. When present a condensed tree is produced.

=item B<--condense2>
  As --condense, but adds a suffix _n<N> where N is the number of leaves that
  were collapsed to the current leaf on the condensed tree.

=item B<--show-leaf-freq>
  Prints out frequencies of leaf labels.

=item B<--show-node-labels>
  Show inner node labels.

=item B<--use-edge-length>
  Use edge lengths when plotting the tree. If not present all leafs end up at the same level.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/Banfield_contax/FL

  pdf_tree.pl --show-tree -i Rhodocyclaceae_FL_clade_mothur_algn.tree -l banfield_medoids_FL_lineage_padded_labels.txt -o Rhodocyclaceae_FL_clade_mothur_algn_lineageLabs_tree.pdf

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd qw(abs_path);
use File::Temp qw/ tempfile /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $treeType = "phylogram";

GetOptions(
  "tree-file|i=s"    => \my $treeFile,
  "pdf-file|o=s"     => \my $pdfFile,
  "labs-file|l=s"    => \my $labsFile,
  "cltr-file|c=s"    => \my $cltrFile,
  "title|t=s"        => \my $title,
  "condense"         => \my $condense,
  "condense2"        => \my $condense2,
  "show-leaf-freq"   => \my $showLeafFreq,
  "show-tree"        => \my $showTree,
  "show-node-labels" => \my $showNodeLabels,
  "tree-type"        => \$treeType,
  "use-edge-length"  => \my $useEdgeLength,
  "igs"              => \my $igs,
  "johanna"          => \my $johanna,
  "debug"            => \my $debug,
  "dry-run"          => \my $dryRun,
  "help|h!"          => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$treeFile)
{
  warn "\n\n\tERROR: Missing tree file";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$pdfFile)
{
  warn "\n\n\tERROR: Missing output file";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $treeFile )
{
  warn "\n\n\tERROR: $treeFile does not exist";
  print "\n\n";
  exit 1;
}


my $readTree       = "/Users/pgajer/organizer/programming/R/libs/read_tree.R";
my $treeLib        = "/Users/pgajer/organizer/programming/R/libs/tree.R";
my $readNewickFile = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";

if ( defined $igs )
{
  $readTree       = "/home/pgajer/devel/MCclassifier/R/read_tree.R";
  $treeLib        = "/home/pgajer/devel/MCclassifier/R/tree.R";
  $readNewickFile = "/home/pgajer/devel/MCclassifier/R/read.newick.R";
}

if ( defined $johanna )
{
  $readTree       = "/Users/jholm/MCclassifier/R/read_tree.R";
  $treeLib        = "/Users/jholm/MCclassifier/R/tree.R";
  $readNewickFile = "/Users/jholm/MCclassifier/R/read.newick.R";
}

if ( !defined $title )
{
  $title = "";
}

if ( $useEdgeLength )
{
  $useEdgeLength = "TRUE";
}
else
{
  $useEdgeLength = "FALSE";
}

if ( $showNodeLabels )
{
  $showNodeLabels = "TRUE";
}
else
{
  $showNodeLabels = "FALSE";
}



####################################################################
##                               MAIN
####################################################################

my $tmpDir = "temp_dir";
my $cmd = "mkdir -p $tmpDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed: $?" if !$dryRun;

$treeFile = abs_path( $treeFile );
$pdfFile  = abs_path( $pdfFile );

if ( $cltrFile )
{
  $cltrFile  = abs_path($cltrFile);
}

if ( $labsFile )
{
  $labsFile  = abs_path($labsFile);
}

my $annTreeFile;
if ( $labsFile && $condense )
{
  $annTreeFile = plot_condensed_tree( $treeFile, $title, $labsFile, $pdfFile);
  print "\nCondensed tree written to $annTreeFile\n\n";
}
elsif ( $labsFile && $condense2 )
{
  $annTreeFile = plot_condensed2_tree( $treeFile, $title, $labsFile, $pdfFile);
  print "\nCondensed tree written to $annTreeFile\n\n";
}
elsif ( $cltrFile && $labsFile )
{
  $annTreeFile = plot_tree_with_cltrs_and_labels( $treeFile, $title, $cltrFile, $labsFile, $pdfFile);
}
elsif ( $labsFile )
{
  $annTreeFile = plot_tree_with_labels( $treeFile, $title, $labsFile, $pdfFile);
}
else
{
  plot_tree_bw( $treeFile, $title, $pdfFile);
}

if ( $showTree && $OSNAME eq "darwin")
{
  my $cmd = "open $pdfFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

if ( $showLeafFreq )
{
  my @cLeaves;
  if ( $annTreeFile )
  {
    @cLeaves = get_leaves( $annTreeFile );
  }
  else
  {
    @cLeaves = get_leaves( $treeFile );
  }

  my %leafFreq; ## table of number of sequences per species
  map { $leafFreq{$_}++ } @cLeaves;
  my @uqLeaves = sort { $leafFreq{$b} <=> $leafFreq{$a} } keys %leafFreq;

  print "\n\nLeaf label frequencies\n";
  print_formated_tbl( \%leafFreq, \@uqLeaves );
}


####################################################################
##                               SUBS
####################################################################

sub cleanup_tmp_files
{
  my $cmd = "rm -f $tmpDir/*";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub get_leaves
{
  my $treeFile = shift;

  my $tmpDir = "temp_dir";
  my ($fh, $leavesFile) = tempfile("leaves.XXXX", SUFFIX => '', OPEN => 1, DIR => $tmpDir);
  close $fh;
  my $cmd = "nw_labels -I $treeFile > $leavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @a = read_array($leavesFile);

  my $cmd = "rm -fr $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  return @a;
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

sub plot_condensed_tree
{
  my ($treeFile, $title, $annFile, $pdfFile) = @_;

  my $annTreeFile = "tmp.tree";
  build_ann_tree( $treeFile, $annFile, $annTreeFile );

  #print "--- Generating condensed tree\n";
  my $condAnnTreeFile = $pdfFile;
  $condAnnTreeFile =~ s/pdf$/tree/;
  condense_tree_only( $annTreeFile, $condAnnTreeFile );

  plot_tree_bw( $condAnnTreeFile, $title, $pdfFile);

  unlink( $annTreeFile );

  return $condAnnTreeFile;
}

sub plot_condensed2_tree
{
  my ($treeFile, $title, $annFile, $pdfFile) = @_;

  my $annTreeFile = "tmp.tree";
  build_ann_tree( $treeFile, $annFile, $annTreeFile );

  #print "--- Generating condensed tree\n";
  my $condAnnTreeFile = $pdfFile;
  $condAnnTreeFile =~ s/pdf$/tree/;
  condense2_tree_only( $annTreeFile, $condAnnTreeFile );

  plot_tree_bw( $condAnnTreeFile, $title, $pdfFile);

  unlink( $annTreeFile );

  return $condAnnTreeFile;
}

## plot tree with colored labels
sub plot_tree_with_labels
{
  my ($treeFile, $title, $labsFile, $pdfFile) = @_;

  if (!defined $title)
  {
    $title = "";
  }

  my $annTreeFile = $pdfFile;
  $annTreeFile =~ s/pdf$/tree/;

  my $Rscript = qq~

require(ape)
source(\"$readTree\")
source(\"$treeLib\")

tr <- tree.read(\"$treeFile\")
nLeaves <- length(tr\$tip.label)

labsTbl <- read.table(\"$labsFile\", header=F, row.names=1, sep=\"\\t\")

labs <- labsTbl[,1]
names(labs) <- rownames(labsTbl)

tr\$node.label <- rep("", length(tr\$node.label))
tr\$tip.label <- as.character(labs[tr\$tip.label])

write.tree(tr, file=\"$annTreeFile\")

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3, family=\"mono\")
plot(tr, type=\"$treeType\", use.edge.length=$useEdgeLength, no.margin=FALSE, show.node.label=T, cex=0.7, tip.color=tip.colors, main=\"$title\")
par(op)
dev.off()
~;

  run_R_script( $Rscript );

  return $annTreeFile;
}


## plot tree with cluster IDs and colored labels
sub plot_tree_with_cltrs_and_labels
{
  my ($treeFile, $title, $cltrFile, $labsFile, $pdfFile) = @_;

  if (!defined $title)
  {
    $title = "";
  }

  my $annTreeFile = $pdfFile;
  $annTreeFile =~ s/pdf$/tree/;

  my $Rscript = qq~

require(ape)
source(\"$readTree\")
source(\"$treeLib\")

##tr <- tree.read(\"$treeFile\")
tr <- read.tree(\"$treeFile\")
nLeaves <- length(tr\$tip.label)

cltrTbl <- read.table(\"$cltrFile\", header=F, row.names=1, sep=\"\\t\")
cltrs <- cltrTbl[,1]
names(cltrs) <- rownames(cltrTbl)

tip.cltr <- as.character(cltrs[tr\$tip.label])

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
        tip.colors[i] <- \"brown\"
    }
}

uqLabs <- unique(tip.cltr)

cltr.node <- c()
for ( lb in uqLabs )
{
    idx <- tip.cltr==lb
    ids <- tip.cltr[idx]
    p <- first.common.ancestor( tr, ids )
    cltr.node[lb] <- p
}

labsTbl <- read.table(\"$labsFile\", header=F, row.names=1, sep=\"\\t\")
labs <- labsTbl[,1]
names(labs) <- rownames(labsTbl)

##tr\$node.label <- rep("", length(tr\$node.label))
tr\$tip.label <- as.character(labs[tr\$tip.label])

write.tree(tr, file=\"$annTreeFile\")

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3, family=\"mono\")
plot(tr, type=\"$treeType\", use.edge.length=$useEdgeLength, no.margin=FALSE, show.node.label=$showNodeLabels, cex=0.7, tip.color=tip.colors, main=\"$title\")
##nodelabels(names(cltr.node), as.vector(cltr.node))
par(op)
dev.off()
~;

  run_R_script( $Rscript );

  return $annTreeFile;
}

## plot tree with without any colors
sub plot_tree_bw
{
  my ($treeFile, $title, $pdfFile) = @_;

  if (!defined $title)
  {
    $title = "";
  }

  my $Rscript = qq~

require(ape)
source(\"$readTree\")
source(\"$treeLib\")

tr <- tree.read(\"$treeFile\")
nLeaves <- length(tr\$tip.label)

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3, family=\"mono\")
plot(tr, type=\"$treeType\", use.edge.length=$useEdgeLength, no.margin=FALSE, show.node.label=$showNodeLabels, cex=0.8, main=\"$title\")
par(op)
dev.off()
~;

  run_R_script( $Rscript );
}

sub old_plot_tree_bw
{
  my ($treeFile, $pdfFile, $title) = @_;

  if (!defined $title)
  {
    $title = "";
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
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3, family=\"mono\")
plot(tr1,type=\"phylogram\", no.margin=FALSE, show.node.label=$showNodeLabels, cex=0.8, main=\"$title\")
par(op)
dev.off()
~;

  run_R_script( $Rscript );
}


sub plot_tree
{
  my ($treeFile, $pdfFile, $title) = @_;

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
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3, family=\"mono\")
plot(tr, tip.color=tip.colors, type=\"phylogram\", no.margin=FALSE, show.node.label=$showNodeLabels, cex=0.8, main=\"$title\")
par(op)
dev.off()
~;

  run_R_script( $Rscript );
}

  # execute an R-script
sub run_R_script
{
  my ($Rscript, $noErrorCheck) = @_;

  my ($fh, $inFile) = tempfile("rTmpXXXX", SUFFIX => '.R', OPEN => 1, DIR => $tmpDir);
  print $fh "$Rscript";
  close $fh;

  my $outFile = $inFile . "out";
  my $cmd = "R CMD BATCH --no-save --no-restore-data $inFile $outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

  if (!$noErrorCheck)
  {
    open IN, "$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
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

exit 0;
