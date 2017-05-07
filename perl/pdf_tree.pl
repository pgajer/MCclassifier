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

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_MC_models_dir

  pdf_tree.pl --show-tree -i model.tree -o model_tree_upTx.pdf

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

GetOptions(
  "tree-file|i=s"   => \my $treeFile,
  "pdf-file|o=s"    => \my $pdfFile,
  "title|t=s"       => \my $title,
  "show-tree"       => \my $showTree,
  "igs"             => \my $igs,
  "johanna"         => \my $johanna,
  "debug"           => \my $debug,
  "dry-run"         => \my $dryRun,
  "help|h!"         => \my $help,
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

if ( !defined $title )
{
  $title = "";
}

my $readNewickFile = "/Users/pgajer/.Rlocal/read.newick.R";

if ( defined $igs )
{
  $readNewickFile = "??";
}

if ( defined $johanna )
{
  $readNewickFile = "/Users/jholm/MCclassifier/perl/read.newick.R";
}

####################################################################
##                               MAIN
####################################################################

my $tmpDir = "temp_dir";
my $cmd = "mkdir -p $tmpDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed: $?" if !$dryRun;

$pdfFile = abs_path($pdfFile);

plot_tree(abs_path($treeFile), $pdfFile, $title);

if ( $showTree && $OSNAME eq "darwin")
{
  my $cmd = "open $pdfFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


####################################################################
##                               SUBS
####################################################################

sub plot_tree
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
