#!/usr/bin/env perl

=head1 NAME

  phylo_split.pl

=head1 DESCRIPTION

  The script generates a split of all bacterial (and Archeal) sequences into
  phylogenetically-based phylo-groups.

  1. Run cut_tree on family cluster tree.

  2. Generate stats of cluster content and sizes

=head1 SYNOPSIS

  phylo_split.pl -m <max phGr size> -t <tree file> -s <cltr size file> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--tree-file, -t>
  Tree file.

=item B<--cltr-size-file, -s>
  Cluster sizes file.

=item B<--max-phGr-size, -m>
  Maximum size of a phylo-group.

=item B<--output-dir, -o>
  Output dir with the following files.

=item B<--report>
  Prints report.

=item B<--run-all>
  Ignore if ( ! -e ... ) statements.

=item B<--verbose, -v>
  Prints content of some output files. Default value: 5000.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/microcontax

  phylo_split.pl -m 4100 -s banfield_medoids_FL_dir/family_cltr_sizeR.txt -t banfield_medoids_FL_dir/family_cond_cltr.tree -o banfield_medoids_FL_dir

  phylo_split.pl -m 4100 -s banfield_medoids_FL_may30_dir/genus_cltr_sizeR.txt -t banfield_medoids_FL_may30_dir/genus_cond_cltr.tree -o banfield_medoids_FL_may30_dir

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

#my $maxCltrSize = 4100;

GetOptions(
  "tree-file|t=s"       => \my $treeFile,
  "cltr-size-file|s=s"  => \my $cltrSizeFile,
  "max-phGr-size|m=s"   => \my $maxCltrSize,
  "output-dir|o=s"      => \my $outDir,
  "report"              => \my $report,
  "run-all"             => \my $runAll,
  "quiet"               => \my $quiet,
  "igs"                 => \my $igs,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
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
elsif ( !$treeFile )
{
  print "ERROR: Missing tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$cltrSizeFile )
{
  print "ERROR: Missing cluster size file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $treeFile || ! -s $treeFile )
{
  warn "ERROR: Tree file $treeFile does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

if ( ! -e $cltrSizeFile || ! -s $cltrSizeFile )
{
  warn "ERROR: Cluter size file does not exist (or has size 0)";
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

my $nw_labels             = "nw_labels";
my $nw_order              = "nw_order";
my $nw_condense           = "nw_condense";
my $nw_rename             = "nw_rename";
my $nw_prune              = "nw_prune";
my $nw_reroot             = "nw_reroot";
my $nw_clade              = "nw_clade";
my $uc2clstr2             = "uc2clstr2.pl";
my $extract_seq_IDs       = "extract_seq_IDs.pl";
my $select_seqs           = "select_seqs.pl";
my $rmGaps                = "rmGaps";
my $FastTree              = "FastTree";
my $R                     = "R";
my $fix_fasta_headers     = "fix_fasta_headers.pl";
my $mothur                = "mothur";
my $usearch6              = "usearch6.0.203_i86osx32";
my $vicut                 = "vicut";
my $readNewickFile        = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";
my $ginsi                 = "/usr/local/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)

if ( defined $igs )
{
  $nw_labels             = "/usr/local/projects/pgajer/bin/nw_labels";
  $nw_order              = "/usr/local/projects/pgajer/bin/nw_order";
  $nw_condense           = "/usr/local/projects/pgajer/bin/nw_condense";
  $nw_rename             = "/usr/local/projects/pgajer/bin/nw_rename";
  $nw_prune              = "/usr/local/projects/pgajer/bin/nw_prune";
  $nw_reroot             = "/usr/local/projects/pgajer/bin/nw_reroot";
  $nw_clade              = "/usr/local/projects/pgajer/bin/nw_clade";
  $uc2clstr2             = "/home/pgajer/devel/MCclassifier/perl/uc2clstr2.pl";
  $extract_seq_IDs       = "/home/pgajer/devel/MCclassifier/perl/extract_seq_IDs.pl";
  $select_seqs           = "/home/pgajer/devel/MCclassifier/perl/select_seqs.pl";
  $rmGaps                = "/usr/local/projects/pgajer/bin/rmGaps";
  $FastTree              = "/home/pgajer/bin/FastTree_no_openMP";
  $R                     = "/home/pgajer/bin/R";
  $fix_fasta_headers     = "/home/pgajer/devel/MCclassifier/perl/fix_fasta_headers.pl";
  $mothur                = "/usr/local/projects/pgajer/bin/mothur";
  $usearch6              = "/local/projects/pgajer/bin/usearch6.0.203_i86linux32";
  $vicut                 = "/usr/local/projects/pgajer/bin/vicut";
  $readNewickFile        = "/local/projects/pgajer/devel/MCclassifier/perl/read.newick.R";
  $ginsi                 = "/home/pgajer/bin/mafft --maxiterate 1000 --globalpair"; # MAFFT v7.310 (2017/Mar/17)
}

####################################################################
##                               MAIN
####################################################################

print "--- Running cut_tree\n";
$cmd = "cut_tree -m $maxCltrSize -s $cltrSizeFile -t $treeFile -o $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Parsing cut_tree results\n";
my $cutFile = "$outDir/thld_" . $maxCltrSize . "_tree_cut_tbl.txt";
my ($rphGr, $rphGrSize) = parse_cut_file( $cutFile );

my %phGrTbl  = %{$rphGr};     # phGr => ref to array of vicut clusters of that phylo-group
my %phGrSize = %{$rphGrSize}; # phGr => size

print "--- Generating phGr condensed tree with size labels\n";
#print "--- Generating phGr condensed tree\n";
my $phGrTxFile = $outDir . "/thld_" . $maxCltrSize . "_phGr_size.tx";
open OUT, ">$phGrTxFile" or die "Cannot open $phGrTxFile for writing: $OS_ERROR";
for my $phGr ( keys %phGrTbl )
{
  my $size = $phGrSize{$phGr};
  my @ids = @{ $phGrTbl{$phGr} };
  for ( @ids )
  {
    print OUT "$_\t" . $phGr . "__n$size\n";
  }
}
close OUT;

my $condPhGrSizeTreeFile = $outDir . "/cond_thld_" . $maxCltrSize . "_phGr_size.tree";
condense_tree( $treeFile, $phGrTxFile, $condPhGrSizeTreeFile );

print "--- Parsing vicut family minNodeCut.cltrs file\n";
## I assume that $outDir is vicut_and_plot.pl script's output dir as well
## I also assume that the processing is taking place at the family level
my $vCltrFile = "$outDir/family_vicut_dir/minNodeCut.cltrs";
my %cltrGeTbl = parse_cltr_file( $vCltrFile ); # cltrGeTbl{$cltr} is a ref to array of genera in that cluster

print "--- Parsing microcontax pkg genus-lineage tbl\n";
my $geLiFile = "/Users/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this file has a header !
my ($rgeFaTbl, $rgeOrTbl, $rgeClTbl, $rgePhTbl) = ge_li_tbl( $geLiFile );

my %geFaTbl = %{$rgeFaTbl};
my %geOrTbl = %{$rgeOrTbl};
my %geClTbl = %{$rgeClTbl};
my %gePhTbl = %{$rgePhTbl};

print "--- Untangling phylo-group taxonomy\n";
my $ident = "\t";
for my $phGr ( sort { $phGrSize{$b} <=> $phGrSize{$a} } keys %phGrTbl )
{
  my @cltrs = @{ $phGrTbl{$phGr} };

  #print "\n\nphGr: $phGr\ncltrs: @cltrs\n";

  my @ges;
  map { push @ges, @{$cltrGeTbl{$_}} } @cltrs; # these are genera names and Banfield taxonomies

  # print "ges\n";
  # for ( sort @ges )
  # {
  #   print "\t$_\n";
  # }
  # print "\n";

  my %faTbl;
  my %faGeTbl;
  my %orTbl;
  my %orFaTbl;
  my %clTbl;
  my %clOrTbl;
  my %phTbl;
  my %phClTbl;

  my @banfield;
  for my $ge ( @ges )
  {
    if ( exists $geFaTbl{$ge} )
    {
      my $fa = $geFaTbl{$ge};
      $faTbl{$fa}++;
      $faGeTbl{$fa}{$ge}++;

      my $or = $geOrTbl{$ge};
      $orTbl{$or}++;
      $orFaTbl{$or}{$fa}++;

      my $cl = $geClTbl{$ge};
      $clTbl{$cl}++;
      $clOrTbl{$cl}{$or}++;

      my $ph = $gePhTbl{$ge};
      $phTbl{$ph}++;
      $phClTbl{$ph}{$cl}++;
    }
    else
    {
      push @banfield, $ge;
    }
  }

  my @fas = keys %faTbl;
  my @ors = keys %orTbl;
  my @cls = keys %clTbl;
  my @phs = keys %phTbl;

  if (0)
  {
    print "fas\n";
    for ( sort { $faTbl{$b} <=> $faTbl{$a} } @fas )
    {
      print "\t$_\t" . $faTbl{$_} . "\n";
    }
    print "\n";

    print "ors\n";
    for ( sort { $orTbl{$b} <=> $orTbl{$a} } @ors )
    {
      print "\t$_\t" . $orTbl{$_} . "\n";
    }
    print "\n";

    print "cls\n";
    for ( sort { $clTbl{$b} <=> $clTbl{$a} } @cls )
    {
      print "\t$_\t" . $clTbl{$_} . "\n";
    }
    print "\n";

    print "phs\n";
    for ( sort { $phTbl{$b} <=> $phTbl{$a} } @phs )
    {
      print "\t$_\t" . $phTbl{$_} . "\n";
    }
    print "\n";
  }


  print "\n-------------------------------\n";
  print "\n$phGr   (" . $phGrSize{$phGr} .")\n\n";
  for my $ph ( @phs )
  {
    print "$ph\n";
    my @clss = keys %{ $phClTbl{$ph} };
    for my $cl ( @clss )
    {
      print "$ident$cl\n";
      if ( !exists $clOrTbl{$cl} )
      {
	warn "\n\n\tERROR: $cl does not exist in clOrTbl";
	print "\n\n";
	exit 1;
      }
      my @orss = keys %{ $clOrTbl{$cl} };
      for my $or ( @orss )
      {
	print "$ident$ident$or\n";
      }
    }
  }

  print "\n";
  for my $tx ( @banfield )
  {
    print "$tx\n";
  }
  print "\n";

}

#plot_color_tree();

cleanup_tmp_files();

####################################################################
##                               SUBS
####################################################################

sub ge_li_tbl
{
  my $file = shift;

  my %geFaTbl;
  my %geOrTbl;
  my %geClTbl;
  my %gePhTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  my $header = <IN>;
  for ( <IN> )
  {
    chomp;
    my @f = split /\s+/;

    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;

    $geFaTbl{$ge} = $fa;
    $geOrTbl{$ge} = $or;
    $geClTbl{$ge} = $cl;
    $gePhTbl{$ge} = $ph;
  }
  close IN;

  return (\%geFaTbl, \%geOrTbl, \%geClTbl, \%gePhTbl);
}


# Parsing vicut clustering results on genus tree with family annotation coming
# from the medoids data
sub parse_cltr_file
{
  my $file = shift;

  # file format

  # Id            	clstr	annot
  # Acanthopleuribacter	0	Acanthopleuribacteraceae
  # Bacteria_FibrobacteresAcidobacteria_group_Acidobacteria_unclassified_Acidobacteria_Thermoanaerobaculum_aquaticum_MP_01	1	NA
  # Bacteria_unclassified_bacteria_CG_uncultured_03	2	NA
  # Chroococcidiopsis	3	Xenococcaceae
  # Symploca	4	Phormidiaceae

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  my $header = <IN>;
  for ( <IN> )
  {
    chomp;
    my ($id, $cltr, $ann) = split /\s+/;
    $cltr = "c$cltr";
    push @{ $tbl{$cltr} }, $id;
  }
  close IN;

  return %tbl;
}

# Parsing the output of vicut_and_plot.pl script
sub parse_cut_file
{
  my $cutFile = shift;

  my %phGr;
  my %phGrSize;

  open IN, "$cutFile" or die "Cannot open $cutFile for reading: $OS_ERROR";
  for ( <IN> )
  {
    chomp;
    my ($id, $cltr, $size) = split /\s+/;
    $cltr = "phGr$cltr";
    push @{ $phGr{$cltr} }, $id;
    $phGrSize{$cltr} = $size;
  }
  close IN;

  return (\%phGr, \%phGrSize);
}

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}


# extract sequence count of a fasta file
sub seq_count
{
  my $file = shift;

  my $wcline = qx/ grep -c '>' $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
}

sub build_tree
{
  my ($algnFile, $treeFile) = @_;

  $cmd = "$FastTree -nt $algnFile > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub reroot_tree
{
  my ($treeFile, $rog, $rrTreeFile ) = @_;

  my @ogs = @{$rog};

  $cmd = "$nw_reroot $treeFile @ogs | $nw_order -  > $rrTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub build_ann_seqID_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile2 | $nw_order -  > $annTreeFile";
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

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile2 | $nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}



sub build_ann_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile | $nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub plot_color_tree
{
  my ($treeFile, $clFile, $pdfFile, $title) = @_;

  my $showBoostrapVals = "F";

  if (!defined $title)
  {
    $title = "";
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
        tip.colors[i] <- "brown"
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

  run_R_script( $Rscript );
}

## plot tree with clade colors
sub plot_tree
{
  my ($treeFile, $title, $pdfFile) = @_;

  my $showBoostrapVals = "T";

  if (!defined $title)
  {
    $title = "";
  }

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
plot(tr, type=\"phylogram\", no.margin=FALSE, show.node.label=F, cex=0.8, main=\"$title\")
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
  my $cmd = "$R CMD BATCH --no-save --no-restore-data $inFile $outFile";
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

# print elements of a hash table
sub print_tbl
{
  my ($rTbl, $r) = @_;

  map {print "$_\t" . $rTbl->{$_} . "\n"} @{$r};
}

# this is a version of print_formated_tbl() where sorting w/r values are
# performed within this routine
sub print_formated_freq_tbl
{
  my $rTbl = shift;

  my @args = sort { $rTbl->{$b} <=> $rTbl->{$a} } keys %{$rTbl};
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
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  print "\n";

  return @args;
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
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  #print "\n";
}

sub mothur_align_and_add
{
  my ($candidateFile, $templateFile, $nProc) = @_;

  my $seqCountBefore = seq_count( $templateFile );

  my @tmp;
  push (@tmp,"align.seqs(candidate=$candidateFile, template=$templateFile, flip=T, processors=$nProc)");

  print_array( \@tmp, "mothur commands" ) if ( $debug || $verbose );

  my $scriptFile = create_mothur_script( \@tmp );
  my $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  my $mothurAlgnFile = $candidateFile;
  $mothurAlgnFile =~ s/fa$/align/;
  if ( ! -e $mothurAlgnFile )
  {
    warn "\n\n\tERROR: Count not find $mothurAlgnFile";
    print "\n\n";
    exit 1;
  }

  $cmd = "cat $mothurAlgnFile >> $templateFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $seqCountAfter = seq_count( $templateFile );

  if ( !$quiet )
  {
    print "\nAlignment file line count BEFORE: $seqCountBefore\n";
    print   "Alignment file line count AFTER:  $seqCountAfter\n\n";
  }

  return ($seqCountBefore, $seqCountAfter);
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

# write hash table to a file
sub write_tbl
{
  my ($rTbl, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# write hash table to a file
sub write_sorted_tbl
{
  my ($rTbl, $r, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} @{$r};
  close OUT;
}

sub run_vicut_no_query
{
  my ($treeFile, $annFile, $vicutDir) = @_;

  my $cmd = "$vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub cleanup_tmp_files
{
  my $cmd = "rm -f $tmpDir/*";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## Create a condensed tree given a tree, translation table and a file name of the
## condensed tree
sub condense_tree
{
  my ($treeFile, $txFile, $condTreeFile) = @_;

  my %txTbl = read_tbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  my $treeLeavesFile = "tmp_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; $nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @treeLeaves = read_array ( $treeLeavesFile );
  unlink ( $treeLeavesFile );

  # if ( !setequal( \@txs, \@treeLeaves ) )
  # {
  #   warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
  #   print "\n\n";
  #   exit 1;
  # }

  my $sppTreeFile = "tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; $nw_rename $treeFile $txFile > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $cmd = "rm -f $condTreeFile; $nw_condense $sppTreeFile > $condTreeFile";
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

exit 0;
