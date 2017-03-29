#!/usr/bin/env perl

=head1 NAME

  plot_tx_mismatches.pl


=head1 DESCRIPTION

  This script generates a pdf file with the following three figures for each
  misclassification of a ref seq, generated by taxonomy_cleanup.pl.

  - density plot of ref/sib random seq's and all ref sequences marking in red
    those of mismatches. If the ref has enough seq's show also its density.

  - phylogenetic spp_seqIDs tree of the clade showing ref species and the
    mismatch taxon with mismatched seq's marked with red.

  - fragment of the taxon tree showing ref spp and its mismatch taxon

  The title is the summary of the mismatch.

=head1 SYNOPSIS

  plot_tx_mismatches.pl -i <group prefix> [Options]

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Groups prefix.

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

  plot_tx_mismatches.pl -i Firmicutes_group_6_V3V4

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
#use File::Basename;
use Cwd 'abs_path';

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "group-prefix|i=s"  => \my $grPrefix,
  "verbatim|v"      => \my $verbatim,
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
  exit;
}

if (!$grPrefix)
{
  warn "\n\n\tERROR: Missing group prefix";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

# [/usr/local/bin]$ sudo ln -s "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" PDFconcat

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

my $showPDFs = 1;

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "\n\n\tERROR: $grDir does not exist";
  print "\n\n";
  exit 1;
}

chdir $grDir;
print "--- Changed dir to $grDir\n";

my $mcDir        = $grPrefix . "_MC_models_dir";
my $errorDir     = $grPrefix . "_MC_models_clError_dir";
my $lineageFile  = $grPrefix . "_final.lineage";
my $csppTreeFile = $grPrefix . "_final_spp_condensed.tree"; # abs_path( $grPrefix . "_final_spp_condensed.tree" );
my $treeFile     = $grPrefix . ".tree";


# Use cmp_tx.pl output to identify mismatches

print "--- Parsing lineage table\n";
my %lineageTbl = read2colTbl($lineageFile);
my %txSeqIDs = get_tx_seq_ids(\%lineageTbl);

# Parsing comb.tx file containing ref and classified taxonomy of each sequence
print "--- Parsing comb.tx file\n";
my $cmpTxFile = "$mcDir/comb.tx";
my ($rRefTx, $rClTx) = readCmbTxTbl($cmpTxFile);

my %refTx = %{$rRefTx};
my %clTx  = %{$rClTx};

my %mmTbl;    # mismatch table
my %mmSpp;    # table of species for which mismatch was detected
my @goodSpp;  # species with perfect classification results
my %mmSeqIDs; # sequence IDs of mm species
my %spSeqIDs; # sequence IDs of each species
for my $id (keys %refTx)
{
  my $b = int($refTx{$id} eq $clTx{$id});
  if ($b==0)
  {
    $mmSpp{$refTx{$id}}++;
    push @{$mmSeqIDs{$refTx{$id}}}, $id;
  }
  $mmTbl{$refTx{$id}}{$clTx{$id}}++;
  push @{$spSeqIDs{$refTx{$id}}}, $id;
}

if ($debug)
{
  print "\nmmSpp:\n";
  printTbl(\%mmSpp);
  print "\n";
}

# Creating for each mismatch species mismatch record as in cmb_tx.pl's
# mismatch_spp.summary file that will be used as a title of each mm species plot
print "--- Generating mismatch records\n";
my %mmRec;
my %mmTxs; # mismatch taxons
for my $sp (sort keys %mmSpp)
{
  my $size = 0;
  for (keys %{$mmTbl{$sp}})
  {
    $size += $mmTbl{$sp}{$_};
  }
  my $mmStr = "$sp\t$size\n";
  for ( sort { $mmTbl{$sp}{$b} <=> $mmTbl{$sp}{$a}} keys %{$mmTbl{$sp}})
  {
    $mmStr .= "\t" . sprintf("%-30s",$_) . "\t" . sprintf("%s",$mmTbl{$sp}{$_}) . "\t" . sprintf("%.1f%%", 100 * $mmTbl{$sp}{$_} / $size ) . "\n";
    push @{$mmTxs{$sp}}, $_;
  }
  $mmRec{$sp} = $mmStr;
}


if ($debug)
{
  my @mmRecVals = values %mmRec;
  printArray(\@mmRecVals, "\nmm Records:");

  print "\nmmTxs:\n";
  printFormatedArrayValuedTbl(\%mmTxs);
  print "\n";
}

print "--- Exctracting ref and sib posterior probability data\n";

opendir(DIR, $errorDir) or die "ERROR: Cannot open $errorDir: $!";
#my $count = 0;
my %mmFile;
while (my $file = readdir(DIR))
{
  next if ($file =~ m/^\./);

  my $prefix = $file;
  $prefix =~ s/\.txt//;
  #print "file: $file\tprefix: $prefix\n"; last if $count==10;
  if ( exists $mmSpp{$prefix} )
  {
    $mmFile{$prefix} = abs_path( $file );
  }
  #$count++;
}
closedir(DIR);


print "--- Checking if mmFile has the samy keys as mmSpp\n";
my @mmSpps = keys %mmSpp;
my @mmFiles = keys %mmFile;

if ( !setequal(\@mmSpps, \@mmFiles) )
{
  warn "\n\n\tERROR: mmSpps and mmFiles do not have the same keys";
  printArray(\@mmSpps, "\nmmSpps:\n");
  printArray(\@mmFiles, "\nmmFiles:\n");
  exit 1;
}

# looping over mmSpps and generating one figure per species

my @pdfFiles;
for my $sp (@mmSpps)
{
  print "\nProcessing $sp\n";

  print "--- Writing mm seq IDs to a file\n";
  my @mmIDs = @{$mmSeqIDs{$sp}};
  my $mmIDsFile = "$mcDir/$sp" . "_mm" . ".seqIDs";
  writeArray(\@mmIDs, $mmIDsFile);

  printArray(\@mmIDs, "mm seqIDs:");
  print "\n";

  print "--- Selecting from $sp fasta file only mm seq's\n";
  my $spFaFile = "$mcDir/$sp" . ".fa";
  my $spMMfaFile = "$mcDir/$sp" . "_mm" . ".fa";
  my $cmd = "select_seqs.pl --quiet -s $mmIDsFile -i $spFaFile -o $spMMfaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


  print "--- Extracting log posterior probabilities of the mismatched sequences\n";
  print "    coming from their species' corresponding MC model\n";
  my $spPPfile = "$mcDir/$sp" . "_mm" . ".postProbs";
  $cmd = "pp_wr_selected_models -d $mcDir -i $spMMfaFile -s $sp -o $spPPfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "--- Generating pp plot\n";
  my $rFile = "$mcDir/$sp" . "_mm_pp.R";
  my $pdfFile = "$mcDir/$sp" . "_mm_pp.pdf";
  plot_pp_density($sp, $spPPfile, $rFile, $mmRec{$sp}, $pdfFile);

  if ( 0 && $showPDFs &&  $OSNAME eq "darwin")
  {
    $cmd = "open $pdfFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  push @pdfFiles, $pdfFile;

  # Building $sp and its mm taxon's clade's tree
  my @ids;
  for my $tx (@{$mmTxs{$sp}})
  {
    push @ids, @{$txSeqIDs{$tx}};
  }

  @ids = unique(\@ids);

  # my $idsFile = "$mcDir/$sp" . "_mm.ids";
  # writeArray(\@ids, $idsFile);

  my $cladeTree = "$mcDir/$sp" . "_mm_clade.tree";
  $cmd = "rm -f $cladeTree; nw_clade $treeFile @ids > $cladeTree";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


  # Extracting leaves of clade tree to compare the size of the tree to the size
  # of @ids array

  my $cladeLeaves = "$mcDir/$sp" . "_mm_clade.leaves";
  $cmd = "rm -f $cladeLeaves; nw_labels -I $cladeTree > $cladeLeaves";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @leaves = readArray($cladeLeaves);
  my @c = comm(\@leaves, \@ids);
  #if (@c != @leaves || @c != @ids)
  {
   # warn "\n\n\tERROR: leaves of $modelTreeFile and species found in the newTx table do not match";
    print "\n\n\tNumber of leaves of the clade tree: " . @leaves . "\n";
    print "\tNumber of sequences of mm taxons: " . @ids . "\n";
    print "\tNumber of sequences: " . @c . "\n\n";
  }


  print "--- Generating tree with <species name>_<seqID> labels at leaves\n";

  my %mmIDsTbl = map { $_ => 1 } @mmIDs;
  my %spIDsTbl = map { $_ => 1 } @{$txSeqIDs{$sp}};

  my %spSeqID;
  my %colorTbl;
  for my $id (@ids)
  {
    my $spID = $refTx{$id} . "_" . $id;
    $spSeqID{$id} = $spID;
    if ( exists $mmIDsTbl{$id} )
    {
      $colorTbl{$spID} = 2; # red
    }
    elsif ( exists $spIDsTbl{$id} )
    {
      $colorTbl{$spID} = 3; # green
    }
    else
    {
      $colorTbl{$spID} = 1; # black
    }
  }

  my $sppSeqIDsFile = "$mcDir/$sp" . "_mm_clade_spp.seqIDs";
  writeTbl(\%spSeqID, $sppSeqIDsFile);

  my $sppSeqIdTreeFile = "$mcDir/$sp" . "_mm_sppSeqIDs.tree";
  $cmd = "rm -f $sppSeqIdTreeFile; nw_rename $cladeTree $sppSeqIDsFile | nw_order -c n  -  > $sppSeqIdTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  # print "--- Generating clade condensed tree\n";
  # my $condSppSeqIdTreeFile = "$grPrefix" . "_mm_sppSeqIDs_cond.tree";
  # $cmd = "rm -f $condSppSeqIdTreeFile; nw_condense $sppSeqIdTreeFile > $condSppSeqIdTreeFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $colorSppSeqIDsFile = "$mcDir/$sp" . "_mm_clade_spp_seqIDs.color";
  writeTbl(\%colorTbl, $colorSppSeqIDsFile);

  print "--- Generating mm clade tree plot\n";
  my $rTreeFile = "$mcDir/$sp" . "_mm_tree.R";
  my $pdfTreeFile = "$mcDir/$sp" . "_mm_tree.pdf";
  plot_mm_clade_tree($sp, $sppSeqIdTreeFile, $colorSppSeqIDsFile, $rTreeFile, $pdfTreeFile);

  if ( 0 && $showPDFs &&  $OSNAME eq "darwin")
  {
    $cmd = "open $pdfTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  push @pdfFiles, $pdfTreeFile;
  #print "rFile: $rFile\n";
}

my $pdfFile = "$mcDir/mm.pdf";
my $cmd = "PDFconcat -o $pdfFile @pdfFiles";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

if ( $showPDFs &&  $OSNAME eq "darwin")
{
  $cmd = "open $pdfFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

####################################################################
##                               SUBS
####################################################################

# Parse cmp_tx.pl comb.tx file
sub readCmbTxTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  # File format
  # S002227729	Lactobacillus_pontis	Lactobacillus_pontis
  # S000607313	Lactobacillus_animalis_apodemi	Lactobacillus_animalis_apodemi
  # S003239402	Lactobacillus_animalis_apodemi	Lactobacillus_animalis_apodemi
  # S001112348	Lactobacillus_reuteri_vaginalis	Lactobacillus_reuteri_vaginalis
  # S000824733	Lactobacillus_brevis	Lactobacillus_brevis

  my %refTbl;
  my %clTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $refTx, $clTx) = split /\s+/,$_;
    $refTbl{$id} = $refTx;
    $clTbl{$id} = $clTx;
  }
  close IN;

  return (\%refTbl, \%clTbl);
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}


# print elements of a hash table whose values are reference to a hash table so
sub printTableValuedTbl{

  my ($rTbl, $rSub) = @_; # the second argument is a subarray of the keys of the table

  my %tbl = %{$rTbl};

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %tbl;
  }

  for (@args)
  {
    print "$_\n";
    for my $e ( keys %{$tbl{$_}} )
    {
      print "\t$e\t" . $tbl{$_}{$e} . "\n";
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

    # writeArray(\@a, "a.txt");
    # writeArray(\@b, "b.txt");
    #print "\n\tNew taxon keys and fasta IDs written to a.txt and b.txt, respectively\n\n";

    if (@a > @b)
    {
      my @d = diff(\@a, \@b);
      print "\nElements a but not b:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    if (@b > @a)
    {
      my @d = diff(\@b, \@a);
      print "\nElements in b that are not a:\n";
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

# extract unique elements from an array
sub unique
{
  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist";
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

sub plot_pp_density
{
  my ($sp, $spPPfile, $rFile, $title, $pdfFile) = @_;

  my $errorFile = "$errorDir/$sp" . ".txt";

  my $Rscript = qq~

## Reading log posterior probabilities files is more complicated as for
## higher taxons the number of these log posteriors is greater than 1000, but
## for siblings its a 1000, so we have a 'table' with different number of
## columns.

file <- \"$errorFile\"
nFields <- count.fields(file, sep = \'\\t\')
maxNumFields <- max(nFields)

tbl <- read.table(file, sep=\"\\t\", col.names = paste0(\"V\",seq_len(maxNumFields)), fill = TRUE, stringsAsFactors=FALSE)
dim(tbl)
ids <- tbl[,1]
refID <- ids[1]
idx <- ids==refID
sum(idx)
x.ref <- as.numeric(as.matrix(tbl[idx,2:nFields[1]]))
x.sib <- as.numeric(as.matrix(tbl[!idx,2:ncol(tbl)]))

if ( length(x.ref)==0 ){
   stop(\"length(x.ref) is 0\")
} else if ( length(x.sib)==0 ) {
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
p0 <- r\$root

## probability of a FP error = integral of d.sib from p0 to +inf
fpError <- 0
if ( p0 < max(d.sib\$x) ) {
  fpError <- integrate(d.sib.fun, p0, max(d.sib\$x))[[1]]
  fnError <- integrate(d.ref.fun, min(d.ref\$x), p0)[[1]]
} else {
  fpError <- 0
  fnError <- 0
}

# read mismatched seq's posterior probabilities for the $sp model
mmPPs <- read.table(\"$spPPfile\", header=F, row.names=1)

ymax <- max(c(d.ref\$y, d.sib\$y))

pdf(\"$pdfFile\", width=10, height=6)
op <- par(mar=c(4,4,6,0.5), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(d.ref, xlim=c(min(x.sib),0), xlab=\"log10[ p(x|M) ]\", ylim=c(0, ymax), las=1, main="")
title(paste(\"$title\", sprintf(\"thld=%.2f fpError=%.2f fnError=%.2f\",p0, fpError, fnError)), cex=0.6)
lines(d.sib, col=2)
if ( length(ids) > 5 ){
 legend(\"topleft\",legend=c(ids[1], \"Siblings\"), fill=c(1,2), inset=0.05, title=\"MC model (M)\", cex=1)
} else {
 legend(\"topleft\",legend=ids, fill=c(1,rep(2,length(ids)-1)), inset=0.05, title=\"MC model (M)\", cex=0.7)
}
abline(v=p0, col='gray80')
points(mmPPs[,1], rep(0, nrow(mmPPs)), pch=20, col=2)
rug(x.sib)
rug(x.ref, col=2)
dev.off()
~;

  runRscript( $Rscript, $rFile, "noErrorCheck" );
}

  # execute an R-script
sub runRscript{

  my ($Rscript, $inFile, $noErrorCheck) = @_;

  open OUT, ">$inFile",  or die "Cannot write to $inFile: $OS_ERROR";
  print OUT "$Rscript";
  close OUT;

  my $outFile = $inFile . "out";
  my $cmd = "R CMD BATCH $inFile $outFile";
  system($cmd) == 0 or die "system($cmd) failed with exit code: $OS_ERROR";

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

sub get_tx_seq_ids
{
  my $rlineageTbl = shift;

  my %lineageTbl = %{$rlineageTbl};

  my %tbl;
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

    push @{$tbl{$sp}}, $id;
    push @{$tbl{$subGe}}, $id;
    push @{$tbl{$ge}}, $id;
    push @{$tbl{$fa}}, $id;
    push @{$tbl{$or}}, $id;
    push @{$tbl{$cl}}, $id;
    push @{$tbl{$ph}}, $id;
  }

  return %tbl;
}

sub plot_mm_clade_tree
{
  my ($sp, $treeFile, $colorFile, $rTreeFile, $pdfFile) = @_;

  my $showBoostrapVals = "F";

  my $Rscript = qq~

source(\"$readNewickFile\")
require(phytools)
library(ade4)

clTbl <- read.table(\"$colorFile\", header=F)
str(clTbl)

cltr <- clTbl[,2]
names(cltr) <- clTbl[,1]

source(\"$readNewickFile\")
require(phytools)

tr1 <- read.newick(file=\"$treeFile\")
tr1 <- collapse.singles(tr1)

tip.colors <- cltr[tr1\$tip.label]

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
plot(tr1,type=\"phylogram\", no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8, tip.color=tip.colors, main=\"$sp\")
par(op)
dev.off()
~;

  runRscript( $Rscript, $rTreeFile );
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

exit 0;
