#!/usr/bin/env perl

=head1 NAME

  new_OG.pl

=head1 DESCRIPTION

  Add seq IDs to outgroup.seqIDs and update the corresponding taxon, lineage and
  alignment files. The original taxon, lineage, outgroup and alignment files are moved
  to a backup directory.

=head1 SYNOPSIS

  new_OG.pl -l <OG lineage file> -f <OG fasta file> -i <group prefix> [Options]

=head1 OPTIONS

=over

=item B<--group-prefix, -i>
  Phylo group prefix.

=item B<--og-li-file, -l>
  Lineage file of sequences that are to become new OG seq's.

=item B<--og-fa-file, -f>
  Fasta file containing sequences that are to become new OG seq's.

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

  cd /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/Proteobacteria_dir/tmp1

  new_OG.pl --debug -f ../OG_for_P10.fa -l ../OG_for_P10.lineage -i Proteobacteria_group_10_V3V4

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
  "group-prefix|i=s" => \my $grPrefix,
  "og-li-file|l=s"   => \my $ogLiFile,
  "og-fa-file|f=s"   => \my $ogFaFile,
  "keep-tmp-files"   => \my $keepTmpFiles,
  "quiet"            => \my $quiet,
  "igs"              => \my $igs,
  "verbose|v"        => \my $verbose,
  "debug"            => \my $debug,
  "dry-run"          => \my $dryRun,
  "help|h!"          => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$grPrefix )
{
  print "\n\nERROR: Missing group prefix\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$ogLiFile )
{
  print "\n\nERROR: Missing OG lineage file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$ogFaFile )
{
  print "\n\nERROR: Missing OG fasa file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $debugStr = "";
my $quietStr = "--quiet";
if ($debug)
{
  $debugStr  = "--debug";
  $quietStr  = "";
}

my $verboseStr = "";
if ($verbose)
{
  $verboseStr  = "--verbose";
}

my $nw_labels      = "nw_labels";
my $nw_order       = "nw_order";
my $nw_condense    = "nw_condense";
my $nw_rename      = "nw_rename";
my $nw_prune       = "nw_prune";
my $nw_reroot      = "nw_reroot";
my $nw_clade       = "nw_clade";
my $FastTree       = "FastTree";
my $select_seqs    = "select_seqs.pl";
my $select_tx      = "select_tx.pl";
my $mothur         = "/Users/pgajer/bin/mothur";
my $readNewickFile = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";
my $R              = "R";

if ( defined $igs )
{
  $nw_labels      = "/usr/local/projects/pgajer/bin/nw_labels";
  $nw_order       = "/usr/local/projects/pgajer/bin/nw_order";
  $nw_condense    = "/usr/local/projects/pgajer/bin/nw_condense";
  $nw_rename      = "/usr/local/projects/pgajer/bin/nw_rename";
  $nw_prune       = "/usr/local/projects/pgajer/bin/nw_prune";
  $nw_reroot      = "/usr/local/projects/pgajer/bin/nw_reroot";
  $nw_clade       = "/usr/local/projects/pgajer/bin/nw_clade";
  $FastTree       = "/home/pgajer/bin/FastTree_no_openMP";
  $select_seqs    = "/home/pgajer/devel/MCclassifier/perl/select_seqs.pl";
  $select_tx      = "/home/pgajer/devel/MCclassifier/perl/select_tx.pl";
  $mothur         = "/usr/local/projects/pgajer/bin/mothur";
  $readNewickFile = "/local/projects/pgajer/devel/MCclassifier/perl/read.newick.R";
  $R              = "/home/pgajer/bin/R";
}

####################################################################
##                               MAIN
####################################################################

my $grDir = $grPrefix . "_dir";

if ( ! -d $grDir )
{
  warn "ERROR: $grDir does not exist";
  exit 1;
}

my $tmpDir = $grDir . "/temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

$tmpDir   = abs_path( $tmpDir );
$ogLiFile = abs_path( $ogLiFile );
$ogFaFile = abs_path( $ogFaFile );

my %ogLiTbl = read_tbl( $ogLiFile );
my @ogs     = keys %ogLiTbl;

print "\nAdding " . @ogs . " sequences to the phylo-group\n\n";
print_array( \@ogs );
print "\n";

print "--- Checking consistency of the OG files\n" if !$quiet;
my $ogLiCount = line_count( $ogLiFile );
my $ogFaCount = seq_count( $ogFaFile );
if ( $ogLiCount != $ogFaCount)
{
  warn "\n\n\tERROR: $ogLiFile and $ogFaFile do not have the same number of seq's";
  print "Number of seq's in the OG lineage file: $ogLiCount\n";
  print "Number of seq's in the OG fasta file:   $ogFaCount\n\n";
  exit 1;
}


print "--- Changing dir to $grDir\n";
chdir $grDir;
##$grPrefix = "$grDir/$grPrefix";

my $liFile    = $grPrefix . ".lineage";
my $txFile    = $grPrefix . ".tx";
my $treeFile  = $grPrefix . ".tree";
my $algnFile  = $grPrefix . "_ginsi_algn.fa";
my $ogFile    = $grPrefix . "_outgroup.seqIDs";


print "--- Checking if the OG seq's are already in the lineage file\n" if !$quiet;

my %liTbl = read_tbl( $liFile );
my @seqs  = keys %liTbl;

my @c = comm( \@ogs, \@seqs );
if ( @c )
{
  print "\n\n\tWARNING: Some of new OG seq's are already in the phylo-group\n";
  print "Please remove to continue\n";
  print "Here they are\n";
  print_array( \@c );
  print "\n\n";
  exit 1;
}


print "--- Creating directory for storing the original files\n" if !$quiet;

my $bkpDir = mk_timeStamp_dir( "orig_data_dir" );

print "--- Moving the original files to the backup dir\n" if !$quiet;

my $cmd = "mv $liFile $bkpDir/$liFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "mv $txFile $bkpDir/$txFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "cp $algnFile $bkpDir/$algnFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "mv $ogFile $bkpDir/$ogFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

$cmd = "mv $treeFile $bkpDir/$treeFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;



print "--- Adding OG seq's to the lineage file\n" if !$quiet;

my $liCountBefore = line_count( "$bkpDir/$liFile" );
$cmd = "cat $ogLiFile $bkpDir/$liFile > $liFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
my $liCountAfter = line_count( $liFile );

if ( !$quiet )
{
  print "\nLineage file line count BEFORE: $liCountBefore\n";
  print   "Lineage file line count AFTER:  $liCountAfter\n\n";
}

print "--- Adding OG seq's to the taxonomy file\n" if !$quiet;
my $txCountBefore = line_count( "$bkpDir/$txFile" );

my %ogTx = map{ $_ => "OG" } @ogs;

my ($fh, $ogTxFile) = tempfile("tmp.XXXX", SUFFIX => 'tx', OPEN => 0, DIR => $tmpDir);
write_tbl( \%ogTx, $ogTxFile );

$cmd = "cat $ogTxFile $bkpDir/$txFile > $txFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
my $txCountAfter = line_count( $txFile );

if ( !$quiet )
{
  print "\nTaxonomy file line count BEFORE: $txCountBefore\n";
  print   "Taxonomy file line count AFTER:  $txCountAfter\n\n";
}

print "--- Writing OG seq's to the OG file\n" if !$quiet;
write_array( \@ogs, $ogFile );

print "--- Aligning OG seq's to the phGr alignment file\n" if !$quiet;
my $nProc = 8;
my ($seqCountBefore, $seqCountAfter) = mothur_align_and_add( $ogFaFile, $algnFile, $nProc );

if ( $seqCountBefore != $seqCountAfter )
{
  print "--- Rebuilding tree using modified alignment\n" if !$quiet;
  my ($tmpFH, $tmpTreeFile) = tempfile("tmp.XXXX", SUFFIX => 'tree', OPEN => 0, DIR => $tmpDir);
  build_tree( $algnFile, $tmpTreeFile );

  print "--- Rerooting the tree using new outgroup sequences\n" if !$quiet;
  reroot_tree( $tmpTreeFile, \@ogs, $treeFile );

  print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
  my $ssTreeFile  = $grPrefix . "_sppSeqIDs.tree";
  build_spp_seqID_tree( $treeFile, $txFile, $ssTreeFile );

  print "--- Testing if OG seq's form a monophylectic clade at the top or bottom of the tree\n";
  if ( test_OG( $treeFile, \%ogLiTbl ) != 0 )
  {
    print "--- Generating pdf of the condensed tree\n";
    my $treeAbsPath = abs_path( $ssTreeFile );
    my $pdfTreeFile = $grPrefix . "_sppSeqIDs_tree.pdf";
    plot_tree( $treeAbsPath, $grPrefix, $pdfTreeFile );

    $cmd = "open $pdfTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    warn "\n\n\tERROR: There is an issue with the tree outgroup seq's";
    print "\n\tPlease check out the following trees\n";
    print "\t$treeFile\n";
    print "\t$ssTreeFile\n";
    print "\n\n";

    exit 1;
  }

  print "\n\n\tSuccessfully added OG seq's\n";
  print "\n\tPlease check out the following trees\n";
  print "\t$treeFile\n";
  print "\t$ssTreeFile\n";
  print "\n\n";

}
else
{
  warn "\n\n\tERROR: Abandoning rebuilding of the tree as the alignment file was not modified";
  print "\n\n";
}


####################################################################
##                               SUBS
####################################################################

# extract line count of a file
sub line_count
{
  my $file = shift;

  my $wcline = qx/ wc -l $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
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

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;

  #local $" = ', ';
  local $" = "\n";
  print "$header:\n" if $header;
  print "@{$a}\n";
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify
{
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
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
sub read2colTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
    exit 1;
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

  print_arrayTbl(\%grTx, "grTx") if $debug;

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
    print_arrayTbl(\%grTx, "grTbl");
    print "\n";
    exit 1;
  }

  print_array(\@t, "t=grTx{cmin}") if $debug;

  # for (@t)
  # {
  #   $pTbl{$parent{$_}}++ if exists $parent{$_};
  # }
  # ##$p = $parent{$t[0]};# its enough to take the parent of the first taxonomic rank as all taxonomic ranks at the lowest level have to have the same parent

  # print_tbl(\%pTbl, "pTbl") if $debug;

  # for my $p (keys %pTbl)
  # {
  #   push @ch, keys %{$children{$p}};
  # }

  @ch = keys %{$children{$gp}};

  print_array(\@ch, "ch") if $debug;

  @t = ($p);
  @d = diff(\@ch, \@t);

  print_array(\@d, "d") if $debug;

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
      print_tbl(\%pTbl, "pTbl");
      print "\n";
    }
    else
    {
      my $p = $parent{$rTaxons->[0]};
      print "p: $p\n";
    }

    print "ch: @ch\n";
    print "d: @d\n";
    exit 1;
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
sub print_tbl
{
  my ($rTbl, $header) = @_;
  print "$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
  print "\n";
}

# write array to a file (one column format)
sub write_array
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT "$_\n"} @{$a};
  close OUT
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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
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

# write hash table to a file
sub write_tbl
{
  my ($rTbl, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
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

# test if OG seq's form one or more clusters in the tree
sub test_OG
{
  my ($treeFile, $rogInd) = @_;

  my %ogInd = %{$rogInd};

  my $debug_test_OG = 0;

  my $ret = 0;

  print "\t--- Extracting leaves from $treeFile\n" if $debug_test_OG;
  my $treeLeavesFile = "$grPrefix" . "_sppSeqIDs.leaves";
  my $cmd = "rm -f $treeLeavesFile; $nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\t--- Reading leaves\n" if $debug_test_OG;
  my @leaves = read_array($treeLeavesFile);

  print "\t--- Checking the number of clusters formed by OG seqs\n" if $debug_test_OG;
  my @ogIdx;
  for my $i (0..$#leaves)
  {
    if ( exists $ogInd{$leaves[$i]} )
    {
      push @ogIdx, $i;
    }
  }

  print_array(\@ogIdx, "\nPositions of OG seq's") if ($debug_test_OG);

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
	print_array(\@start, "start");
	print_array(\@end, "end");
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

  print_array( \@og, "\nOutgroup elements" ) if $debug;

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

    my $imax = argmax( \@rangeSize );
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
	$cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$grPrefix" . "_clade.leaves";
	$cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Reading the leaves\n" if $debug_test_OG;
	my @ogCladeLeaves = read_array($ogCladeTreeLeavesFile);

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

    print_array(\@og, "og");
    print "\n";

    my $maxOGbigSize = 100;
    if ( @ogBig < $maxOGbigSize )
    {
      print_array(\@ogBig, "Leaves from first to last OG seq");
    }

    print "\t--- Extracting the clade of OG sequences\n";
    my $ogCladeTreeFile = "$grPrefix" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug_test_OG;
    my @ogCladeLeaves = read_array($ogCladeTreeLeavesFile);

    my $maxCladeSize = 100;
    if ( @ogCladeLeaves < $maxCladeSize )
    {
      print_array( \@ogCladeLeaves, "OG Clade Leaves" );
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
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$grPrefix" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug_test_OG;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug_test_OG;
    my @ogCladeLeaves = read_array( $ogCladeTreeLeavesFile );

    if ( @ogCladeLeaves != @og )
    {
      warn "\n\n\tERROR: The outgroup sequences do not form a monophyletic clade!";

      my $maxCladeSize = 100;
      if ( @ogCladeLeaves < $maxCladeSize )
      {
	print_array(\@ogCladeLeaves, "OG Clade Leaves");
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

  print "\n\tOG seq's form a monophylectic clade at the top or bottom of the tree\n\n" if $debug_test_OG;

  return $ret;
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

sub mk_timeStamp_dir
{
  my $baseName = shift;

  my @now = localtime();
  my $timeStamp = sprintf("%04d%02d%02d_%02d%02d%02d",
			  $now[5]+1900, $now[4]+1, $now[3],
			  $now[2],      $now[1],   $now[0]);

  my $tmpDir = $baseName . "_$timeStamp";
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  return $tmpDir;
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

sub build_spp_seqID_tree
{
  my ($treeFile, $annFile, $ssTreeFile ) = @_;

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $ssTreeFile; $nw_rename $treeFile $annFile2 | $nw_order -  > $ssTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
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

exit 0;
