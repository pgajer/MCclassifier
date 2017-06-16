#!/usr/bin/env perl

=head1 NAME

  genus_size.pl  -o <output dir> [Options]

=head1 DESCRIPTION

  Generate a two column table of our db genus sizes

  <genus> <size>

=head1 SYNOPSIS

  genus_size.pl [Options]

=head1 OPTIONS

=over

=item B<--output-dir, -o>
  Output dir

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

  cd ~/projects/PECAN/data/Banfield_contax

  genus_size.pl -o genus_size_dir

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

my $maxCltrSize = 3000;
my $nProc       = 8;
my $spSizeThld  = 100;
GetOptions(
  "output-dir|o=s"      => \my $outDir,
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

my $baseDir               = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
my $spGeFile              = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/species_genus_tbl_may19_2017.txt";
my $geLiFile              = "/Users/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this file has a header !
# my $ogFaFile              = "/Users/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4.fa";
# my $ogTxFile              = "/Users/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4_v2.tx";

my $txFile                = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/master_V3V4_no_outliers.tx";

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

  $spGeFile              = "/usr/local/projects/pgajer/devel/MCextras/data/species_genus_tbl_may19_2017.txt";
  $geLiFile              = "/usr/local/projects/pgajer/devel/MCextras/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this one has a header !
  #$ogFaFile              = "/usr/local/projects/pgajer/devel/MCextras/data/phylo_split/Archaea_Caldococcus_noboribetus_S000414080.fa";

  $baseDir               = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
  $spGeFile              = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_groups/v0.2/species_genus_tbl_may19_2017.txt";
  $geLiFile              = "/usr/local/projects/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this file has a header !
  # $ogFaFile              = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4.fa";
  # $ogTxFile              = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4_v2.tx";
}

####################################################################
##                               MAIN
####################################################################

# creating output directory
my $cmd = "mkdir -p $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Generating our db genus size tables\n";
my ($rsizeTbl, $rsizeTblR, $rspSeqIDs, $rspSeqIDsR) = get_size_tbls();

my %sizeTbl  = %{$rsizeTbl};  # taxon => size
my %sizeTblR = %{$rsizeTblR}; # taxon => size with the uppler limit on the number of representatives of at the species level

my %spSeqIDs  = %{$rspSeqIDs};
my %spSeqIDsR = %{$rspSeqIDsR};

my $sizeFile = $outDir . "/genus_size.txt";
write_tbl( \%sizeTbl, $sizeFile );

my $sizeFileR = $outDir . "/genus_sizeR.txt";
write_tbl( \%sizeTblR, $sizeFileR );

my $outFile  = $outDir . "/all.tx";
my $outFileR = $outDir . "/size_restricted.tx";
open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
open OUTR, ">$outFileR" or die "Cannot open $outFileR for writing: $OS_ERROR\n";
for my $sp ( keys %spSeqIDs )
{
  my @ids = @{ $spSeqIDs{$sp} };
  for ( @ids )
  {
    print OUT "$_\t$sp\n";
  }

  @ids = @{ $spSeqIDsR{$sp} };
  for ( @ids )
  {
    print OUTR "$_\t$sp\n";
  }
}
close OUT;
close OUTR;

print "\nGenus size tables written to\n";
print "$sizeFile\n";
print "$sizeFileR\n";

print "\nTaxon and restricted taxon files\n";
print "$outFile\n";
print "$outFileR\n";
print "\n\n";

####################################################################
##                               SUBS
####################################################################

## constructing tables reporting the number of sequences we have in our db at
## each taxonomic rank
sub get_size_tbls
{
  print "--- Parsing seqID => species table\n";
  my %spTbl = read_tbl( $txFile );

  print "--- Parsing species-genus tbl\n";
  my %spGeTbl = read_tbl( $spGeFile );

  # print "--- Building genus tables\n";
  # my %geTbl;
  # foreach my $id ( keys %spTbl )
  # {
  #   my $sp = $spTbl{$id};

  #   if ( ! exists $spGeTbl{$sp} )
  #   {
  #     my ($g, $suffix) = split "_", $sp;
  #     print "WARNING: $sp not found in spGeTbl - using $g\n";
  #     $spGeTbl{$sp} = $g;
  #   }
  #   my $ge = $spGeTbl{$sp};
  #   $geTbl{$id} = $id;
  # }

  print "--- Building species-genus tables\n";
  my %spSeqIDs;
  my %geSpFreq;
  foreach my $id ( keys %spTbl )
  {
    my $sp = $spTbl{$id};

    if ( ! exists $spGeTbl{$sp} )
    {
      my ($g, $suffix) = split "_", $sp;
      print "WARNING: $sp not found in spGeTbl - using $g\n";
      $spGeTbl{$sp} = $g;
    }
    my $ge = $spGeTbl{$sp};
    $geSpFreq{$ge}{$sp}++;
    push @{ $spSeqIDs{$sp} }, $id;
  }

  print "--- Computing genus sizes (with and without capping of species sizes)\n";

  my %sizeTbl;  # genus => size
  my %sizeTblR; # genus => size
  my %spSeqIDsR; # seqID => species (with random selection of seqIDs for species with more than $spSizeThld seq's)

  for my $ge ( keys %geSpFreq )
  {
    print "\rProcessing $ge                                      ";
    my %spFreq = %{ $geSpFreq{$ge} };
    my @spp = keys %spFreq;

    my $geSize = 0;
    map { $geSize += $spFreq{$_} } @spp;
    $sizeTbl{$ge} = $geSize;

    for ( @spp )
    {
      if ( $spFreq{$_} > $spSizeThld )
      {
	$spFreq{$_} = $spSizeThld;
	my @ids = @{ $spSeqIDs{$_} };
	fisher_yates_shuffle( \@ids );
	@{ $spSeqIDsR{$_} } = @ids[0..($spSizeThld-1)];
      }
      else
      {
	@{ $spSeqIDsR{$_} } = @{ $spSeqIDs{$_} };
      }
    }

    $geSize = 0;
    map { $geSize += $spFreq{$_} } @spp;
    $sizeTblR{$ge} = $geSize;
  }
  print "\r                                                       \r";
  return (\%sizeTbl, \%sizeTblR, \%spSeqIDs, \%spSeqIDsR);
}

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle
{
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

sub get_li_files
{
  my @liFiles0 = ("Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4.lineage",
		  "Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4.lineage",
		  "Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4.lineage",
		  "Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4.lineage",
		  "Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4.lineage",
		  "Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4.lineage",
		  "Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4.lineage",
		  "Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4.lineage",
		  "Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4.lineage",
		  "Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4.lineage",
		  "Chloroflexi_V3V4_dir/Chloroflexi_V3V4.lineage",
		  "Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4.lineage",
		  "Fusobacteria_V3V4_dir/Fusobacteria_V3V4.lineage",
		  "Nitrospirae_V3V4_dir/Nitrospirae_V3V4.lineage",
		  "Planctomycetes_V3V4_dir/Planctomycetes_V3V4.lineage",
		  "Spirochaetes_V3V4_dir/Spirochaetes_V3V4.lineage",
		  "Tenericutes_V3V4_dir/Tenericutes_V3V4.lineage",
		  "Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4.lineage",
		  "phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4.lineage",
		  "Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4.lineage",
		  "Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4.lineage",
		  "Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4.lineage",
		  "Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4.lineage",
		  "Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4.lineage",
		  "Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4.lineage",
		  "Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4.lineage",
		  "Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4.lineage",
		  "Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4.lineage",
		  "Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4.lineage",
		  "Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4.lineage",
		  "Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4.lineage",
		  "Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4.lineage",
		  "Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4.lineage",
		  "Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4.lineage",
		  "Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4.lineage",
		  "Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4.lineage",
		  "Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4.lineage",
		  "Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4.lineage",
		  "Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4.lineage",
		  "Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4.lineage",
		  "Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4.lineage",
		  "Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4.lineage",
		  "Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4.lineage");

  my @liFiles = map{ $baseDir . $_ } @liFiles0;

  ## testing for existence
  for my $file ( @liFiles )
  {
    if ( ! -e $file )
    {
      warn "\n\n\tERROR: $file does not exist";
      print "\n\n";
      exit 1;
    }
  }

  return @liFiles;
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

# test if OG seq's form one or more clusters in the tree
sub test_OG
{
  my ($treeFile, $rogInd) = @_;

  my %ogInd = %{$rogInd};

  my $debug_test_OG = 0;

  my $ret = 0;

  print "\t--- Extracting leaves from $treeFile\n" if $debug_test_OG;
  my $treeLeavesFile = "$outDir" . "/sppSeqIDs.leaves";
  my $cmd = "rm -f $treeLeavesFile; $nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\t--- Reading leaves\n" if $debug;
  my @leaves = read_array($treeLeavesFile);

  print "\t--- Checking the number of clusters formed by OG seqs\n" if $debug;
  my @ogIdx;
  for my $i (0..$#leaves)
  {
    if ( exists $ogInd{$leaves[$i]} )
    {
      push @ogIdx, $i;
    }
  }

  print_array(\@ogIdx, "\nPositions of OG seq's") if ($debug);

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

      if (0 && $debug)
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
    warn "\n\n\tERROR: start and end arrays have different lengths!";
    print "length(start): " . @start . "\n";
    print "length(end): " . @end . "\n\n";
    $ret = 1;
  }

  my @rangeSize;
  for my $i (0..$#start)
  {
    push @rangeSize, ($end[$i] - $start[$i]+1);
  }

  if ($debug)
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
	my $ogCladeTreeFile = "$outDir/genus" . "_clade.tree";
	$cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$outDir/genus" . "_clade.leaves";
	$cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Reading the leaves\n" if $debug;
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
    my $ogCladeTreeFile = "$outDir/genus" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$outDir/genus" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug;
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
    print "\t--- Extracting the clade of OG sequences\n" if $debug;
    my $ogCladeTreeFile = "$outDir/genus" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$outDir/genus" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug;
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

  print "\n\tOG seq's form a monophylectic clade at the top or bottom of the tree\n\n" if $debug;

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

# print elements of a hash table whose values are reference to a hash table so
sub print_tbl_valued_tbl
{
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

## Create a condensed tree given a tree with certain labels at the leaves
sub condense_tree_only
{
  my ($treeFile, $condTreeFile) = @_;

  $cmd = "rm -f $condTreeFile; $nw_condense $treeFile > $condTreeFile";
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

sub get_tbls_from_lineage
{
  my $file = shift;

  my %geTbl;   # seqID => genus
  my %faTbl; # genus => family
  my %orTbl; # genus => order
  my %clTbl; # genus => class
  my %phTbl; # genus => phylum

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    chomp;
    my @f = split /\s+/;

    my $id = shift @f;
    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;
    #my $do = shift @f;

    $geTbl{$id} = $ge;
    $faTbl{$id} = $fa;
    $orTbl{$id} = $or;
    $clTbl{$id} = $cl;
    $phTbl{$id} = $ph;
  }
  close IN;

  return (\%geTbl, \%phTbl, \%clTbl, \%orTbl, \%faTbl);
}

exit 0;
