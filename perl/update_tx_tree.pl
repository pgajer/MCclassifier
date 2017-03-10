#!/usr/bin/env perl

=head1 NAME

  update_tx_tree.pl

=head1 DESCRIPTION

  Given a vicut output directory, vicutDir and a tree used with vicut, use
  vicutDir/minNodeCut.cltrs, to update taxonomy - NA sequences have taxonomy
  assigned based on the majority vote - using only non-NA sequences.

  the only case NAs stay as such is when their cluster contains only NAs


=head1 SYNOPSIS

  update_tx_tree.pl -a <tx file> -t <tree file> -d <vicut directory> [Options]

=head1 OPTIONS

=over

=item B<--tx-file, -i>
  Tab or space delimited taxonomy file.

=item B<--tree-file, -t>
  A tree fila that was used with vicut commant.

=item B<--vicut-dir, -d>
  vicut output directory, that will be used also as an output directory.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  cd ~/devel/MCclassifier/data/RDP
  update_tx_tree.pl --debug -a test_dir2/Lactobacillaceae_nr_good_wOG.tx -t test_dir2/Lactobacillaceae_nr_SATe_aln.good_w_ecoli_57_1472_trimming_rm_L_sp_good95_tf_wOG_rr.tree -d test_vicut_dir

  cd /Users/pgajer/devel/MCclassifier/data/SILVA_123/Lactobacillaceae_ss100_dir
  update_tx_tree.pl -a Lactobacillaceae.good.tx -t Lactobacillaceae_full_sate_rr.tree -d Lactobacillaceae_vicut_full_sate_dir

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Cwd 'abs_path';

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "tx-file|a=s"      => \my $txFile,
  "tree-file|t=s"    => \my $treeFile,
  "vicut-dir|d=s"    => \my $vicutDir,
  "dry-run"          => \my $dryRun,
  "debug"            => \my $debug,
  "help|h!"          => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

my $newTxFile = "$vicutDir/minNodeCut.cltrs"; # was minNodeCut_NAge1_TXge1_querySeqs.taxonomy";

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$txFile)
{
  print "\n\nERROR: Missing taxonomy file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$treeFile)
{
  print "\n\nERROR: Missing tree file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$vicutDir)
{
  print "\n\nERROR: Missing vicut directory\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if ( ! -f $txFile )
{
  print "\n\nERROR: $txFile does not exist\n\n\n";
  exit;
}

if ( ! -f $treeFile )
{
  print "\n\nERROR: $treeFile does not exist\n\n\n";
  exit;
}

if ( ! -f $newTxFile )
{
  print "\n\nERROR: $newTxFile does not exist\n\n\n";
  exit;
}

####################################################################
##                               MAIN
####################################################################



print "\t- Updating min-node-cut taxonomy using majority vote\n";

my $cltrFile = "$vicutDir/minNodeCut.cltrs";
if ( ! -f $cltrFile )
{
  print "\n\nERROR: $cltrFile does not exist\n\n\n";
  exit;
}

my ($rid2clTbl, $rclTbl, $rclFreqTbl, $rtxTbl)  = readCltrTbl($cltrFile);

my %clFreqTbl = %{$rclFreqTbl}; # cltrID => table of taxon frequencies within the given cluster
my %clTbl     = %{$rclTbl};     # cltrID => ref to array of IDs in the given cluster
my %txTbl     = %{$rtxTbl};     # seqID => taxonomy
my %id2clTbl  = %{$rid2clTbl};  # seqID => cltrID

for my $cl ( keys %clFreqTbl )
{
  my @txs = sort { $clFreqTbl{$cl}{$b} <=> $clFreqTbl{$cl}{$a} || $a cmp $b } keys %{$clFreqTbl{$cl}};
  # NOTE the sorting by given the number of sequences of the given species in
  # different clusters (starting from the largest).  If there are two clusters
  # with the same number of seq's of the given species they are going to be
  # sorted alphabetically.
  # This is quite arbitrary choice of representative species for the given cluster.
  # An alternative could be a name that is a concatenation of the two species of the same number of sequences.
  if ($debug)
  {
    print "Cluster $cl taxons\n";
    map {print "$_\t" . $clFreqTbl{$cl}{$_} . "\n"} @txs;
    print "\n";
  }

  if ( @txs > 1 && $txs[0] ne "NA") # && $clFreqTbl{$cl}{$txs[0]} > $clFreqTbl{$cl}{$txs[1]} )
    # NOTE: if there is more than on species in a vicut cluster and there are two
    # species with the same number of seq's and all other species have less
    # sequences, then the species which is alphabetically first among the ones
    # with the largest number of sequences is going to be used to propagate its
    # species name to all other sequences.
  {
    # changing taxonomy of NA sequences to the majority taxonomy
    for my $id ( @{$clTbl{$cl}} )
    {
      $txTbl{$id} = $txs[0];
    }
  }
  elsif ( @txs > 1 && $txs[0] eq "NA" )
  {
    # changing taxonomy of NA sequences to the second abundant (after NA) taxonomy
    for my $id ( @{$clTbl{$cl}} )
    {
      $txTbl{$id} = $txs[1];
    }
  }
}
## parsing tax table from before vicut taxonomy modifications

my %tx = read2colTbl($txFile);

## Generating updated taxonomy file
## changing NAs to taxonomy in the original taxonomy table
for my $id ( keys %txTbl )
{
  if ( $txTbl{$id} eq "NA" )
  {
    if ( exists $tx{$id} )
    {
      $txTbl{$id} = $tx{$id};
    }
    else
    {
      print "\nWARNING: $id not present in $txFile table\n\n";
    }
  }
}

my %spFreqTbl;
for my $id ( keys %txTbl )
{
  my $sp = $txTbl{$id};
  my $cl = $id2clTbl{$id};
  $spFreqTbl{$sp}{$cl}++;
}

my %spFreq;
for my $id ( keys %txTbl )
{
  my $sp = $txTbl{$id};
  $spFreq{$sp}++;
}

printTbl(\%spFreq, "spFreq") if $debug;

my %spIDs;
for my $id ( keys %txTbl )
{
  my $sp = $txTbl{$id};
  push @{$spIDs{$sp}}, $id;
}

print "--- Changing taxonomy of species found in more than one cluster\n";
print "    if there is a dominating cluster (size > size of others)\n";
print "    then set all small cluster's taxonomies to <genus>_sp\n";
print "    and keep the taxonomy of the largest cluster\n";
for my $sp ( keys %spFreqTbl )
{
  my @cls = sort { $spFreqTbl{$sp}{$b} <=> $spFreqTbl{$sp}{$a} } keys %{$spFreqTbl{$sp}};
  if ( @cls > 1 && $spFreqTbl{$sp}{$cls[0]} > $spFreqTbl{$sp}{$cls[1]} )
  {
    ##my ($g, $s) = split "_", $sp;

    if ($debug)
    {
      print "\n\nProcessing $sp\tnCltrs: " . @cls . "\n";
      print "Cluster sizes: ";
      map { print $spFreqTbl{$sp}{$_} . ", "} @cls;
    }

    my @f = split "_", $sp;
    my $g = shift @f;
    #my $s = shift @f;
    print "Genus: $g\n" if $debug;
    my $cmax = shift @cls;
    my $spSp = $g . "_sp";
    print "sp species: $spSp\n" if $debug;
    for my $cl (@cls)
    {
      my @spSeqIDs = comm($clTbl{$cl}, $spIDs{$sp});
      print "Processing seq's of cluster $cl of size " . @spSeqIDs . "\n" if $debug;
      for my $id ( @spSeqIDs )
      {
	$txTbl{$id} = $spSp;
      }
    }
  } # end of if ( @cls > 1
}

if ($debug)
{
  my %spFreq2;
  for my $id ( keys %txTbl )
  {
    my $sp = $txTbl{$id};
    $spFreq2{$sp}++;
  }

  printTbl(\%spFreq2, "spFreq2");
}

my $updatedTxFile = "$vicutDir/updated.tx";
writeTbl(\%txTbl, $updatedTxFile);

#print "\tUpdated taxonomy  written to $updatedTxFile\n\n";


####################################################################
##                               SUBS
####################################################################

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read2colTbl{

  my $file = shift;

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

# read cltrs and new species assignment
# elements of the first column to the second column
# remove the header line
sub readCltrTbl
{
  my $file = shift;

  my %txTbl;     # seqID => taxonomy
  my %id2clTbl;  # seqID => cltrID
  my %clTbl;     # cltrID => ref to array of IDs in the given cluster
  my %clFreqTbl; # cltrID => table of taxon frequencies within the given cluster
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  my $header = <IN>;
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    push @{$clTbl{$cl}}, $id;
    $clFreqTbl{$cl}{$tx}++;
    $txTbl{$id} = $tx;
    $id2clTbl{$id} = $cl;
  }
  close IN;

  return (\%id2clTbl, \%clTbl, \%clFreqTbl, \%txTbl);
}

# write hash table to a file
sub writeTbl{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} keys %tbl;
  close OUT;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
  print "\n\n";
}


# print elements of a hash table
sub printTbl
{
  my ($rTbl, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
  print "\n\n";
}

# count L crispatus seq's
sub countLc
{
  my $rTbl = shift;
  my $i = 0;
  for (keys %$rTbl)
  {
    if( $rTbl->{$_} =~ /crispatus/ )
    {
      $i++;
    }
  }
  return $i;
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

exit;
