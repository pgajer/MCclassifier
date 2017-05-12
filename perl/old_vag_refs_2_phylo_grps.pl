#!/usr/bin/env perl

=head1 NAME

  old_vag_refs_2_phylo_grps.pl

=head1 DESCRIPTION

  Mapping old vaginal V3V4 reference sequences to the appropriate phylo-groups

=head1 SYNOPSIS

  old_vag_refs_2_phylo_grps.pl -o <out dir> [Options]

=head1 OPTIONS

=over

=item B<--out-dir, -o>
  Output directory.

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

  cd /Users/pgajer/devel/MCextras/data/vaginal_old_vs_PECAN_May11_2017

  old_vag_refs_2_phylo_grps.pl -o vag_expanded_V3V4_phylo_groups_dir

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "out-dir|o=s"     => \my $outDir,
  "verbose|v"       => \my $verbose,
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

if ( !$outDir )
{
  print "\n\n\tERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $outDir )
{
  my $cmd = "mkdir -p $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my $debugStr = "";
my $quietStr = "--quiet";
if ($debug)
{
  $debugStr = "--debug";
  $quietStr = "";
}

####################################################################
##                               MAIN
####################################################################

my $baseDir        = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";


my @algnFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_ginsi_algn.fa",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_ginsi_algn.fa",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_ginsi_algn.fa",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4_ginsi_algn.fa",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4_ginsi_algn.fa",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4_ginsi_algn.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4_ginsi_algn.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4_ginsi_algn.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_ginsi_algn.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4_ginsi_algn.fa",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_ginsi_algn.fa",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4_ginsi_algn.fa",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_ginsi_algn.fa");

my @algnFiles = map{ $baseDir . $_ } @algnFiles0;
## print "algnFiles: @algnFiles\n";

## all that follows is for the old V3V4 vaginal ref seq's

## in ~/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013

## $ faSize vaginal_319_806_v2.fa
## 1660

## $ wc -l vaginal_319_806_v2.tx
##     1660 vaginal_319_806_v2.tx

## $ time classify -v --skip-err-thld -d /Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/all_bacteria_V3V4_MC_models_dir -i vaginal_319_806_v2.fa -o vaginal_319_806_v2_pecan_dir

print "--- Parsing PECAN tx table of old ref seq's\n";
my $oRefPeTxFile = "/Users/pgajer/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013/vaginal_319_806_v2_pecan_dir/MC_order7_results.txt";
my %oRefPeTx = read_tbl( $oRefPeTxFile );

print "--- Parsing the origninal tx table of old ref seq's\n";
my $oRefTxFile = "/Users/pgajer/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013/vaginal_319_806_v2.tx";
my %oRefTx = read_tbl( $oRefTxFile );

my %oTxIDs; # sp => ref to array of seq IDs of the old ref seq's of that species
for my $id ( keys %oRefTx )
{
  push @{ $oTxIDs{$oRefTx{$id}} }, $id;
}


print "--- Parsing species => phGr table\n";
my $sp2phGrFile = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mmDir_May5/master_sp_to_phGr.txt";
my %sp2phGr = parse_sp_2_phGr_tbl( $sp2phGrFile );

## use classifiers results to map  <old ref's> => <phylo-group>

print "--- Generating phGr => old ref seqIDs table\n";
my %phGr2oRefs; # phGr => ref to array of oRefs with PECAN classification from that phGr
my %phGr2oPeTx; # phGr => freq table of PECAN classifications
my %phGr2oTx;   # phGr => freq table the original taxonomy
my %tx2pe;      # original taxonomy => freq table of PECAN classifications
my %tx2phGr;    # original taxonomy => freq tbl of phGr's
for my $id ( keys %oRefPeTx )
{
  my $sp = $oRefPeTx{$id};
  if ( !exists $sp2phGr{$sp} )
  {
    warn "\n\n\tERROR: $sp does not exist in sp2phGr";
    print "\n\n";
    exit 1;
  }
  my $phGr = $sp2phGr{$sp};
  push @{ $phGr2oRefs{$phGr} }, $id;
  $phGr2oPeTx{$phGr}{$sp}++;

  my $origSp = $oRefTx{$id};
  $phGr2oTx{$phGr}{$origSp}++;
  $tx2pe{$origSp}{$sp}++;
  $tx2phGr{$origSp}{$phGr}++;
}

if ( $verbose )
{
  print "\nOld vs PECAN ref vaginal taxonomy grouped by phylo-groups\n";
  my $nTabs = 2;
  my $nPerfectMatches = 0; # number of original species whose are seq's are classified by PECAN to the same species
  for my $phGr ( sort keys %phGr2oTx )
  {
    print "$phGr\n";
    my %origSpFreq = %{ $phGr2oTx{$phGr} };
    my @origSpp = sort { $phGr2oTx{$phGr}{$b} <=> $phGr2oTx{$phGr}{$a} } keys %origSpFreq;
    for my $origSp ( @origSpp )
    {
      print "\t$origSp\n";
      my %peSpFreq = %{ $tx2pe{$origSp} };
      my @spp = sort { $peSpFreq{$b} <=> $peSpFreq{$a} } keys %peSpFreq;
      $nPerfectMatches++ if ( @spp == 1 && $origSp eq $spp[0] );
      print_formated_tbl(\%peSpFreq, \@spp, $nTabs);
    }
  }
  print "\n";

  print "\n\n\tNumber of perfect matches: $nPerfectMatches\n\n";

  print "Species found in more than one phylo-group\n";
  for my $sp ( keys %tx2phGr )
  {
    my %phGrFreq = %{ $tx2phGr{$sp} };
    if ( keys %phGrFreq > 1 )
    {
      print "$sp\n";
      for my $phGr ( keys %phGrFreq )
      {
	print "\t$phGr\t" . $phGrFreq{$phGr} . "\n";
      }
    }
  }
  print "\n\n";
}


print "--- Removing old ref species from phylogroup if its present in two phylogroups and its count in this one is 1\n";
for my $origSp ( keys %tx2phGr )
{
  my %phGrFreq = %{ $tx2phGr{$origSp} };
  if ( keys %phGrFreq > 1 )
  {
    for my $phGr ( keys %phGrFreq )
    {
      if ( $phGrFreq{$phGr} == 1 )
      {
	delete $phGr2oTx{$phGr}{$origSp};
	delete $tx2phGr{$origSp}->{$phGr};

	if ( keys %{ $phGr2oTx{$phGr} } == 0 )
	{
	  delete $phGr2oTx{$phGr};
	}
      }
    }
  }
}


## dereplicate old ref's with ref of phylo-groups, to see which sequences are
## already present there and if they are what taxonomy was assigned to them

# print "--- Testing if any old ref seq's are a part of the given phyl-group\n";
# for my $phGr ( keys %phGr2oTx )
# {
#     my %origSpFreq = %{ $phGr2oTx{$phGr} };
#     my @origSpp = sort { $phGr2oTx{$phGr}{$b} <=> $phGr2oTx{$phGr}{$a} } keys %origSpFreq;
#     my @oSeqIDs;
#     for my $origSp ( @origSpp )
#     {
#     }
# }

## use species linage table

# $ head /Users/pgajer/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013/vaginal_319_806_v2.fullTx

#  Bifidobacterium_longum	g_Bifidobacterium	f_Bifidobacteriaceae	o_Bifidobacteriales	c_Actinobacteria	p_Actinobacteria	d_Bacteria
#  Lactobacillus_helveticus	g_Lactobacillus	f_Lactobacillaceae	o_Lactobacillales	c_Bacilli	p_Firmicutes	d_Bacteria
#  Pseudomonas_aeruginosa	g_Pseudomonas	f_Pseudomonadaceae	o_Pseudomonadales	c_Gammaproteobacteria	p_Proteobacteria	d_Bacteria
#  Acinetobacter_baumannii	g_Acinetobacter	f_Moraxellaceae	o_Pseudomonadales	c_Gammaproteobacteria	p_Proteobacteria	d_Bacteria

## to create old ref lineage table of the form

# S001125253	Root;Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Kandleria;Kandleria_sp
# S002739411	Root;Bacteria;Firmicutes;Negativicutes;Selenomonadales;Veillonellaceae;Veillonella;Veillonella_parvula
# S002945661	Root;Bacteria;Firmicutes;Thermolithobacteria;Thermolithobacterales;Thermolithobacteraceae;Thermolithobacter;Thermolithobacter_sp

print "--- Parsing old ref's species lineage table\n";
my $oRefSppLiFile = "/Users/pgajer/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013/vaginal_319_806_v2.fullTx";
my %oRefSppLi = parse_spp_li_tbl( $oRefSppLiFile );

if ( $verbose )
{
  print "\nOld ref species lineage table\n";
  print_tbl( \%oRefSppLi );
  print "\n\n";
}


## generate alignment of old ref's to final algn of the phGr

my $msrVagFaFile = "/Users/pgajer/projects/16S_rRNA_pipeline/vaginal_species_oct18_2013/vaginal_319_806_v2.fa";

print "--- Generating extended alignment, tree, lineage and taxon files\n";

##
## main loop
##

my $phGr = "Bacteroidetes_group_2_V3V4";
#for my $phGr ( keys %phGr2oTx )
{
  print "Processing $phGr\n";

  # files to generate

  # my $lineageFile     = $grPrefix . ".lineage";
  # my $algnFile	= $grPrefix . "_ginsi_algn.fa";
  # my $trimmedAlgnFile = $grPrefix . "_ginsi_algn.fa"; # "_algn_trimmed.fa";
  # my $outgroupFile    = $grPrefix . "_outgroup.seqIDs";
  # my $treeFile	= $grPrefix . ".tree";
  # my $txFile          = $grPrefix . ".tx";

  ##
  print "\tCreating phylo-group dir\n";
  my $phGrDir = $outDir . "/" . $phGr . "_dir";
  ## print "\n\nphGrDir: $phGrDir\n"; exit;

  #next if -e $phGrDir;
  #next if $phGr eq "Actinobacteria_group_0_V3V4";

  my $cmd = "mkdir -p $phGrDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my $logFile = $phGrDir . "/log_file.txt";
  open LOUT, ">$logFile" or die "Cannot open $logFile for writing: $OS_ERROR\n";

  ##
  print "\tGetting seq IDs of old ref's that were assigned to the given phGr\n";
  my %origSpFreq = %{ $phGr2oTx{$phGr} };
  my @origSpp = sort { $phGr2oTx{$phGr}{$b} <=> $phGr2oTx{$phGr}{$a} } keys %origSpFreq;
  my @oSeqIDs;
  for my $origSp ( @origSpp )
  {
    #print "origSp: $origSp\toTxIDs{$origSp}: " . @{$oTxIDs{$origSp}} . "\n";
    push @oSeqIDs, @{ $oTxIDs{$origSp} };
  }

  if ( @oSeqIDs == 0 )
  {
    warn "\n\n\tERROR: array oSeqIDs is empty";
    print "\n\n";
    exit 1;
  }

  print "\nNumber of old vag ref seq's: " . @oSeqIDs . "\n" if $debug;
  print LOUT "\nNumber of old vag ref seq's: " . @oSeqIDs . "\n" if $debug;

  ##
  print "\tIdentifying algn file of the given phylo-group\n";
  my @f = grep { $_ =~ /$phGr/ } @algnFiles;
  my $phGrAlgnFile = $f[0];

  print "phGrAlgnFile: $phGrAlgnFile\n" if $debug;
  ## Final alignment has OG seq's

  if ( ! -e $phGrAlgnFile )
  {
    warn "\n\n\tERROR: $phGrAlgnFile does not exist";
    print "\n\n";
    exit 1;
  }

  my $baseName = $phGrAlgnFile;
  $baseName =~ s/_ginsi_algn.fa//;
  print "\nbaseName: $baseName\n" if $debug;


  ##
  print "\tIdentifying old ref seq's present in the phGr\n";

  my $phGrAlgnSeqIDsFile = $phGrDir . "/$phGr" . "_algn.seqIDs";
  $cmd = "extract_seq_IDs.pl -i $phGrAlgnFile -o $phGrAlgnSeqIDsFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my @phGrSeqIDs = read_array( $phGrAlgnSeqIDsFile );

  my @c = comm( \@oSeqIDs, \@phGrSeqIDs );
  if ( @c )
  {
    print "\nDetected " . @c . " seq's common with old ref in $phGr\n";

    print LOUT "Old vag ref sequences detected in the current phylo-group\n";
    for ( @c )
    {
      print LOUT "$_\t" . $oRefTx{$_} . "\n";
      print "$_\t" . $oRefTx{$_} . "\n";
    }
    print LOUT "\n";
    print "\n";

    print "Removing them from old refs\n\n";
    @oSeqIDs = diff( \@oSeqIDs, \@c );
  }

  ##
  print "\tGenerating fasta file of these seq's\n";
  my $idsFile = $phGrDir . "/$phGr" . "_vag_refs.seqIDs";
  write_array( \@oSeqIDs, $idsFile);

  my $faFile = $phGrDir . "/$phGr" . "_vag_refs.fa";
  $cmd = "select_seqs.pl $quietStr -s $idsFile -i $msrVagFaFile -o $faFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


  ##
  print "\tAligning old vag ref's to the phGr's alignment\n";
  #my $oAlgnFile = $phGrDir . "/$phGr" . "_vag_refs.algn";

  my @tmp;
  push (@tmp,"align.seqs(candidate=$faFile, template=$phGrAlgnFile, processors=8, flip=T)");

  printArray(\@tmp, "mothur commands") if ($debug || $verbose);

  my $scriptFile = create_mothur_script(\@tmp);
  $cmd = "mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  ##
  print "\tMerging candidate and template\n";
  my $mothurAlgnFile = $faFile;
  $mothurAlgnFile =~ s/fa$/align/;
  print "mothurAlgnFile: $mothurAlgnFile\n" if $debug;

  my $extAlgnFile = $phGrDir . "/$phGr" . "_ginsi_algn.fa";
  $cmd = "rm -f $extAlgnFile; cat $mothurAlgnFile $phGrAlgnFile > $extAlgnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ##
  print "\tCopying the original outgroup file\n";
  my $ogFile = $baseName . "_outgroup.seqIDs";
  print "ogFile: $ogFile\n" if $debug;

  my $extOGfile = $phGrDir . "/$phGr" . "_outgroup.seqIDs";
  $cmd = "cp $ogFile $extOGfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my @ogSeqs = read_array( $extOGfile );

  ##
  print "\tGenerating old vag refs lineage file\n";
  my $oLiFile = $phGrDir . "/$phGr" . "_vag_refs.lineage";
  open OUT, ">$oLiFile" or die "Cannot open $oLiFile for writing: $OS_ERROR";
  for my $id ( @oSeqIDs )
  {
    my $sp = $oRefTx{$id};
    print OUT "$id\t" . $oRefSppLi{$sp} . "\n";
  }
  close OUT;

  ##
  print "\tGenerating extended lineage file\n";
  my $liFile = $baseName . ".lineage";
  print "liFile: $liFile\n" if $debug;

  my $extLiFile = $phGrDir . "/$phGr" . ".lineage";
  $cmd = "rm -f $extLiFile; cat $oLiFile $liFile > $extLiFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;


  ##
  print "\tGenerating old vag refs taxon file\n";
  my $oTxFile = $phGrDir . "/$phGr" . "_vag_refs.tx";
  open OUT, ">$oTxFile" or die "Cannot open $oTxFile for writing: $OS_ERROR";
  for my $id ( @oSeqIDs )
  {
    my $sp = $oRefTx{$id};
    print OUT "$id\t$sp\n";
  }
  close OUT;

  ##
  print "\tGenerating extended taxon file\n";
  my $txFile = $baseName . ".tx";
  print "txFile: $txFile\n" if $debug;

  my $extTxFile = $phGrDir . "/$phGr" . ".tx";
  $cmd = "rm -f $extTxFile; cat $oTxFile $txFile > $extTxFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  ##
  print "\tBuilding tree\n";
  my $unrootedTreeFile = $phGrDir . "/$phGr" . "_unrooted.tree";
  $cmd = "rm -f $unrootedTreeFile; FastTree -nt $extAlgnFile > $unrootedTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "\tExtracting the clade of OG sequences\n";
  my $ogCladeTreeFile = $phGrDir . "/$phGr" . "_og_clade.tree";
  $cmd = "rm -f $ogCladeTreeFile; nw_clade $unrootedTreeFile @ogSeqs > $ogCladeTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  print "\tExtracting leaves of the OG clade\n";
  my $ogCladeTreeLeavesFile = $phGrDir . "/$phGr" . "_og_clade.leaves";
  $cmd = "rm -f $ogCladeTreeLeavesFile; nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my @ogClade = read_array( $ogCladeTreeLeavesFile );

  ## setting these seq's to be outgroup seq's
  @ogSeqs = @ogClade;
  write_array( \@ogSeqs, $extOGfile );

  my $extTreeFile = $phGrDir . "/$phGr" . ".tree";
  $cmd = "nw_reroot $unrootedTreeFile @ogSeqs > $extTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # print "--- Generating tree with <species name>_<seqID> labels at leaves\n";
  # my $sppSeqIDsFile  = $phGrDir . "/$phGr" . ".sppSeqIDs";
  # $cmd = "rm -f $sppSeqIDsFile; awk '{print \$1\"\\t\"\$2\"__\"\$1}' $extTxFile > $sppSeqIDsFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  # $cmd = "rm -f $sppSeqIDsTreeFile; nw_rename $grTreeFile $sppSeqIDsFile | nw_order -  > $sppSeqIDsTreeFile";
  # print "\tcmd=$cmd\n" if $dryRun || $debug;
  # system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  close LOUT;
  #exit;
}

## run tx cleanup

## build master models

## classify MM

## run validation scripts to validate the tx assignments

## assign taxonomy to MM seq's

####################################################################
##                               SUBS
####################################################################

## turning

#  Bifidobacterium_longum	g_Bifidobacterium	f_Bifidobacteriaceae	o_Bifidobacteriales	c_Actinobacteria	p_Actinobacteria	d_Bacteria
#  Lactobacillus_helveticus	g_Lactobacillus	f_Lactobacillaceae	o_Lactobacillales	c_Bacilli	p_Firmicutes	d_Bacteria
#  Pseudomonas_aeruginosa	g_Pseudomonas	f_Pseudomonadaceae	o_Pseudomonadales	c_Gammaproteobacteria	p_Proteobacteria	d_Bacteria

## into

## Veillonella_parvula 	Root;Bacteria;Firmicutes;Negativicutes;Selenomonadales;Veillonellaceae;Veillonella;Veillonella_parvula

## table
sub parse_spp_li_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    next if /^$/;
    my ($sp, $ge, $fa, $or, $cl, $ph, $do) = split /\s+/, $_;

    if ( !$do )
    {
      warn "\n\n\tERROR: do undefined; probably wrong number of level in the lineage line $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    if ( $ge =~ /g_/ )
    {
      $ge =~ s/g_//;
    }
    else
    {
      warn "\n\n\tERROR: genus does not have g_ prefix in $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    if ( $fa =~ /f_/ )
    {
      $fa =~ s/f_//;
    }
    else
    {
      warn "\n\n\tERROR: family does not have f_ prefix in $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    if ( $or =~ /o_/ )
    {
      $or =~ s/o_//;
    }
    else
    {
      warn "\n\n\tERROR: order does not have o_ prefix in $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    if ( $cl =~ /c_/ )
    {
      $cl =~ s/c_//;
    }
    else
    {
      warn "\n\n\tERROR: class does not have c_ prefix in $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    if ( $ph =~ /p_/ )
    {
      $ph =~ s/p_//;
    }
    else
    {
      warn "\n\n\tERROR: phylum does not have p_ prefix in $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    if ( $do =~ /d_/ )
    {
      $do =~ s/d_//;
    }
    else
    {
      warn "\n\n\tERROR: domain does not have d_ prefix in $_";
      print "\n\n";
      print "lineage file: $file\n\n";
      exit 1;
    }

    $tbl{$sp} = "Root;$do;$ph;$cl;$or;$fa;$ge;$sp";
  }
  close IN;

  return %tbl;

}

# print elements of a hash table so that arguments are aligned
sub print_formated_tbl
{
  my ($rTbl, $rSub, $nTabs) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ( $rSub )
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  my $tabs = "";
  if ( $nTabs )
  {
    for (1..$nTabs)
    {
      $tabs .= "\t";
    }
  }

  for ( @args )
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "$tabs$_$pad" . $rTbl->{$_} . "\n";
  }
  print "\n";
}


# parse sp_to_phGr.txt table

# Ilumatobacter_sp_4	Actinobacteria_group_0_V3V4
# Aciditerrimonas_sp_9	Actinobacteria_group_0_V3V4
# Aciditerrimonas_sp_1	Actinobacteria_group_0_V3V4
# Aciditerrimonas_sp_13	Actinobacteria_group_0_V3V4

sub parse_sp_2_phGr_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    next if /^$/;
    my ($sp, $phGr ) = split /\s+/, $_;
    $tbl{$sp} = $phGr;
  }
  close IN;

  return %tbl;
}

# print elements of a hash table
sub print_tbl
{
  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# read two column table
sub read_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in read_tbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    next if /^$/;
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}


sub create_mothur_script{

    my (@arr) = @{$_[0]};
    my $file = "mothur_script.txt"; ##tmpnam();
    open OUT, ">$file" or die "Cannot open file $file to write: $!\n";
    foreach my $c (@arr){
        print OUT $c . "\n";
    }
    print OUT "quit()\n";

    return $file;
}

# write array to a file (one column format)
sub write_array
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

# read table with one column
sub read_array
{
  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in read_array(): $file does not exist";
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

exit 0;
