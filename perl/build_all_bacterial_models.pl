#!/usr/bin/env perl

=head1 NAME

  build_all_bacterial_models.pl

=head1 DESCRIPTION

  Generating MC models for all bacterial (and Archaeal) species using 16S rRNA
  sequences split into phylo-groups whose list is in an input file.

=head1 SYNOPSIS

  build_all_bacterial_models.pl -i <phGr list file> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--phGr-file, -i>
  Input fasta file.

=item B<--out-dir, -o>
  Output dir.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--quiet>
  Do not print progress messages.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/phylo_groups/v0.3/cx_hb_rdp_FL_5500_phGr_dir

  build_all_bacterial_models.pl -i V3V4_dirs -o V3V4_MC_models

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::Util qw( sum );
use List::MoreUtils qw(uniq);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $offsetCoef   = 0.9;
my $txSizeThld   = 10;

GetOptions(
  "phGr-file|i=s"       => \my $phGrFile,
  "out-dir|o=s"         => \my $mcDir,
  "offset-coef|c=f"     => \$offsetCoef,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
  "skip-err-thld"       => \my $skipErrThld,
  "quiet"               => \my $quiet,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$phGrFile )
{
  print "\n\nERROR: Missing phylo-group list file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$mcDir )
{
  print "\n\nERROR: Missing output dir\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $phGrFile )
{
  warn "\n\n\tERROR: $phGrFile does not exist";
  print "\n\n";
  exit 1;
}


my $quietStr = "";
if ( $quiet  )
{
  $quietStr = "--quiet";
}

my $debugStr = "";
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

my $skipErrThldStr = "";
if ($skipErrThld)
{
  $skipErrThldStr = "--skip-err-thld";
}

my $ppEmbeddingStr = "";
if ($ppEmbedding)
{
  $ppEmbeddingStr = "--pp-embedding"
}

my $cmd = "rm -rf $mcDir; mkdir -p $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;


####################################################################
##                               MAIN
####################################################################

print "--- Building fasta file\n";
my $faFile = $mcDir . "/all.fa";
build_fa_file( $phGrFile, $faFile );

print "--- Identifying duplicate sequences\n";
my @seqIDs = get_seqIDs_from_fa( $faFile );
my @uniqSeqIDs = uniq( @seqIDs );
my $NRseqIDFile = $mcDir . "/all_nr.seqIDs";
open OUT, ">$NRseqIDFile" or die "Cannot open $NRseqIDFile for appending: $OS_ERROR";
foreach my $s (@uniqSeqIDs)
  {
    print OUT "$s\n";
  }
close OUT;

print "--- Building species lineage file\n";
my $spLiFile = $mcDir . "/spp.lineage";
my %delSpp = build_spLi_file( $phGrFile, $spLiFile );

print "\nDeleting the following species appearing in multiple phylo-groups\n";
my @del = keys %delSpp;
print_array( \@del );
print "\n";

print "--- Building taxon file\n";
my $txFile = $mcDir . "/spp.tx";
build_tx_file( $phGrFile, $txFile, \%delSpp );

print "--- Removing duplicate sequences from taxonomy file\n";
my $NRtxFile = $mcDir . "/all_nr.tx";
$cmd = "select_tx.pl -i $txFile -s $NRseqIDFile -o $NRtxFile; mv $NRtxFile $txFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Removing duplicate sequences from fasta file\n";
my $NRfaFile = $mcDir . "/all_nr.fa";
$cmd = "select_seqs.pl -i $faFile -s $txFile -o $NRfaFile; mv $NRfaFile $faFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Removing duplicate sequences from lineage file\n";
my $NRLinFile = $mcDir . "/all_nr.lineage";
$cmd = "select_fullTx.pl -t $txFile -f $spLiFile -o $NRLinFile; mv $NRLinFile $spLiFile";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Building model tree and creating taxon's reference fasta files\n";
$cmd = "buildModelTree $quietStr -l $spLiFile -i $faFile -t $txFile -o $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Building MC models\n";
$cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Estimating error thresholds\n";
$cmd = "est_error_thlds -v --offset-coef $offsetCoef --tx-size-thld $txSizeThld -d $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Running classify on the ref fasta file\n";
$cmd = "classify $skipErrThldStr -d $mcDir -i $faFile -o $mcDir"; # $ppEmbeddingStr
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Comparing ref seq's taxonomy with the classification results\n";
$cmd = "cmp_tx.pl --verbose $quietStr -i $txFile -j $mcDir/MC_order7_results.txt -o $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "--- Comparing ref seq's taxonomy with the classification results\n";
$cmd = "get_lineage.pl -a $txFile -l $spLiFile -o $mcDir/sppSeqID.lineage";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

my $LineageFile = "$mcDir/sppSeqID.lineage";
my %Lineage = read2colTbl($LineageFile);
my $clTxFile = "$mcDir/MC_order7_results.txt";
my %clTx = readTbl($clTxFile);
my @clIDs   = keys %clTx;

print "--- Scoring ref seq's taxonomy with the classification results\n";
my %clScore = score_classification( \@clIDs, \%Lineage, \%clTx);
my $clScore_sum = sum values %clScore;
my $clScore_perc = 100.0 * $clScore_sum /  scalar @clIDs; 
print "\n\nUsing error thresholds, the classification score is: " . $clScore_perc . "\n";
my $clScoreFile = "$mcDir/clScore.txt";
open OUT, ">$clScoreFile" or die "Cannot open $clScoreFile for appending: $OS_ERROR";
print OUT "seqID\tcorrect_lineage\tclassification_score\n";
foreach my $key (keys %clScore)
{
  print OUT $key . "\t". $Lineage{$key} . "\t". $clScore{$key} . "\n";
}
close OUT;

print "\n\t fasta file written to $faFile\n";
print "\t taxon file written to $txFile\n";
print "\t spLineage written to $spLiFile\n";
print "\t MC master models written to $mcDir\n\n";
print "\t Classification scores written to $clScoreFile\n\n";


####################################################################
##                               SUBS
####################################################################


# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

# print elements of a hash table
sub print_tbl
{
  my ($rTbl, $r) = @_;

  if ( !$r )
  {
    my @k = keys %{$rTbl};
    $r = \@k;
  }
  map {print "$_\t" . $rTbl->{$_} . "\n"} @{$r};
}

# read table with one column
sub read_array{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
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

sub build_fa_file
{
  my ($phGrFile, $faFile) = @_;

  unlink( $faFile );

  my @dirs = read_array( $phGrFile );
  for my $dir ( @dirs )
  {
    my $phGr = $dir;
    $phGr =~ s/_dir//;
    my @p = split "_", $phGr;
    my $idx = shift @p;
    $idx =~ s/phGr//;
    my $file = $dir . "/$phGr" .  "_final.fa";
    my $cmd = "cat $file >> $faFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed: $?" if !$dryRun;
  }
}

sub build_tx_file
{
  my ($phGrFile, $txFile, $rdelSpp) = @_;

  unlink( $txFile );

  my %delSpp = %{$rdelSpp};
  my @dirs = read_array( $phGrFile );
  for my $dir ( @dirs )
  {
    my $phGr = $dir;
    $phGr =~ s/_dir//;
    my @p = split "_", $phGr;
    my $idx = shift @p;
    $idx =~ s/phGr//;
    my $file = $dir . "/$phGr" . "_final.lineage";
    open OUT, ">>$txFile" or die "Cannot open $txFile for writing: $OS_ERROR\n";
    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach ( <IN> )
    {
      chomp;
      my ($id, $li) = split /\s+/;
      my @f = split ";", $li;
      my $sp = pop @f;
      if ( !exists $delSpp{$sp} )
      {
        print OUT "$id\t$sp\n";
      }
    }
    close IN;
    close OUT;
  }
}

sub build_spLi_file
{
  my ($phGrFile, $spLiFile) = @_;

  unlink( $spLiFile );

  my %spLiTbl;

  ## If a family, order, classe or phylum has a numerical suffix '_\d+' it will
  ## made unique by attaching to it a suffix p<N>, where N is the phylo-group
  ## index. For example, if the taxon is present in phGr9_V3V4 (among other
  ## groups) it will have the suffix _p9 attached to it.

  my @dirs = read_array( $phGrFile );
  for my $dir ( @dirs )
  {
    my $phGr = $dir;
    $phGr =~ s/_dir//;
    my @p = split "_", $phGr;
    my $idx = shift @p;
    $idx =~ s/phGr//;
    my $file = $dir . "/$phGr" . "_final.lineage";
    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach ( <IN> )
    {
      chomp;
      my ($id, $li) = split /\s+/;
      my $origLi = $li;
      my @f = split ";", $li;
      my $n = @f;

      my $dom = shift @f;
      my $sp = pop @f;
      my $sge = pop @f;

      my $ge;
      if ( $n == 8 && $dom eq "Root" )
      {
        $ge = $sge;
      }
      else
      {
        $ge = pop @f;
      }

      my $fa = pop @f;
      my $or = pop @f;
      my $cl = pop @f;
      my $ph = pop @f;
      if ( @f )
      {
        $dom = pop @f;
      }

      my @t = split "_", $fa;
      if ( @t > 3 && $t[$#t] ne "etal" )
      {
        $fa = join "_", @t[0..2];
        $fa .= "_etal";
      }

      @t = split "_", $or;
      if ( @t > 3 && $t[$#t] ne "etal" )
      {
        $or = join "_", @t[0..2];
        $or .= "_etal";
      }

      @t = split "_", $cl;
      if ( @t > 3 && $t[$#t] ne "etal" )
      {
        $cl = join "_", @t[0..2];
        $cl .= "_etal";
      }

      @t = split "_", $ph;
      if ( @t > 3 && $t[$#t] ne "etal" )
      {
        $ph = join "_", @t[0..2];
        $ph .= "_etal";
      }

      $fa .= "_p$idx" if $fa =~ /\d$/;
      $or .= "_p$idx" if $or =~ /\d$/;
      $cl .= "_p$idx" if $cl =~ /\d$/;
      $ph .= "_p$idx" if $ph =~ /\d$/;

      $sge = "sg_$sge";
      $ge  = "g_$ge";
      $fa  = "f_$fa";
      $or  = "o_$or";
      $cl  = "c_$cl";
      $ph  = "p_$ph";
      $dom = "d_$dom";

      $li = "$sp\t$sge\t$ge\t$fa\t$or\t$cl\t$ph\t$dom";

      $spLiTbl{$sp}{$li}++;
      #$spLiTbl{$sp}{$li} = $phGr;
    }
    close IN;
  }

  my %delSpp; # indicator table for species appearing in multiple phylo-groups,
              # to be used in taxonomy table construction and deletion of the
              # corresponding sequences
  open OUT, ">$spLiFile" or die "Cannot open $spLiFile for writing: $OS_ERROR\n";
  for my $sp ( keys %spLiTbl )
  {
    my @lis = keys %{$spLiTbl{$sp}};

    if ( @lis > 1 )
    {
      $delSpp{$sp} = 1;
      # for my $li ( @lis )
      #   {
      #     print $spLiTbl{$sp}{$li} . ": $li\n";
      #   }
    }
    else
    {
      my $li = shift @lis;
      print OUT "$li\n";
    }
  }
  close OUT;

  return %delSpp;
}

sub score_classification{ 
  my $rclIDs = shift;
  my $rtestLineage = shift;
  my $rclTx = shift;

  my %testLineage = %$rtestLineage;
  my %clTx = %$rclTx;

  my %clScore;

  my @clIDs = @$rclIDs;
  foreach my $s (@clIDs) ## For each classified sequence
  {
    chomp $s;
    ##print "\nHere is the lineage for $s: ". $testLineage{$s} . "\n";

    ## t is equal to the classified taxa
    my $t = $clTx{$s};
    if (!$t)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in $clTxFile\n";
      }
    } 
    else
    {
      ##print "The value of $s is $clTx{$s}.\n";
    }
    ## a is equal to the reference lineage
    my $a = $testLineage{$s};
    if (!$a)
    {
      if ($debug)
      {
        print "\nWARNING: $s not found in testLineage hash table\n;"
      }
      next;
    }
    #print "\nThe value of $s is $testLineage{$s}\n";
    ## if the classified taxa is found in the reference lineage
    my $value = 0;
    if ( grep /$t/, $a) 
    {   
      print "\n $t is part of $a.\n";
      $value++;
      if ($debug)
      {
        print "\n $t is part of $a.\n";
      }
      if ( $value > 0) ##If the classifed taxonomy is found in the array
      {
        if ($clTx{$s} =~ /^d_/) 
        {
          $clScore{$s} = 0;
        }
        elsif ($clTx{$s} =~ /^p_/) 
        {
          $clScore{$s} = 0.25;
        }
        elsif ($clTx{$s} =~ /^c_/) 
        {
          $clScore{$s} = 0.325;
        }
        elsif ($clTx{$s} =~ /^o_/) 
        {
          $clScore{$s} = 0.5;
        }
        elsif ($clTx{$s} =~ /^f_/) 
        {
          $clScore{$s} = 0.625;
        }
        elsif ($clTx{$s} =~ /^g_/) 
        {
          $clScore{$s} = 0.75;
        }          
        elsif ($clTx{$s} =~ /^sg_/) 
        {
          $clScore{$s} = 0.825;
        }
        else 
        {
          $clScore{$s} = 1;
        }
      }
      elsif ($value = 0) ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
      {
        $clScore{$s} = 0;
      }   
      #print "\nThe score for $s is $clScore{$s}\n";
    }
    else ## And of course, these have to be correct with the lineage. If they are not correct, they get a zero.
    {
      $clScore{$s} = 0;
    } 
  }
return %clScore;
}

sub get_seqIDs_from_fa
{
  my $file = shift;

  my $quiet = 1;
  my $startRun = time();
  my $endRun = time();

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR\n";
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

sub read2colTbl{

  my $file = shift;

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
    if ($debug)
    {
    print "\nMaking a hash table with $id connected to $tbl{$id}\n"
    }
  }
  close IN;

  return %tbl;
}

sub readTbl{

  my $file = shift;

  if ( ! -f $file )
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

exit 0;
