#!/usr/bin/env perl

=head1 NAME

  classify_mm_spp.pl

=head1 DESCRIPTION

  Get the master report generated by get_master_mm_report.pl and classify call
  species (species to which M&M sequences were classified)

  Goal 1. Identify call species with one cluster consisting of that species and
  the NA group (query seq's).

=head1 SYNOPSIS

  classify_mm_spp.pl

=head1 OPTIONS

=over

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

  classify_mm_spp.pl

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
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


####################################################################
##                               MAIN
####################################################################

my $mmDir = "/Users/pgajer/projects/M_and_M/new_16S_classification_data";

print "--- Parsing file with PECAN generated taxonomy on M&M's sequences\n";
my $qTxFile = "/Users/pgajer/devel/MCextras/data/mm_pecan_tx2_no_controls.txt";
my %spSize = readQtxTbl($qTxFile); # sp => number of seq's classified to this species

my $nSeqs = 0;
for (keys %spSize)
{
  $nSeqs += $spSize{$_};
}

print "--- Parsing table of species detected in the M&M dataset\n";
my $sppTblFile = $mmDir . "/mm_uq_spp_report_no_controls.txt";
my %spToPhGr = parseSpTbl($sppTblFile); # spToPhGr: <sp> => <sp's phylo-group>

my @allSpp = keys %spToPhGr;
my $nAllSpp = @allSpp;


## Extracting species records

my %rTbl; # <species name> => <record with cluster info text>
my %recComplete; # $recComplete{$sp} = 1 if the record contains 'Cluster' and 'range(pp) strings in it. Otherwise recComplete is not defined on $sp

my $sp;
my $startRec = 0;
my $foundCluster = 0;
my $inPhGrHeader = 0;

my $masterReportFile = $mmDir . "/master_mm_validate_report.txt";

my $file = $masterReportFile;
open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
for (<IN>)
{
  if (/^---/)
  {
    $startRec = 1;
    $foundCluster = 0;
  }
  elsif ( /^sp: (\w+)/ && $startRec )
  {
    ## checking if the old sp has a complete record
    if ( defined $sp && exists $rTbl{$sp} && !exists $recComplete{$sp} )
    {
      # print "\n\n\tWARNING: $sp has incomplete record\n";
      # print $rTbl{$sp};
      # print "\tDeleting it from rTbl table\n";
      delete $rTbl{$sp};
    }

    chomp;
    $sp = $1;
    $startRec = 0;
  }
  else
  {
    if (/Cluster/)
    {
      $foundCluster = 1;
    }
    elsif (/range/ && $foundCluster)
    {
      $recComplete{$sp} = 1;
    }
    elsif (/^===/ && !$inPhGrHeader)
    {
      $inPhGrHeader = 1;
    }
    elsif (/^===/ && $inPhGrHeader)
    {
      $inPhGrHeader = 0;
    }

    if ( defined $sp && !/^===/ && !$inPhGrHeader )
    {
      $rTbl{$sp} .= $_;
    }
  }
}
close IN;

my %spCltrs; # <sp> => <ref to array of cluster hash tables>
my %spMatch; # <sp> => 1 if one of the clusters contains sp
my $clStr;
my $clID;
my $sizeStr;
my %clElSize; # cluster hash table: <element> => <size of the cluster element>
my $cltrBlankLineCounter = 0;
my $nSpp = 0;
undef $sp;

open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
for ( <IN> )
{
  chomp;
  print "line: $_\n" if $debug;

  if ( /^sp: (\w+)/ )
  {
    $nSpp++;

    if ( defined $sp ) # this means that this is not the first record; get previous record's cluster data
    {
      %{$spCltrs{$sp}} = %clElSize;
    }

    if ( $debug )
    {
      my $nCltrs = keys %clElSize;
      print "in sp: $sp\tnCltrs: $nCltrs\n";
      for my $cl (keys %clElSize)
      {
	print "cl: $cl\n";
	my %elts = %{ $clElSize{$cl} };
	for my $el (keys %elts)
	{
	  print "\t$el\t" . $elts{$el} . "\n";
	}
	print "\n";
      }
      print "\n";
    }

    $sp = $1;
    undef %clElSize;
  }
  elsif ( /^n/ )
  {
    next;
  }
  elsif ( /range/ )
  {
    $foundCluster = 0;
  }
  elsif ( /Cluster/ )
  {
    $foundCluster = 1;
    $spMatch{$sp} = 0;
    ($clStr, $clID, $sizeStr) = split /\s+/;
  }
  elsif ( /^$/ && $foundCluster )
  {
    $cltrBlankLineCounter++;
    if ( $cltrBlankLineCounter== 2 )
    {
      $cltrBlankLineCounter = 0;
      $foundCluster = 0;
    }
  }
  elsif ( !/^$/ && !/^---/ && $foundCluster )
  {
    my ($el, $size) = split /\s+/;

    if ( !defined $size )
    {
      print "\n\n\tWARNING: size undefined: line: $_\n\n";
      exit;
    }

    $el =~ s/://;
    $clElSize{$clID}{$el} = $size+0;
    if ( $el eq $sp )
    {
      $spMatch{$sp} = 1;
    }
  }
}
close IN;

my $nType1spNAspp      = 0;
my $nType1spNAseqs     = 0;

my $nType1spNAplusSpp  = 0;
my $nType1spNAplusSeqs = 0;

my $nType1NAseqs        = 0;
my $nType1spNAminusSeqs = 0;
my $nType2spNAplusSeqs  = 0;
my $nType2spNAminusSeqs = 0;

my @type1spNAspp;
my @type1spNAplusSpp;
my @type1NAspp;
my @type1spNAminusSpp;
my @type2spNAplusSpp;
my @type2spNAminusSpp;

my $type1spNAfileSpp = $mmDir . "/type_1spNA_spp.txt";
my $type1spNAplusFileSpp = $mmDir . "/type_1spNA_plus_spp.txt";
open OUT1, ">$type1spNAfileSpp" or die "Cannot open $type1spNAfileSpp for writing: $OS_ERROR";
open OUT2, ">$type1spNAplusFileSpp" or die "Cannot open $type1spNAplusFileSpp for writing: $OS_ERROR";
print OUT1 "sp\tspSize\tnaSize\n"; # header
for my $sp ( keys %spCltrs )
{
  my %cltrs = %{$spCltrs{$sp}};
  my $nCltrs = keys %cltrs;
  my @clIDs = keys %cltrs;

  if ( !exists $spMatch{$sp} )
  {
    print "\n\n\tERROR: $sp not in spMatch\n\n";
    exit;
  }

  if ( $spMatch{$sp}==1 && $nCltrs==1 )
  {
    %clElSize = %{ $cltrs{$clIDs[0]} };
    my @elts = keys %clElSize;
    my @spNA = ($sp, "NA");
    if ( setequal(\@elts, \@spNA) )
    {
      my $spSize = $clElSize{$sp};
      my $naSize = $clElSize{"NA"};
      print OUT1 "$sp\t$spSize\t$naSize\n";
      push @type1spNAspp, $sp;
      $nType1spNAspp++;
      $nType1spNAseqs += $spSize{$sp};
    }
    else
    {
      print OUT2 "$sp";
      for (@elts)
      {
	print OUT2 "\t$_";
      }
      print OUT2 "\n";
      push @type1spNAplusSpp, $sp;
      $nType1spNAplusSpp++;
      $nType1spNAplusSeqs += $spSize{$sp};
    }
  }
  elsif ( $spMatch{$sp}==0 && $nCltrs==1 )
  {
    %clElSize = %{ $cltrs{$clIDs[0]} };
    my @elts = keys %clElSize;
    if ( @elts==1 )
    {
      # type 1NA
      push @type1NAspp, $sp;
      $nType1NAseqs += $spSize{$sp};
    }
    else
    {
      # type 1spNA-
      push @type1spNAminusSpp, $sp;
      $nType1spNAminusSeqs += $spSize{$sp};
    }
  }
  elsif ( $spMatch{$sp}==1 && $nCltrs>1)
  {
    # type 2spNA+
    push @type2spNAplusSpp, $sp;
    $nType2spNAplusSeqs += $spSize{$sp};
  }
  else
  {
    # $spMatch{$sp}==0 && $nCltrs>1
    # type 2spNA-
    push @type2spNAminusSpp, $sp;
    $nType2spNAminusSeqs += $spSize{$sp};
  }
}
close OUT1;
close OUT2;


my @allTypes;
push @allTypes, @type1spNAspp;
push @allTypes, @type1spNAplusSpp;
push @allTypes, @type1NAspp;
push @allTypes, @type1spNAminusSpp;
push @allTypes, @type2spNAplusSpp;
push @allTypes, @type2spNAminusSpp;

my @d = diff(\@allSpp, \@allTypes);
if (@d)
{
  print "\n\n\td: @d\n";
}

my $type1spNAfile = $mmDir . "/type_1spNA_rec.txt";
open OUT, ">$type1spNAfile" or die "Cannot open $type1spNAfile for writing: $OS_ERROR";
for my $sp ( @type1spNAspp )
{
  if ( !exists $rTbl{$sp} )
  {
    print "\n\n\tERROR: $sp not in rTbl\n\n";
    next;
  }
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
close OUT;

my $type1spNAplusFile = $mmDir . "/type_1spNA_plus_rec.txt";
open OUT, ">$type1spNAplusFile" or die "Cannot open $type1spNAplusFile for writing: $OS_ERROR";
for my $sp ( @type1spNAplusSpp )
{
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
close OUT;


my $type1NAfile = $mmDir . "/type_1NA_rec.txt";
open OUT, ">$type1NAfile" or die "Cannot open $type1NAfile for writing: $OS_ERROR";
for my $sp ( @type1NAspp )
{
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
close OUT;


my $type1spNAminusFile = $mmDir . "/type_1spNA_minus_rec.txt";
open OUT, ">$type1spNAminusFile" or die "Cannot open $type1spNAminusFile for writing: $OS_ERROR";
for my $sp ( @type1spNAminusSpp )
{
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
close OUT;

my $type2spNAplusFile = $mmDir . "/type_2spNA_plus_rec.txt";
open OUT, ">$type2spNAplusFile" or die "Cannot open $type2spNAplusFile for writing: $OS_ERROR";
for my $sp ( @type2spNAplusSpp )
{
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
close OUT;

my $type2spNAminusFile = $mmDir . "/type_2spNA_minus_rec.txt";
open OUT, ">$type2spNAminusFile" or die "Cannot open $type2spNAminusFile for writing: $OS_ERROR";
for my $sp ( @type2spNAminusSpp )
{
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
close OUT;


my $nType1NAspp        = @type1NAspp;
my $pType1NAspp        = sprintf("%.1f%%", 100.0 * $nType1NAspp / $nSpp);
my $pType1NAseqs        = sprintf("%.1f%%", 100.0 * $nType1NAseqs / $nSeqs);

my $nType1spNAminusSpp = @type1spNAminusSpp;
my $pType1spNAminusSpp = sprintf("%.1f%%", 100.0 * $nType1spNAminusSpp / $nSpp);
my $pType1spNAminusSeqs = sprintf("%.1f%%", 100.0 * $nType1spNAminusSeqs / $nSeqs);

my $nType2spNAplusSpp  = @type2spNAplusSpp;
my $pType2spNAplusSpp  = sprintf("%.1f%%", 100.0 * $nType2spNAplusSpp / $nSpp);
my $pType2spNAplusSeqs  = sprintf("%.1f%%", 100.0 * $nType2spNAplusSeqs / $nSeqs);

my $nType2spNAminusSpp = @type2spNAminusSpp;
my $pType2spNAminusSpp = sprintf("%.1f%%", 100.0 * $nType2spNAminusSpp / $nSpp);
my $pType2spNAminusSeqs = sprintf("%.1f%%", 100.0 * $nType2spNAminusSeqs / $nSeqs);

my $pType1spNAspp      = sprintf("%.1f%%", 100.0 * $nType1spNAspp / $nSpp);
my $pType1spNAseqs     = sprintf("%.1f%%", 100.0 * $nType1spNAseqs / $nSeqs);

my $pType1spNAplusSpp  = sprintf("%.1f%%", 100.0 * $nType1spNAplusSpp / $nSpp);
my $pType1spNAplusSeqs = sprintf("%.1f%%", 100.0 * $nType1spNAplusSeqs / $nSeqs);

$nSeqs = commify($nSeqs);

print "\n\n\tNo. of sequences:           $nSeqs\n";
print     "\tNo. of species:             $nSpp\n";
print     "\tNo. of type 1spNA species:  $nType1spNAspp ($pType1spNAspp  cov: $pType1spNAseqs)\n";
print     "\tNo. of type 1spNA+ species: $nType1spNAplusSpp ($pType1spNAplusSpp cov: $pType1spNAplusSeqs)\n";
print     "\tNo. of type 1NA species:    $nType1NAspp ($pType1NAspp cov: $pType1NAseqs)\n";
print     "\tNo. of type 1spNA- species: $nType1spNAminusSpp ($pType1spNAminusSpp cov: $pType1spNAminusSeqs)\n";
print     "\tNo. of type 2spNA+ species: $nType2spNAplusSpp ($pType2spNAplusSpp cov: $pType2spNAplusSeqs)\n";
print     "\tNo. of type 2spNA- species: $nType2spNAminusSpp ($pType2spNAminusSpp cov: $pType2spNAminusSeqs)\n";

my @suffixes = (".txt");

print "\n\tType 1spNA  species table written to " . basename($type1spNAfileSpp, @suffixes) . ".txt\n";
print   "\tType 1spNA+ species table written to " . basename($type1spNAplusFileSpp, @suffixes) . ".txt\n";

print "\n\tType 1spNA  species records written to " . basename($type1spNAfile, @suffixes) . ".txt\n";
print   "\tType 1spNA+ species records written to " . basename($type1spNAplusFile, @suffixes) . ".txt\n";
print   "\tType 1NA    species records written to " . basename($type1NAfile, @suffixes) . ".txt\n";
print   "\tType 1spNA- species records written to " . basename($type1spNAminusFile, @suffixes) . ".txt\n";
print   "\tType 2spNA+ species records written to " . basename($type2spNAplusFile, @suffixes) . ".txt\n";
print   "\tType 2spNA- species records written to " . basename($type2spNAminusFile, @suffixes) . ".txt\n\n";

####################################################################
##                               SUBS
####################################################################

## Parsing file with the following three columns
## Generating a hash table <sp> => <sp's phylo-group>

## Akkermansia_muciniphila	2413	Verrucomicrobia_V3V4
## Akkermansia_sp_5	14	Verrucomicrobia_V3V4
## Akkermansia_sp_7	13	Verrucomicrobia_V3V4
## Coraliomargarita_akajimensis	2	Verrucomicrobia_V3V4
## Akkermansia_sp_4	1	Verrucomicrobia_V3V4
## Fibrobacter_sp	51504	phyla_lessthen_1k_wOG_V3V4

sub parseSpTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($sp, $size, $phGr) = split /\s+/,$_;
    $tbl{$sp} = $phGr;
  }
  close IN;

  return %tbl;
}

# difference of two arrays
sub diff{

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

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTbl{

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
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  #print "\n";
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTblToFile{

  my ($rTbl, $rSub, $fh) = @_; # the second argument is a subarray of the keys of the table

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
    print $fh "$_$pad" . $rTbl->{$_} . "\n";
  }
  print $fh "\n";
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify {
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
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

##
## parse 4 column taxon table
##

# file format

#                                                tx   pp                        phGr
# 1642.V1_0                     Lactobacillus_iners 0.93     Firmicutes_group_6_V3V4
# 0980.V2_1                     Lactobacillus_iners 0.97     Firmicutes_group_6_V3V4
# 1670.V2_2 Lactobacillus_crispatus_kefiranofaciens 0.98     Firmicutes_group_6_V3V4
# 0711.V3_3                       Atopobium_vaginae 0.56 Actinobacteria_group_0_V3V4
# 1149.V1_4                     Lactobacillus_iners 0.94     Firmicutes_group_6_V3V4
# 1386.V1_5                     Prevotella_buccalis 0.85  Bacteroidetes_group_2_V3V4
# ...
sub readQtxTbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readQtxTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %spSize;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    next if /^$/;
    chomp;
    my ($id, $sp, $pp, $gr) = split /\s+/,$_;
    $spSize{$sp}++;
  }
  close IN;

  return %spSize;
}

exit 0;
