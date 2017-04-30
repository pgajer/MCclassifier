#!/usr/bin/env perl

=head1 NAME

  get_master_mm_report.pl

=head1 DESCRIPTION

  mm_validate_pecan.pl and mm_validate2_pecan.pl generated vicut_cltrs_report.txt files.
  Each species constitutes a record, which is of the form

  ------------------------------------------------

  Stenotrophomonas_maltophilia

  n:     13
  n(nr): 6


  Cluster c71 (457)

 Stenotrophomonas_maltophilia: 451
 NA:                           6

 Coverage: 100.0% (13 out of 13 seq's)
  Size ranks: 1 2 3 4 5 6
  Size %'s: 61.54 7.69 7.69 7.69 7.69 7.69
  pp's: 0.94 0.67 0.56 0.95 0.92 0.85
  range(pp): [0.56, 0.95]


  Thus the record start with the hyphen line '------------------------------------------------'
  and goes to the next line (or the end of file).

  The purpose of this script is to detect duplicate records and also incomplete records.

  All together there should be 1022 species/records in all reports.

  The script is going to produce a master report with all records.

  Also missing_species.txt file fill be produced.

  Maybe later other types of output will be implemented.

=head1 SYNOPSIS

  get_master_mm_report.pl

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

  get_master_mm_report.pl

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

## parsing table of species detected in the M&M dataset
my $sppTblFile = $mmDir . "/mm_uq_spp_report_no_controls.txt";
my %spToPhGr = parseSpTbl($sppTblFile); # spToPhGr: <sp> => <sp's phylo-group>

my @allSpp = keys %spToPhGr;
my $nAllSpp = @allSpp;

my %rTbl; # <species name> => record with cluster info text
my $sp;
my $startRec = 0;
my %recComplete; # $recComplete{$sp} = 1 if the record contains 'Cluster' and 'range(pp) strings in it. Otherwise recComplete is not defined on $sp
my $foundCluster = 0;
my $inPhGrHeader = 0;

opendir(DIR, $mmDir) or die $!;
while (my $file = readdir(DIR))
{
  next if ($file !~ /^mm_validate_reports_dir/);
  #print "\nFround file: $file\n";
  my $dir = $file;
  $file .= "/vicut_cltrs_report.txt";
  #print "Redefined to $file\n"; exit;

  if ( ! -e $file )
  {
    my $cmd = "rm -rf $dir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    next;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for (<IN>)
  {
    if (/^------------------------------------------------/)
    {
      $startRec = 1;
      $foundCluster = 0;
    }
    elsif ( !/^$/ && $startRec )
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
      $sp = $_;
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
      elsif (/^=======/ && !$inPhGrHeader)
      {
	$inPhGrHeader = 1;
      }
      elsif (/^=======/ && $inPhGrHeader)
      {
	$inPhGrHeader = 0;
      }

      if ( defined $sp && !/^=======/ && !$inPhGrHeader )
      {
	$rTbl{$sp} .= $_;
      }
    }
  }
  close IN;
}
closedir(DIR);

##
## Deleting control (non-vaginal) species
##
my @controlsOnlySpp;
for my $sp (keys %rTbl)
{
  if ( !exists $spToPhGr{$sp} )
  {
    delete $rTbl{$sp};
    if (exists $recComplete{$sp} )
    {
      delete $recComplete{$sp};
    }
    push @controlsOnlySpp, $sp;
  }
}

if (@controlsOnlySpp)
{
  print "\n\nDeleted the following controls only species\n\n";
  printArray(\@controlsOnlySpp);
}
print "\n\n";

##
## Summary
##
my @recSpp = keys %rTbl;
my $nRecSpp = @recSpp;
my @complRecSpp = keys %recComplete;
my $nComplRecSpp = @complRecSpp;
my @incompRecSpp = diff(\@recSpp, \@complRecSpp);
my $nIncompRecSpp = @incompRecSpp;

if (@incompRecSpp)
{
  printArray(\@incompRecSpp, "\nSpecies with incomplete records\n");
}

if (0)
{
  print "\nDetected records\n";
  for my $sp (@complRecSpp)
  {
    print "------------------------------------------------\n";
    print "\n$sp\n";
    print $rTbl{$sp};
  }
  print "\n";
}

my $masterFile = $mmDir . "/master_mm_validate_report.txt";
open OUT, ">$masterFile" or die "Cannot open $masterFile for writing: $OS_ERROR";
for my $sp (@complRecSpp)
{
  print OUT "------------------------------------------------\n";
  print OUT "\n$sp\n";
  print OUT $rTbl{$sp};
}
print OUT "\n";
close OUT;


## Generating a list of species missing record
my @spMissingRec = diff(\@allSpp, \@complRecSpp);

my $missingRecFile = $mmDir . "/spp_missing_record.txt";
open OUT, ">$missingRecFile" or die "Cannot open $missingRecFile for writing: $OS_ERROR\n";
for my $sp (@spMissingRec)
{
  print OUT "$sp\t-1\t" . $spToPhGr{$sp} . "\n";
}
close OUT;

print "\n\n\tNo. of species detected in the M&M project: $nAllSpp\n";
print     "\tNo. of all detected species:                $nRecSpp\n";
print     "\tNo. of complete record species:             $nComplRecSpp\n";
print     "\tNo. of species with incomplete record:      $nIncompRecSpp\n";
print     "\tNo. of species with no record:              " . @spMissingRec . "\n";

print "\n\tComplete records written to $masterFile\n";
print "\tSpecies with missing record written to $missingRecFile\n\n";


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

exit 0;
