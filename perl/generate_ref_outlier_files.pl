#!/usr/bin/env perl

=head1 NAME

  generate_ref_outlier_files.pl

=head1 DESCRIPTION

  Given refParentName.txt and ref_outlier_list.txt paths
  parse the last and write seqIDs to the appropriate locations

=head1 SYNOPSIS

  generate_ref_outlier_files.pl [Options]

=head1 OPTIONS

=over

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

  generate_ref_outlier_files.pl --debug

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "verbatim|v"      => \my $verbatim,
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



####################################################################
##                               MAIN
####################################################################

my $pFile = "/Users/pgajer/devel/MCextras/data/refParentName.txt";
my %pTbl = readTbl($pFile);

my $rFile = "/Users/pgajer/devel/MCextras/data/ref_outliers_list.txt";
if ( ! -e $rFile )
{
  warn "\n\n\tERROR: $rFile does not exist";
  print "\n\n";
  exit 1;
}

my $mainDir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";
my @ids;
my %bTbl;
open IN, "$rFile" or die "Cannot open $rFile for reading: $OS_ERROR\n";
foreach (<IN>)
{
  next if /^#/;
  next if /^\s+$/;

  chomp;

  if (/^S00/)
  {
    my ($id, $pp) = split /\s+/, $_;
    push @ids, $id;
  }
  else
  {
    my $prefix = $_;
    #print "line: $_\n";
    $prefix =~ s/_dir\/bad_seqIDs\.pp//;
    #print "prefix: $prefix\n";
    my $dir;
    if (exists $pTbl{$prefix})
    {
      $dir = $pTbl{$prefix};
      my $key = $mainDir . $dir . "/$_";
      push @{$bTbl{$key}}, @ids;
      @ids = ();
    }
    else
    {
      warn "\n\n\tpTbl{$prefix} does not exit";
      print "\n\n";
      exit 1;
    }
  }
}
close IN;

if ($debug)
{
  print "init bTbl:\n";
  printArrayValuedTblCounts(\%bTbl);
  print "\n";

  # print "bTbl:\n";
  # printArrayValuedTbl(\%bTbl);
  # print "\n\n";
}



# $ grep -c S00 /Users/pgajer/devel/MCextras/data/ref_outliers_list.txt
# 826

# count number of elements in all arrays of bTbl
if (0)
{
  my $n = 0;
  for (keys %bTbl)
  {
    $n += @{$bTbl{$_}};
  }
  print "n: $n\n\n";
}

## Reading outlier sequences from another file generated in merging_MC_models.R
## here max.sib for species with max.sib > min.ref has been identified for deletion

$rFile = "/Users/pgajer/devel/MCextras/data/ref_outliers_list2.txt";
if ( ! -e $rFile )
{
  warn "\n\n\tERROR: $rFile does not exist";
  print "\n\n";
  exit 1;
}

open IN, "$rFile" or die "Cannot open $rFile for reading: $OS_ERROR\n";
foreach (<IN>)
{
  next if /^#/;
  next if /^\s+$/;

  chomp;

  if (/^S00/)
  {
    my ($id, $pp) = split /\s+/, $_;
    push @ids, $id;
  }
  else
  {
    my $prefix = $_;
    #print "line: $_\n";
    $prefix =~ s/_dir\/bad_seqIDs\.pp//;
    #print "prefix: $prefix\n";
    my $dir;
    if (exists $pTbl{$prefix})
    {
      $dir = $pTbl{$prefix};
      my $key = $mainDir . $dir . "/$_";
      push @{$bTbl{$key}}, @ids;
      @ids = ();
    }
    else
    {
      warn "\n\n\tpTbl{$prefix} does not exit";
      print "\n\n";
      exit 1;
    }
  }
}
close IN;

if ($debug)
{
  print "\n\nAfter update 1 bTbl:\n";
  printArrayValuedTblCounts(\%bTbl);
  print "\n";

  # print "bTbl:\n";
  # printArrayValuedTbl(\%bTbl);
  # print "\n\n";
}


## ref_outliers_list3.txt was generated in modeling_ref_pp_p2.R

$rFile = "/Users/pgajer/devel/MCextras/data/ref_outliers_list3.txt";
if ( ! -e $rFile )
{
  warn "\n\n\tERROR: $rFile does not exist";
  print "\n\n";
  exit 1;
}

open IN, "$rFile" or die "Cannot open $rFile for reading: $OS_ERROR\n";
foreach (<IN>)
{
  next if /^#/;
  next if /^\s+$/;

  chomp;

  if (/^S00/)
  {
    my ($id, $pp) = split /\s+/, $_;
    push @ids, $id;
  }
  else
  {
    my $prefix = $_;
    #print "line: $_\n";
    $prefix =~ s/_dir\/bad_seqIDs\.pp//;
    #print "prefix: $prefix\n";
    my $dir;
    if (exists $pTbl{$prefix})
    {
      $dir = $pTbl{$prefix};
      my $key = $mainDir . $dir . "/$_";
      push @{$bTbl{$key}}, @ids;
      @ids = ();
    }
    else
    {
      warn "\n\n\tpTbl{$prefix} does not exit";
      print "\n\n";
      exit 1;
    }
  }
}
close IN;


if ($debug)
{
  print "\n\nAfter update 2 bTbl:\n";
  printArrayValuedTblCounts(\%bTbl);
  print "\n";

  # print "bTbl:\n";
  # printArrayValuedTbl(\%bTbl);
  # print "\n\n";
}


## ref_outliers_list4.txt was generated in modeling_ref_pp_p2.R

$rFile = "/Users/pgajer/devel/MCextras/data/ref_outliers_list4.txt";
if ( ! -e $rFile )
{
  warn "\n\n\tERROR: $rFile does not exist";
  print "\n\n";
  exit 1;
}

open IN, "$rFile" or die "Cannot open $rFile for reading: $OS_ERROR\n";
foreach (<IN>)
{
  next if /^#/;
  next if /^\s+$/;

  chomp;

  if (/^S00/)
  {
    my ($id, $pp) = split /\s+/, $_;
    push @ids, $id;
  }
  else
  {
    my $prefix = $_;
    #print "line: $_\n";
    $prefix =~ s/_dir\/bad_seqIDs\.pp//;
    #print "prefix: $prefix\n";
    my $dir;
    if (exists $pTbl{$prefix})
    {
      $dir = $pTbl{$prefix};
      my $key = $mainDir . $dir . "/$_";
      push @{$bTbl{$key}}, @ids;
      @ids = ();
    }
    else
    {
      warn "\n\n\tpTbl{$prefix} does not exit";
      print "\n\n";
      exit 1;
    }
  }
}
close IN;

if ($debug)
{
  print "\n\nAfter update 3 bTbl:\n";
  printArrayValuedTblCounts(\%bTbl);
  print "\n";

  # print "bTbl:\n";
  # printArrayValuedTbl(\%bTbl);
  # print "\n\n";
}


for my $file (keys %bTbl)
{
  print "Writing to $file\n";
  open OUT, ">$file" or die "Cannot open $file for writing: $OS_ERROR\n";
  for my $e (@{$bTbl{$file}})
  {
    print OUT "$e\n";
  }
  close OUT;
}



####################################################################
##                               SUBS
####################################################################

sub printArrayValuedTbl{

  my $rTbl = shift;

  for (keys %{$rTbl})
  {
    print "$_\n";
    for my $e (@{$rTbl->{$_}})
    {
      print "\t$e\n";
    }
  }
  print "\n";
}

sub printArrayValuedTblCounts{

  my $rTbl = shift;

  for (sort keys %{$rTbl})
  {
    print "$_\t" . @{$rTbl->{$_}} . "\n";
  }
  print "\n";
}


sub writeArrayValuedTbl{

  my ($rTbl, $outFile) = @_;

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
  for (keys %{$rTbl})
  {
    print "$_\n";
    for my $e (@{$rTbl->{$_}})
    {
      print OUT "\t$e\n";
    }
  }
  close OUT;
}


# read two column table; create a table that assigns
# elements of the first column to the second column
sub readTbl
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



exit 0;
