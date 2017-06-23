#!/usr/bin/env perl

=head1 NAME

  clean_master_lineage_v2.pl

=head1 DESCRIPTION

  This script modifies species names so that they are of the form Genus_speciesname
  where Genus is capitilized and speciesname is a lowercase

  '/', '-' and '_' characters (except the '_' between genus and speciesname) are changed to '|'

  _gp1/2/3 changes to |gp1/2/3

  Escherichia-Shigella_sp => Escherichia|Shigella_sp

  sp.suffix => sp

  Genus name and prefix of species name have to be the same

=head1 SYNOPSIS

  clean_master_lineage_v2.pl -i <lineage file> -o <output lineage file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input lineage file.

=item B<--output-file, -o>
  Output lineage file.

=item B<--verbose, -v>
  Prints content of some output files. Default value: 5000.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/RDP/r2.2

  clean_master_lineage_v2.pl -i rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis_nr.lineage -o rdp_Bacteria_curated.lineage

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

GetOptions(
  "input-file|i=s"      => \my $inFile,
  "output-file|o=s"      => \my $outFile,
  "quiet"               => \my $quiet,
  "igs"                 => \my $igs,
  "verbose|v"           => \my $verbose,
  "dry-run"             => \my $dryRun,
  "debug"               => \my $debug,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$inFile )
{
  print "ERROR: Missing input lineage file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outFile )
{
  print "ERROR: Missing output file name\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $debugStr = "";
my $quietStr = "--quiet";
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

####################################################################
##                               MAIN
####################################################################

print "--- Parsing lineage table\n";
my $suffix;
my %liTbl;
open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
for ( <IN> )
{
  chomp;
  my ($id, $li) = split /\s+/;

  my @f = split ";", $li;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;

  my $origSp = $sp;
  #$sp = "AAA-BBB_gp1_sp.sldur-23X";
  #$sp = "Curtobacterium_MG-2011-84-GV";
  #$sp = "Mycobacterium_n";
  #$sp = "Granulicatella_para|adiacens";
  #$sp = "Selenomonas_flueggei|like";
  #print "sp BEFORE: $sp\n";
  $sp =~ s/\//\|/g;
  $sp =~ s/\-/\|/g;
  $sp =~ s/_gp\d//;
  $sp =~ s/\..+//;
  $sp =~ s/_[[:upper:]].+/_sp/;
  $sp =~ s/_\w$/_sp/;
  $sp =~ s/Ruminococcus2/Ruminococcus/;
  $sp =~ s/Armatimonas\|Armatimonadetes/Armatimonas/;
  $sp =~ s/Chthonomonas\|Armatimonadetes/Chthonomonas/;
  $sp =~ s/Escherichia\|Shigella/EscherichiaShigella/;
  $sp =~ s/_para\|/_/;
  $sp =~ s/\|like$//;
  $sp =~ s/lxb\|14/lxb14/;
  $sp =~ s/lxb\|3/lxb3/;
  #print "sp AFTER: $sp\n";  exit 1;

  if ( $sp !~ /[[:upper:]][\||\w]+_[[:lower:]]\w+/ && $sp !~ /BVAB/ )
  {
    warn "\n\n\tERROR: incorrect format sp: $sp";
    print "origSp: $origSp\n\n";
    exit 1;
  }

  if ( $sp !~ /BVAB/ )
  {
    ($ge, $suffix) = split "_", $sp;
  }

  if ( $ge eq "Ruminococcus" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }

  if ( $ge eq $fa )
  {
    next;
  }

  $li = join ";", (@f, $fa, $ge, $sp);
  $liTbl{$id} = $li;
  #print "$id\t$li\n"; exit 1;
}
close IN;

print "--- Checking uniquness of parents\n";
if ( each_tx_has_unique_parent(\%liTbl) )
{
  exit 1;
}

write_tbl( \%liTbl, $outFile );


####################################################################
##                               SUBS
####################################################################

## check if each node of the lineage structure has only one parent
sub each_tx_has_unique_parent
{
  my $r = shift;
  my %liTbl = %{$r};

  my  %prt;
  for my $id ( keys %liTbl )
  {
    my $lineage = $liTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $prt{$sp}{$ge}++;
    $prt{$ge}{$fa}++;
    $prt{$fa}{$or}++;
    $prt{$or}{$cl}++;
    $prt{$cl}{$ph}++;
  }

  my $ret = 0;
  for my $tx (keys %prt)
  {
    my $nPrts = keys %{$prt{$tx}};
    if ( $nPrts > 1 )
    {
      $ret = 1;
      warn "\n\n\tERROR: $tx has more than one parent";
      print "\n\t$tx parents\n";
      for ( keys %{$prt{$tx}} )
      {
        print "\t\t$_\t" . $prt{$tx}{$_} . "\n";
      }
      print "\n\n";

      #print "Attempting to fix the problem. Going with more abundant lineage\n";
      my $tmpLiFile = "tmp.lineage";
      write_tbl( \%liTbl, $tmpLiFile );
      print "\tLineage file written to $tmpLiFile\n\n";
    }
  }

  return $ret;
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

exit 0;
