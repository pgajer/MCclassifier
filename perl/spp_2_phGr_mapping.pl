#!/usr/bin/env perl

=head1 NAME

  spp_2_phGr_mapping.pl

=head1 DESCRIPTION

  Generating a two column table: <species> <phylo-group>

=head1 SYNOPSIS

  spp_2_phGr_mapping.pl -i <phGr list file> -o <sp => phGr table file> [Options]

=head1 OPTIONS

=over

=item B<--phGr-file, -i>
  Input fasta file.

=item B<--out-file, -o>
  Output file.

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

  spp_2_phGr_mapping.pl -i V3V4_dirs -o V3V4_spp_2_phGr_tbl.txt

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
  "phGr-file|i=s"       => \my $phGrFile,
  "out-file|o=s"        => \my $outFile,
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
elsif ( !$outFile )
{
  print "\n\nERROR: Missing output file\n\n";
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

####################################################################
##                               MAIN
####################################################################

print "--- Building species => phGr table\n";
build_sp_2_phGr_tbl( $phGrFile, $outFile );

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

# write hash table to a file
sub write_tbl
{
  my ($rTbl, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

sub build_sp_2_phGr_tbl
{
  my ( $phGrFile, $outFile ) = @_;

  unlink( $outFile );

  my @dirs = read_array( $phGrFile );

  my %spTbl;
  for my $dir ( @dirs )
  {
    my $phGr = $dir;
    $phGr =~ s/_dir//;
    my $file = $dir . "/$phGr" . "_final.lineage";
    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach ( <IN> )
    {
      chomp;
      my ($id, $li) = split /\s+/;
      my @f = split ";", $li;
      my $sp = pop @f;
      $spTbl{$sp} = $phGr; # this assumes that each species is in only one phGr,
                           # which is true for most species but a few ones that
                           # were elimineated from the models
    }
    close IN;
  }

  my @delSpp = ("Geodermatophilus_sp_1", "Inhella_sp", "Hyphomonas_sp",
  "EscherichiaShigella_sp", "Erwinia_sp_1", "Alteromonas_sp_1",
  "Skermania_piniformis", "Alteromonas_sp_2", "Leeuwenhoekiella_sp",
  "Geodermatophilus_sp_2", "Succinimonas_sp", "Rothia_sp", "Erwinia_sp_2");

  delete @spTbl{$delSpp};

  write_tbl( \%spTbl, $outFile );
}

exit 0;
