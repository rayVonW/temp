#!/usr/bin/env perl

use strict;
use warnings;
use Text::CSV;
use Getopt::Long;

my $usage = <<END;
count_tags.pl --seq_dir DIRECTORY --barcodes BARCODE_LIST [--by_tag --five_p_seq STRING --three_p_seq STRING]

Counts barcode tags from a directory of fastq files. Sample names are extracted from the 
fastq file names and a table is generated with counts for each of the barcodes in the
barcodes list (onyl if there is atleast one occurrence) for each of the samples (libraries).

required:
--seq_dir      directory of fastq files - one per sample
--barcodes     file of mappings of genes to barcodes format:
               gene_id,barcode
               PBANKA_000010,ttcgccgggcc
               PBANKA_000030,caggcaatcgg

optional:
--by_tag       Usually, each tag corresponds to one gene and counts are combined if 
               we have more than one tag for the same gene. With this option, counts
               are separated by tag even if they map to the same gene.
--five_p_seq
--three_p_seq  These options allow setting the sequence context around the tag if it
               isn't the standard BA primer/R2-amp97 amplicon.
               They must be 3' and 5' of the barcode in the orientation as given #
               in the barcodes file.
               Only the 5nt directly adjacant to the barcode are used by default.
--ignore_missing_tag ignore is if there is a row in the barcodes file where there is
                     no entry for the barcode. If not set, this is an error.
                     This can be used to employ design data files as barcodes lists
                     which may have empty barcodes in some cases.

--nomatch_out_file FILE  if given, write (in fastq) the sequences that don't match any 
                         oligo (or match more than once) to this file

END

my $by_tag ; 

# the sequence of the tag amplification primer
my $ba_primer = 'GTAATTCGTGCGCGTCAG';
# sequence from cassette primer R2 to primer binding site for arg97
# this is the same for all constructs
my $r2_to_amp97 = 'CCGCCTACTGCGACTATAGAGATATCAACCACTTTGTACAAGAAAGCTGGGTGGTACCCATCGAAATTGAAGG';
my $input_dir;
my $barcode_table;
my $ignore_missing_tag;
my $nomatch_out_file;
my $iswalkupservice;

GetOptions(
  'seq_dir=s' => \$input_dir,
  'ignore_missing_tag' => \$ignore_missing_tag,
  'barcodes=s' => \$barcode_table,
  'by_tag' => \$by_tag,
  'five_p_seq=s' => \$ba_primer,
  'three_p_seq=s' => \$r2_to_amp97,
  'nomatch_out_file=s' => \$nomatch_out_file,
  'help|h' => sub{ die $usage },
);

die "no directory of sequences (--seq_dir) given\n$usage" unless $input_dir; 
die "$input_dir not readable\n" unless -r $input_dir;
die "$input_dir not a directory\n" unless -d $input_dir;
die "no barcodes table (--barcodes) given\n$usage" unless $barcode_table;
die "$barcode_table is not readable" unless -r $barcode_table && -f $barcode_table;

# extract the relevant bits of the context sequence
die "five_p_seq (BA-primer) must be > 5nt long" unless $ba_primer && length($ba_primer)>=5;
die "three_p_seq (R2-amp97) must be > 5nt long" unless $r2_to_amp97 && length($r2_to_amp97)>=5;
my $ba_primer_end = substr( $ba_primer, -5);
my $ba_primer_end_rc = revcom( $ba_primer_end );
my $r2_start = substr( $r2_to_amp97, 0, 5);
my $r2_start_rc = revcom( $r2_start );

my $nomatch_out_fh;
if ( $nomatch_out_file ){
  open $nomatch_out_fh, '>', $nomatch_out_file or die "could not open $nomatch_out_file for writing\n";
}

my $csv = Text::CSV->new ( { binary => 1 } ) or die Text::CSV->error_diag;
open my $fh, '<', $barcode_table or die "could not open $barcode_table\n";

my $colref = $csv->getline( $fh );
$csv->column_names( @$colref );

# read barcodes
my %barcode2genes;
while (my $colref = $csv->getline_hr( $fh )){
  my $gene_id = $colref->{'gene id'} || $colref->{'gene_id'} or die "could not read gene id";
  my $tag = $colref->{'tag'} || $colref->{'Barcode'} || $colref->{barcode};
  die "could not read barcode tag" if ! $tag and ! $ignore_missing_tag;
  next if $tag=~/none/i;
  if (length($tag) <8 or length($tag)>16){
    die "found a barcode outside allowed length range (8-16): $tag, length:".length($tag)."\n";
  }
  $barcode2genes{lc $tag} = $gene_id;
}
close $fh;

# get fastq files
my @fastq_files;
opendir(my $dh, $input_dir) || die "can't opendir $input_dir: $!";
@fastq_files = grep { /\.fastq$/ && -f "$input_dir/$_" } readdir($dh);
closedir $dh;

my %counts;
my %samples;

foreach my $fastq_file (@fastq_files){
  warn "reading file $fastq_file\n";
  my $sample;
    ($sample) = ($fastq_file=~/^(.*?\d+).fastq/);
  if ( !$sample){
    die "could not parse fastq file name '$fastq_file'"
  }
  $samples{$sample} = 1;
  open (my $fh, '<',  "$input_dir/$fastq_file") or die "could not open $input_dir/$fastq_file";
  while (my ( $read_id, $seq, $read_id2, $qual) = next_seq_from_fastq($fh)){
    my @found_tags_fwd = ( $seq=~/$ba_primer_end(\w{8,16})$r2_start/ig), 
    my @found_rev_tags = map { revcom( $1  ) }( $seq=~/$r2_start_rc(\w{8,16})$ba_primer_end_rc/ig);
    my @found_tags = (@found_tags_fwd, @found_rev_tags);

    # there could be more than one pattern match (we use very short context sequences)
    # find the ones that actually contain an existing tag. It is an error if there is more than one
    my @existing_found_tags;
    foreach my $tag (@found_tags) {
      push @existing_found_tags, $tag if defined $barcode2genes{ lc $tag };
    }
    if ( @existing_found_tags == 1){
        ++$counts{ lc $existing_found_tags[0] }{$sample};
    } else {
      if ( @existing_found_tags > 1){
        warn "Found ".@existing_found_tags." matching tag in sequence $seq - counting as 'no match'\n" ;
      }
      ++$counts{no_match}{$sample} ;
      warn "NOMATCH $seq\n";
      if ($nomatch_out_fh){
        print $nomatch_out_fh "$read_id\n$seq\n$read_id2\n$qual\n";
      }
    }  
  }
}

close $nomatch_out_fh if $nomatch_out_fh;

# print results
print "\"barcode\",\"gene\",";

print join(',', sort keys %samples) . "\n";

# the primary key will be barcodes/tags
foreach my $pk ( sort mysort keys %counts ){
  my $gene = $barcode2genes{ $pk } ||'';
  print "\"$pk\",\"$gene\"";
  foreach my $sample ( sort keys %samples ){
      my $count = defined $counts{$pk}{$sample}? $counts{$pk}{$sample} : 0;
      print ",$count";
  }
  print "\n";
}

sub mysort{
  if ($a eq 'no_match'){
    return -1 ;
  } elsif ($b eq 'no_match'){
    return 1;
  } else {
    return $a cmp $b ;
  }
}

sub next_seq_from_fastq{
  my $fh = shift;
  my $id = <$fh>;
  return unless $id ;
  chomp( $id );
  chomp( my $seq = <$fh>);
  chomp( my $id2 = <$fh>);
  chomp( my $qual = <$fh>);
  return ( $id, $seq, $id2, $qual );
}

sub revcom {
  my $s = shift;
  return unless $s;
  $s=~tr/AGCTagct/TCGAtcga/;
  $s = reverse( $s );
  return $s;
}
