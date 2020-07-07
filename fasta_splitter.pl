#!/usr/bin/perl -w

use strict;

my $fasta_file=shift;
my $num_jobs = shift;

die "Usage: fasta_splitter.pl <fasta_file> <num_jobs>" unless $fasta_file && $num_jobs;

my $wc = `grep '>' $fasta_file | wc`;
$wc =~ s/^\s+//;
my ($numlines, @junk) = split /\s+/, $wc;

my $num_per_round = int($numlines/$num_jobs)+1;

print "splitting $numlines line $fasta_file into $num_jobs jobs with $num_per_round in each job\n";

my $suffix = 0;
my $line_num=1;

#my $sh_out = "$scarf_file.sh";
#open (SH_OUT, ">$sh_out") or die "can't open $sh_out: $!\n";

#my $file_base = $scarf_file;
open (IN, $fasta_file) or die "can't open $fasta_file: $!\n";
my @temp = split /\//, $fasta_file;
my $last_part = $temp[$#temp];
my $outfile = "$last_part.$suffix";

#die "have temp @temp\t$outfile\n";
open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
while(<IN>) {
	if (/^>/ && 0 == ($line_num++ % $num_per_round)) {
#		print SH_OUT "./seqMap -q $outfile -D dbfiles/ecoli_for_solexa -w 12 | gzip > $outfile.out.gz\n";
		close OUT;	
		$suffix++;
		$outfile = "$last_part.$suffix";
		open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
	}
	print OUT $_;
}
close OUT;
#$outfile = "$scarf_file.$suffix";
#print SH_OUT "./seqMap -q $outfile -D dbfiles/ecoli_for_solexa -w 12 | gzip > $outfile.out.gz\n";

#close SH_OUT;

