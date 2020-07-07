#!/usr/bin/perl -w

use strict;

my $gbk_file = shift;
my $locus_prefix = shift;

die "Usage: fasta_splitter.pl <genbank file> <locus prefix>" unless $gbk_file;


$locus_prefix = "" unless $locus_prefix;

my $text = slurp_file($gbk_file);
$text =~ s/\/\/$//s; # remove last field separator

my @records = split '\/\/', $text;

my $offset = 0;

my $faa_file = $gbk_file;
my $fna_file = $gbk_file;
my $ptt_file = $gbk_file;
my $rnt_file = $gbk_file;
my $contig_file = $gbk_file;
$faa_file =~ s/gbk$/faa/;
$fna_file =~ s/gbk$/fna/;
$ptt_file =~ s/gbk$/ptt/;
$rnt_file =~ s/gbk$/rnt/;
$contig_file =~ s/gbk$/contig/;

open (FAA, ">$faa_file") or die "can't open $faa_file: $!\n";
open (FNA, ">$fna_file") or die "can't open $fna_file: $!\n";
open (PTT, ">$ptt_file") or die "can't open $ptt_file: $!\n";
open (RNT, ">$rnt_file") or die "can't open $rnt_file: $!\n";

my $write_contig_info = 0;
if (scalar @records > 1) {
	open (CONTIG, ">$contig_file") or die "can't open $contig_file: $!\n";
	$write_contig_info = 1;
}


my $first = 1;
my $gene_count = 0;
my $sequence_offset = 0;
for my $r (@records) {
	my $p = parse_record($r);
	

	if ($first) {
		my $org = $p->{ORGANISM};
		init_ptt($org);
		init_rnt($org);

		if (scalar @records > 1) {
			init_contig($org);
		}

		$first = 0;
	}

	write_features($p, \$gene_count,$sequence_offset);
	write_nucleotides($p, \$gene_count);
	if (scalar @records > 1) {
		printf CONTIG "$p->{LOCUS}\t%d\t%d\t+\t%d\n", $sequence_offset+1, $sequence_offset+$p->{sequence_length},$p->{sequence_length};
	}
	$sequence_offset += $p->{sequence_length};
}

close FAA;
close FNA;


#die "have " . scalar @records . " records\n";
sub init_ptt {
	my $organism = shift;
	print PTT "$organism, whole genome shotgun sequencing project\n";
	print PTT "Product Name\tStart\tEnd\tStrand\tLength\tGi\tGeneID\tLocus\tLocus_tag\tCOG(s)\n";
}

sub init_rnt {
	my $organism = shift;
	print RNT "$organism, whole genome shotgun sequencing project\n";
	print RNT "Product Name\tStart\tEnd\tStrand\tLength\tGeneID\tLocus\tLocus_tag\n";
}

sub init_contig {
	my $organism = shift;
	print CONTIG "$organism, whole genome shotgun sequencing project\n";
	print CONTIG "Accession\tStart\tEnd\tStrand\tLength\n";
}

sub write_nucleotides {
	my $parsed_record = shift;
	print FNA to_fasta("gi|0|ref|$parsed_record->{LOCUS}| $parsed_record->{ORGANISM}",$parsed_record->{sequence});
}

sub write_features {
	my $parsed_record = shift;
	my $gene_count = shift;
	my $seq_offset = shift;

	my $features = $parsed_record->{features};

	for my $f (@$features) {
#		my @k = keys %$f;
#		print "keys @k\n";
		$$gene_count++;
		my $start = $f->{start} + $seq_offset;
		my $stop = $f->{stop} + $seq_offset;
		if ($f->{type} eq 'CDS') {
			print FAA to_fasta("gi|$$gene_count|ref|$f->{locus_tag}| $f->{product}",$f->{translation});
			print PTT "$f->{product}\t$start\t$stop\t$f->{direction}\t0\t$$gene_count\t$$gene_count\t0\t$f->{locus_tag}\t-\n";
		}
		else {
			my $len = ($stop-$start) + 1;
			#my @k = keys %$f; 
#			print "have @k\n";
			#print RNT "$f->{product}\t$start\t$stop\t$f->{direction}\t$len\t$$gene_count\t0\t$f->{locus_tag}\n";
			my $locus_tag = $locus_prefix . "r" . $$gene_count;
			print RNT "$f->{product}\t$start\t$stop\t$f->{direction}\t$len\t$$gene_count\t0\t$locus_tag\n" if $f->{product};
		}
		
	}

}

sub to_fasta {
        my ($seqName, $seq, $len) = @_;

        $len = 80 unless $len;

        my $formatted_seq = ">$seqName\n";
        while (my $chunk = substr($seq, 0, $len, "")) {
                $formatted_seq .= "$chunk\n";
        }

        return $formatted_seq;
}

sub parse_record {
	my $r = shift;

	my @lines = split '\n', $r;
	my %parse_info;
	my $in_sequence = 0;
	my $in_feature = 0;
	my $in_translate=0;
	my $in_product=0;

	my %feature;
	
	for (@lines) {
		if ($in_sequence) {
			s/\d+//;
			s/\s+//g;
			$parse_info{sequence} .= $_;
		}
		elsif ($in_feature) {
			if (/\/locus_tag="(.*)"/) {
#				die "$_ tag $1\n";
				$feature{locus_tag} =  $1;
			}
			elsif (/^\s+\/product=/) {
				s/^\s+\/product=\"//;
				$in_product = 1;
			}
			elsif (/^\s+\/translation=/) {
				s/^\s+\/translation=\"//;
				$in_translate=1;
			}
		
			if ($in_product) {
				s/\s+/ /;
				$feature{product} .= $_;
				if (/\"$/) {
					$feature{product} =~ s/\"$//;
					$in_product=0;
#					print "have product $feature{product}\n";
				}
			}
			elsif ($in_translate) {
				s/\s+//;
				$feature{translation} .= $_;
#				print "trans $feature{translation}\n";
				if (/\"$/) {
					$feature{translation} =~ s/\"$//;
					$in_translate=0;
					$in_feature=0;
#					die "$feature{translation}\n";
				}
			}
#			elsif (/\/product="(.*)"/) {
		}
		
		if (/^LOCUS/) {
			my @pieces = split /\s+/;
#			die "have pieces @pieces\t$pieces[1]\n";
			$parse_info{LOCUS} = $pieces[1];
		}
		elsif (/^SOURCE\s+(.*)/) {
#			die "have organism $1\n";
			$parse_info{ORGANISM} = $1;
		}
		elsif (/^\s+(CDS)/ || /^\s+(tRNA)/ || /^\s+(rRNA)/) {
			my $type = $1;
			$in_feature = 1;	
			my ($start,$stop) = /(\d+)\.\.(\d+)/;
			my $direction = "";
			if (/complement/) {
				$direction = "-";
			}
			else {
				$direction = "+";
			}
			if (%feature) {
				my %copy = %feature;

				my @k = keys %copy;
#				die "adding @k\n";
				push @{$parse_info{features}}, \%copy;
				%feature = ();
			}
			$feature{start} = $start;
			$feature{stop} = $stop;
			$feature{direction} = $direction;
			$feature{type} = $type;
#			print;
#			die "\t\t$start\t$stop\t$direction\t$type\n" if $type eq 'tRNA';

		}
		elsif (/^ORIGIN/) {
			if (%feature) {
				my %copy = %feature;	
				push @{$parse_info{features}}, \%copy;
				%feature = ();
			}
			$in_sequence=1;
		}


	}
	$parse_info{sequence_length} = length($parse_info{sequence});
#	die "have " . scalar @lines . " lines\n";
#	die "have $parse_info{LOCUS}\n$parse_info{sequence_length}\n";
	
	return \%parse_info;
}

sub slurp_file {
	my $gbk_file = shift;

	open (IN, $gbk_file) or die "can't open $gbk_file: $!\n";
	my $text="";
	while (<IN>) {
		$text .= $_;
	}

	close IN;

	return $text;
}
