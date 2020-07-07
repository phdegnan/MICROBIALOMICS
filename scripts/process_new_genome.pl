#!/usr/bin/perl -w

#    Script: process_new_genome.pl
#    ___________________________________________________________________________
#
#    Version 1.0
#
#    Copyright (C) 2008-2009 Brian Muegge and Jeremiah Faith
#		Patrick Degnan revised last 2014
#
#    http://hamlet.wustl.edu/microbialomics_dev/docs_db/
#
#    About: License
#
#       Licensed under the GNU General Public License
#
#        This program is free software; you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation; either version 2 of the License, or
#        (at your option) any later version.
#
#        This program is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with this program; if not, visit http://www.gnu.org/licenses/gpl.txt
#        or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
#        Boston, MA  02111-1307  USA.
#
# Topic: Getting started
# The annotation pipeline can start with either a set of contigs in FASTA format or an accession 
# to an NCBI refseq genome.
#
# (start code) 
#  perl process_new_genome.pl -i <refseq_accession or fasta_file>
# (end)
#
# Topic: Annotating a complete RefSeq genome
# uses the existing gene and RNA calls and determines the functional annotation of each gene
#
# You'll need to download the *.rnt, *.ptt, *.faa, and *.fna files prior to running process_new_genome.pl, and have them all in the same directory.
# process_new_genome.pl will rename to headers in the *.faa file to correspond to the gene symbol name in the *.ptt file (this makes it easier to 
# parse the results and put them in a database later).  The renamed protein sequences will then be run against all of the annotation programs.  The genes,
# tRNAs, and rRNAs are all taken from the *.ptt and *.rnt files rather than called through the relevant software (i.e. Glimmer, GeneMark, tRNA-Scan, and rnammer).
#
# Topic: Annotating a set of contigs
# concatenates contigs, calls genes and RNAs, and determines the functional annotation of each gene
#
# For contigs, the program will contenate the contigs filtering out smaller ones according to the user-set threshold.   
# Then it will run GeneMark and Glimmer to generate an initial gene set.  GeneMark genes are given preference, and Glimmer genes are 
# added in if they don't overlap any GeneMark gene calls by more than 10%.
#
# Topic: To Do
#      *singalP and psortB to predict protein localization
#
#
# Current Annotations:
# HMMSCAN using TIGRFAM and PFAM
# COG and KEGG (from string database) using PMMER
# subcellular localization using CELLO
# Use INFERNAL to search for RFAM ncRNAs and riboswitches
# Find transcriptional terminators with RNIE

## Converted from Omega (by way of Louise) to Symbiont


use strict;
use Getopt::Long;
use FindBin qw($Bin);
## PHD
# EDIT EXECUTABLE FILE LOCATIONS TO MATCH SYSTEM
#my $BLAST = '/usr/bin/blastall';
#my $BLASTP = '/usr//bin/blastp';##udpated June 2014
my $PHMMER = '/usr/bin/phmmer';##udpated June 2014
#my $HMMPFAM = '/usr/bin/hmm2pfam';
my $HMMSCAN  = '/usr/bin/hmmscan';## udpated June 2014
my $BATCH = '/home/pdegnan/Scripts/split_run_check_combine.pl';##udpated June 2014
my $GLIMMER = '/home/pdegnan/Scripts/g3-from-scratch_w_concat.pl';##udpated June 2014
my $PRODIGAL = '/usr/local/bin/prodigal';## udpated June 2014
my $TRNA_SCAN = '/usr/local/bin/tRNAscan-SE';##udpated June 2014
my $RNAMMER = '/usr/local/bin/rnammer';#!!
#my $RNAMMER = '/home/pdegnan/Software/rnammer-1.2.src/rnammer';# alt used when main not working
my $CELLO = '/home/pdegnan/Software/CELLO_libsvm-2.6/CELLO/CELLO';
my $INFERNAL = '/home/pdegnan/Scripts/multi_cmsearch.pl'; ##udpated June 2014
my $RNIE = '/home/pdegnan/Software/RNIE-master/rnie.pl'; ##udpated June 2014

# EDIT DATABASE LOCATIONS TO MATCH SYSTEM
my $STRING = '/data/DB/Stringv9.1.faa';##udpated June 2014
my $KEGG = '/data/DB/KEGG.faa';##udpated June 2014
#my $TIGRFAM = '/data/DB/TIGRFAMs_9.0/TIGRFAMs_9.0_HMM.LIB';
#my $PFAM = '/data/DB/Pfam_ls'; # currently using version 23; 
# use ls = global domain alignment; fs is for local alignment on protein fragments
my $TIGRFAM14 = '/data/DB/TIGRFAMs_14.0/TIGRFAMs_14.0.HMM';##udpated June 2014
my $PFAM27 = '/data/DB/PFAM_27.0/Pfam-A.hmm'; ##udpated June 2014




#my $GENE_MARK = "/home/jglab/bmuegge/Annotation_Pipeline/bin/genemark_suite_linux/gmsn.pl";
#my $GENE_MARK_HMMP = "/home/jglab/bmuegge/Annotation_Pipeline/bin/genemark_suite_linux/gmhmmp";
#my $GLIMMER = '/home/jglab/faithj/glimmer2.13/run-glimmer2';
#my $TRNA_SCAN = '/home/comp/jglab/faithj/tRNAscan-SE-1.23/tRNAscan-SE';


# http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi  transl_table=11
my $TRANS_TABLE = {
	TTT => 'F', TCT => 'S', TAT => 'Y', TGT => 'C', TTC => 'F', TCC => 'S',
	TAC => 'Y', TGC => 'C', TTA => 'L', TCA => 'S', TAA => '*', TGA => '*',
	TTG => 'L', TCG => 'S', TAG => '*', TGG => 'W', CTT => 'L', CCT => 'P',
	CAT => 'H', CGT => 'R', CTC => 'L', CCC => 'P', CAC => 'H', CGC => 'R',
	CTA => 'L', CCA => 'P', CAA => 'Q', CGA => 'R', CTG => 'L', CCG => 'P',
	CAG => 'Q', CGG => 'R', ATT => 'I', ACT => 'T', AAT => 'N', AGT => 'S',
	ATC => 'I', ACC => 'T', AAC => 'N', AGC => 'S', ATA => 'I', ACA => 'T',
	AAA => 'K', AGA => 'R', ATG => 'M', ACG => 'T', AAG => 'K', AGG => 'R',
	GTT => 'V', GCT => 'A', GAT => 'D', GGT => 'G', GTC => 'V', GCC => 'A',
	GAC => 'D', GGC => 'G', GTA => 'V', GCA => 'A', GAA => 'E', GGA => 'G',
	GTG => 'V', GCG => 'A', GAG => 'E', GGG => 'G'
};

##### END OF GLOBALS #####

my $accession;
my $trim_size = 300;
my $outfile;
my $gene_prefix;
my $genome_name;
my $kingdom = "bac";
#my $rename;
my $procs=8;
my $anno_only="NO";
GetOptions ( #"s" => \$start_job,
            "infile=s"   => \$accession,
            "outfile=s" => \$outfile,
	    "gene_prefix=s" => \$gene_prefix,
	    "name_of_genome=s" => \$genome_name,
	    "kingdom=s" => \$kingdom,
	    "trim_size=i" => \$trim_size,
#            "rename=s" => \$rename,
	    "anno=s" => \$anno_only,
	    );  # flag

die "\nUsage: perl process_new_genome.pl
\t-i refseq_accession or fasta_file_with_contigs
\t-o outfile basename (all outfiles will be called basename.some_suffix)
\t-n genome_name
\t-g gene_prefix (e.g. BAC will have gene names BAC0001, BAC0002, etc...)
\t-t trim_size (i.e. discard contigs shorter than this; default 300bp)
\t-k kingdom (bac, arc); for bacteria also use bacn for gram- and bacp for gram+ if you want to run CELLO too
\t-a Find genes only? or Rename GenBank files only? [Yy]es (default = No)\n\n" unless $accession;

#\t-r rename GenBank files only? [Yy]es (default = No)

my $amino_acid_file = "";
if ($accession =~ /NC_/ || $accession =~ /NZ_/) { # an NCBI genome already has a gene annotation
	my $outdir = get_temp_dir();
	print STDERR "results stored in $outdir\n";
	$amino_acid_file = process_ncbi_genome($accession, $outdir);
	if($anno_only =~ /[Yy]/){die "$accession files parsed\n";}

}
else { # create a ptt file, an rnt file and an amino acid file on an unannotated genome
	die "Usage: perl process_new_genome.pl -i <refseq_accession or fasta_file_with_contigs> -o <basename> -g <gene_prefix> -n [genome_name]\n" unless $accession && $outfile && $gene_prefix;
	$genome_name = "genome" unless $genome_name;
	my $outdir = get_temp_dir();
#	my $outdir = 'ann_5143';
#	print STDERR "results stored in $outdir\n";
	$amino_acid_file = process_new_genome($accession, $outdir, $trim_size, $outfile, $gene_prefix, $genome_name, $kingdom);
#	$amino_acid_file = 'ann_5143/NC_TFaecalibacterium_TS28s.faa';
}

#print "\$amino_acid_file\t[$amino_acid_file]\n";

#$amino_acid_file = 'ann_16787/NC_Bbreve.faa';


run_functional_annotation($amino_acid_file, $kingdom);


################### END OF MAIN ###################


# Function: process_new_genome
# organizing function for contig contenation, gene calling, RNA calling, and functional annotation of a new genome.
#
# Raw contigs are first size filtered and concatenated with a linker.  Gene calls (via Glimmer and GeneMark), 
# rRNAs (via rnammer), and tRNAs (via tRNA-scan) are run on the concatenated genome.
# 
# Parameters:
#   $contig_file - file of contigs for a single organism in FASTA format
#   $outdir - the directory where all of the results will be stored
#   $trim_size - eliminate contigs shorter than this amount
sub process_new_genome {
	my ($contig_file, $outdir, $trim_size, $outfile, $prefix, $genome_name, $kingdom) = @_;

	# concatenate the contigs  ***** create .fna file *****
	my ($single_genome, $single_genome_seq, $contigs) = concatenate_contigs($contig_file, $outdir, $trim_size, $outfile);
#	my @keys = keys %$contigs;
#	die "@keys\n";

	# make a copy of the concatenated genome for glimmer	
	system("cp $single_genome $contig_file.concat");

	# call tRNA and rRNA first, they have priority (i.e. don't allow genes to overlap too much with them)
	my $tRNAs = call_tRNAs($single_genome);
	my $rRNAs = call_rRNAs($contig_file, $contigs, $kingdom);
	my $misc_RNAs = call_misc_RNAs($contig_file,$contigs);##udpated June 2014
	my $tt = call_tt($contig_file,$contigs);##udpated June 2014
#	my $masked_sequence = $single_genome_seq;

	# mask out the tRNAs and rRNAs so we don't call genes there
	#my $edge_to_keep = 10;
	#mask_genome(\$masked_sequence, $rRNAs, $edge_to_keep);
	#$edge_to_keep = 5; # tRNAs are short so be even more conservative;  in RNA-seq theres usually nothin but rRNA expressed around rRNAs
	#mask_genome(\$masked_sequence, $tRNAs, $edge_to_keep);

	#my $masked_genome = "$outdir/$outfile" . "_masked.fna";
	#write_fasta_file($outfile , $masked_sequence, $masked_genome);

	#print STDERR "genome rRNAs and tRNAs masked prior to gene calling (see $masked_genome)\n";

	# now call the coding stuff 
	my $genes = call_genes($contig_file, $contigs);

	my $overlap_window = 10;
	$genes = eliminate_genes($genes, $tRNAs, $overlap_window);
	$genes = eliminate_genes($genes, $rRNAs, $overlap_window);

	add_gene_symbol($genes, 'CDS', $prefix);
	add_gene_symbol($tRNAs, 'tRNA', $prefix);
	add_gene_symbol($rRNAs, 'rRNA', $prefix);
	add_gene_symbol($misc_RNAs, 'misc_RNA', $prefix);##udpated June 2014
	add_gene_symbol($tt, 'TransTerm', $prefix);##udpated June 2014
	my $RNAs = [@$tRNAs, @$rRNAs, @$misc_RNAs, @$tt];
	#die "have " . scalar @$RNAs . " RNAs\n";

	my $all_features = [@$genes, @$tRNAs, @$rRNAs, @$misc_RNAs, @$tt];
	create_nucleotide_features_file($single_genome_seq, $all_features, $outdir .'/'. $outfile);

	# create amino acid file  ***** create *.faa file *****
	my $aafile = create_amino_acid_file($single_genome_seq, $genes, $outdir .'/'. $outfile);
	# create .ptt file
	create_ptt_file($genes, $outdir .'/'. $outfile . '.ptt', $genome_name);
	# create .rnt file
	create_rnt_file($RNAs, $outdir .'/'. $outfile . '.rnt', $genome_name);
	if($anno_only =~/[Yy]/){die "Annotation complete.\nFull stop.\nNo homology search\n\n";}
	return $aafile;
}
sub create_rnt_file {
	my ($RNAs, $outfile, $genome_name) = @_;

	my @RNAs = sort  {$a->{start} <=> $b->{start} } @$RNAs;

	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
	print OUT "$genome_name,\n";
	print OUT scalar @RNAs . " RNAs\n";
	print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";

	for my $r (@RNAs) {
		print OUT "$r->{start}..$r->{stop}\t$r->{dir}\t0\t0\t-\t$r->{symbol}\t-\t-\t$r->{description}\n";
	}
	close OUT;
}

sub create_ptt_file {
	my ($genes, $outfile, $genome_name) = @_;

	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
	print OUT "$genome_name,\n";
	print OUT scalar @$genes . " proteins\n";
	print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";

	for my $g (@$genes) {
		print OUT "$g->{start}..$g->{stop}\t$g->{dir}\t0\t0\t-\t$g->{symbol}\t-\t-\t$g->{type}\n";
	}
	close OUT;
}

sub add_gene_symbol {
	my ($features, $type, $prefix) = @_;

	my $size = scalar @$features;
	my $len = length($size);

	for my $i (0 .. $#$features) {
		my $f = $features->[$i];
		my $feature_id = $prefix . get_id($i+1, $len);  # underscore removed for internal consistency, add back later

#		my $feature_id = $prefix . '_' . get_id($i+1, $len);  # underscore is needed for official locus tag (http://www.ebi.ac.uk/embl/Documentation/locus_tag_usage.html)
#		die "have size $size len $len feature id $feature_id\n";
		if ($type eq 'rRNA') {
			$feature_id .= "r";
		}
		elsif ($type eq 'tRNA') {
			$feature_id .= "t";
		}
		elsif ($type eq 'misc_RNA') {
			$feature_id .= "m";
		}
		elsif ($type eq 'TransTerm') {
			$feature_id .= "tt";
		}
		$f->{symbol} = $feature_id;
	}

}

sub get_id {
	my ($num, $len) = @_;

	while (length($num) < $len) {
		$num = '0' . $num;  # put zeros on the front of the number until all numbers are the same length
	}

	return $num;
}

sub eliminate_genes {
	my ($genes, $reserved, $overlap_window) = @_;

	my @passed_genes;
	for my $g (@$genes) {
		if (no_overlap($g, $reserved, $overlap_window)) {
			push @passed_genes, $g;
		}
	}


	return \@passed_genes;
}

sub no_overlap {
	my ($g, $to_check, $overlap_window) = @_;


	my $start = $g->{start};
	my $stop = $g->{stop};


	for my $t (@$to_check) {
		if ($t->{start} < $start && $t->{stop}-$overlap_window > $start || $t->{start}+$overlap_window > $stop && $t->{stop} < $stop || $t->{start}+$overlap_window < $start && $t->{stop}-$overlap_window > $stop) {
			print STDERR "removing gene $g->{start} $g->{stop} that overlaps with $t->{type} $t->{start} $t->{stop}\n";
			return 0;
		}

	}

	return 1;
}


sub create_nucleotide_features_file {
	my ($seq, $features, $outfile) = @_;
	$outfile .= "_features.fna";

	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
	for my $g (@$features) {
		my $subseq = substr($seq, $g->{start}-1, ($g->{stop}-$g->{start})+1);

		if ($g->{dir} eq '-') {
			$subseq = reverse_complement($subseq);
		}
		print OUT toFasta($g->{symbol} . "|$g->{type}|$g->{description}", $subseq);
	}
	close OUT;
#	die "created features file\n";
}

sub create_amino_acid_file {
	my ($seq, $genes, $outfile) = @_;
	
	$outfile .= ".faa";

	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";

	for my $g (@$genes) {
		my $subseq = substr($seq, $g->{start}-1, ($g->{stop}-$g->{start})+1);

		if ($g->{dir} eq '-') {
			$subseq = reverse_complement($subseq);
		}
#		print toFasta("$g->{description} dir  $g->{dir} start $g->{start} stop $g->{stop}", $subseq);
#		print "\n";
		print OUT toFasta($g->{symbol}, translate_coding_sequence($subseq));
	}
	close OUT;

	return $outfile;
}

sub reverse_complement {
	my $dna = shift;
	my $revcomp = reverse($dna);

	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

sub translate_coding_sequence {
	my $cds = shift;

	my $aa_seq = "";
	while (my $codon = substr($cds, 0, 3, "")) {
#		print "$codon ";
		unless ($TRANS_TABLE->{$codon}) {
			print STDERR "error couldn't find codon $codon\n";
		}
		$aa_seq .= $TRANS_TABLE->{$codon};
	}
	$aa_seq =~ s/^./M/; # replace first character with methionine
	$aa_seq =~ s/\*$//; # replace stop codon

	return $aa_seq;
}

# we replace the $features (minus the amount in $to_keep from each side)
# with a series of NNNN's and if possible a $linker to prevent gene calls
# on top of features we don't want genes called for
sub mask_genome {
	my ($seq_ref, $features, $to_keep) = @_;
	
	$to_keep = 0 unless $to_keep;

	my $linker = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
#	my $linker = "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";
#	my $linker = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
	my $linkerLen = length($linker);
	for my $f (@$features) {
		my $len = ($f->{stop} - $f->{start}) + 1;
		my $start = $f->{start} + $to_keep;
		my $adj_len = $len - ($to_keep*2);

		if ($adj_len > 0) { # if we still have sequence left after the trim
			my $NNNNN = get_N($adj_len);

			# if possible put the stop codons in the mask
			if ($adj_len >= $linkerLen) {
				substr($NNNNN, 0, $linkerLen, $linker);
			}

			my $subseqA = substr($$seq_ref, $f->{start}-1, $len);
			substr($$seq_ref, $start-1, $adj_len, $NNNNN); # replace sequence with N's
			my $subseqB = substr($$seq_ref, $f->{start}-1, $len);

#			print "$f->{description}\n$subseqA\n$subseqB\n\n";
		}
	}
}

sub get_N {
	my $num_N = shift;

	my $n = "";
	for (1 .. $num_N) {
		$n .= "N";
	}

	return $n;
}

sub call_rRNAs {
	my $genome_file = shift;
	my $contigs = shift;
	my $kingdom = shift;

	my $prefix = "c" . int(rand()*1000);
	my $outfile = $prefix . "_rRNA.out";


#	if (0) {
	if ($kingdom =~ /^a/) {	#archea
		print STDERR "Running rnammer using archeal model\n";
		system("$RNAMMER -S arc -m lsu,ssu,tsu -gff $outfile < $genome_file");
	}
	else { # default to bacteria
		print STDERR "Running rnammer using bacterial model\n";
		system("$RNAMMER -S bac -m lsu,ssu,tsu -gff $outfile < $genome_file");
	}
#	}
	my @lines = `cat $outfile`;

#	print "rRNAMER\n@lines\n";

	return trim_rRNAs($contigs, @lines);
}

# from the rnammer paper
# As rRNAs within a genome tend to be very similar, usually with at least 99%
# identity, different full model hits within a genome corresponding to actual 
# rRNAs should be expected to have similar scores. However, we found a substantial
# number of hits with far lower scores which we assume to be pseudogenes, truncated 
# rRNAs or otherwise nonfunctional rRNA copies. To ensure that these did not have
# an adverse effect on the analyses, we excluded full model hits having a score 
# less than 80% of the maximal score in that genome. These are listed in 
# Supplementary Table S2.
sub trim_rRNAs {
	my ($contigs, @lines) = @_;

	my ($best16S, $best23S, $best5S) = (0,0,0);
	my @rRNAs;
	for (@lines) {
		next if /^\s*$/;
		next if /^#/;
		my ($header,$b,$type, $start, $stop, $score, $dir, $d, $description) = split /\t/;
		$header = ">$header";

#		print "have |$header|\n";
		my @keys = keys %$contigs;
#		print "contig list is @keys\n";

		# adjust for concatenated contig
		unless ($contigs->{$header}) {
			next;
			print "warning could not find header $header\n";
		}
		$start = $contigs->{$header} + $start;
		$stop = $contigs->{$header} + $stop;
#		print "still have $header\n";

		# use same naming schema as the NCBI *.rna file
		if ($description =~ /23/) {
			$description = '23S ribosomal RNA';
			$best23S = $score if $score > $best23S; # keep the best one
		}
		elsif ($description =~ /16/) {
			$description = '16S ribosomal RNA';
			$best16S = $score if $score > $best16S; # keep the best one
		}
		elsif ($description =~ /5/) {
			$description = '5S ribosomal RNA';
			$best5S = $score if $score > $best5S; # keep the best one
		}

		# store in hash for easy handling
		my %hash = (
			start => $start,
			stop => $stop,
			dir => $dir,
			description => $description,
			type => $type,
			score => $score,
		);
		push @rRNAs, \%hash;		
#		print "keeping $start\t$stop\t$dir\t$description\t$type\n";

	}
	
	print "best are $best5S\t$best16S\t$best23S\n";
	# trim out the results keeping only those within 80% of max score as suggested by rnammer	
	my @trimmed_rRNAs;
	my $threshold = 0.8; # 80%
	for my $r (@rRNAs) {
		if ($r->{description} =~ /23/) {
			push @trimmed_rRNAs, $r if $r->{score} >= $threshold*$best23S;
		}
		elsif ($r->{description} =~ /16/) {
			push @trimmed_rRNAs, $r if $r->{score} >= $threshold*$best16S;
		}
		elsif ($r->{description} =~ /5/) {
			push @trimmed_rRNAs, $r if $r->{score} >= $threshold*$best5S;
		}
	}

	print STDERR "kept " . scalar @trimmed_rRNAs . " rRNAs after 80% trim\n";

	return \@trimmed_rRNAs;
}

sub call_tRNAs {
	my $genome_file = shift;

	my $prefix = "c" . int(rand()*1000);

	my $outfile = $prefix . "_tRNA.out";
	print STDERR "Running tRNA-Scan\n";
	system("$TRNA_SCAN -q -B $genome_file > $outfile");
	my @lines = `cat $outfile`;

	my $type = 'tRNA';
	my @tRNAs;
	shift @lines; shift @lines; shift @lines; # skip the header

	# should I throw away Undet and Pseudo tRNAs?
	for (@lines) {
		my ($name, $n, $start, $stop, $amino_acid, @stuff) = split /\s+/;
		my $dir = '+';
		my $description = $amino_acid . ' tRNA';		

		if ($start > $stop) {
			my $temp = $start;
			$start = $stop;
			$stop = $temp;
			$dir = '-';   # tRNA-scan gives direction by start..stop orientation
		}
		my %hash = (
			start => $start,
			stop => $stop,
			dir => $dir,
			description => $description,
			type => $type,
		);
#		print "have $start\t$stop\t$dir\t$description\t$type\n";
		push @tRNAs, \%hash;		
	}

	print STDERR "kept " . scalar @tRNAs . " tRNAs\n";

	return \@tRNAs;
}
sub call_misc_RNAs {
	my $genome_file = shift;
	my $contigs = shift;
	
	print STDERR "Running Infernal\n";
	my $cmsearch_tbl='cmsearch.txt';
	system("$INFERNAL -r /data/DB/Core_RFAM_list.txt -g $genome_file -o $cmsearch_tbl");
	my @lines = `cat $cmsearch_tbl`;
	
	my @infernal_res;

#	my $offset = 0;  # offset due to contigs
	my $header;
	foreach my $l (@lines) {
		
		unless ($l=~/^\#/){
			my @cols = split(/\s+/,$l);
			$header=">" . $cols[0];
			#print "[$header][$contigs->{$header}]\n";
			unless ($contigs->{$header}) {
				next;
				print "warning could not find header $header\n";
			}
			my $dir = $cols[9];
			my $description = "$cols[3] $cols[2]";	
			my $type = 'misc_RNA';
			my $start=$cols[7];
			my $stop=$cols[8];
			#print "[$start][$stop]\n";
			$start = $start + $contigs->{$header};
			$stop = $stop + $contigs->{$header};
			if($start > $stop){
				my $temp=$start;
				$start=$stop;
				$stop=$temp;
			}				
			#print "[$start][$stop]\n";
			my %hash = (
				start => $start,
				stop => $stop,
				dir => $dir,
				description => $description,
				type => $type,
			);
#			print "have $start\t$stop\t$dir\t$description\t$type for $header\n";
			push @infernal_res, \%hash;	
		}	
	}

	@infernal_res = sort { $a->{start} <=> $b->{start} } @infernal_res;
	`rm $cmsearch_tbl table.txt`;
	print STDERR "kept " . scalar @infernal_res . " misc_RNAs\n";
	return \@infernal_res;
}
sub call_tt {
	my $genome_file = shift;
	my $contigs = shift;
	
	print STDERR "Running RNIE\n";
	my $cmsearch_tbl=$genome_file;
	$cmsearch_tbl=~s/\.(.+)/\-genomeMode\-rnie\.gff/;
	system("$RNIE -g --genome -md /home/pdegnan/Software/RNIE-master/models/ -f $genome_file");
	my @lines = `cat $cmsearch_tbl`;
	
	my @infernal_res;

#	my $offset = 0;  # offset due to contigs
	my $header;
	foreach my $l (@lines) {
		
		unless ($l=~/^\#/){
			my @cols = split(/\s+/,$l);
			$header=">" . $cols[0];
			#print "[$header][$contigs->{$header}]\n";
			unless ($contigs->{$header}) {
				next;
				print "warning could not find header $header\n";
			}
			my $dir = $cols[6];
			my $description = "$cols[1] $cols[2]";	
			my $type = 'TransTerm';
			my $start=$cols[3];
			my $stop=$cols[4];
			#print "[$start][$stop]\n";
			$start = $start + $contigs->{$header};
			$stop = $stop + $contigs->{$header};
			#print "[$start][$stop]\n";
			my %hash = (
				start => $start,
				stop => $stop,
				dir => $dir,
				description => $description,
				type => $type,
			);
#			print "have $start\t$stop\t$dir\t$description\t$type for $header\n";
			push @infernal_res, \%hash;	
		}	
	}

	@infernal_res = sort { $a->{start} <=> $b->{start} } @infernal_res;
	`rm *-genomeMode*`;
	print STDERR "kept " . scalar @infernal_res . " Rho-independent Transcriptional Terminators\n";
	return \@infernal_res;
}	
sub call_genes {
	my $genome_file = shift;
	my $contigs = shift;
	my @global_features;

	#my $glimmer_out = 'run.predict';
	#print STDERR "Running Glimmer\n";
	#system("$GLIMMER $genome_file run");
	# read glimmer
	#my @lines = `cat $glimmer_out`;
	#my $genes = read_glimmer($contigs, @lines);
	
	print STDERR "Running Prodigal\n";
	my $prodigal_gff='prodigal.gff';
	#my $prodigal_faa='prodigal.faa';
	#my $prodigal_ffn='prodigal.ffn';
	#-a $prodigal_faa -d $prodigal_ffn 
	system("$PRODIGAL -f gff -c -m -i $genome_file -o $prodigal_gff");
	#-m	mask regions containing Ns
	#-f	output format gff
	#-a	protein fasta file
	#-d	nucleotide fasta file
	#-i 	input file
	# read prodigal
	my @lines = `cat $prodigal_gff`;
	my $genes = read_prodigal($contigs, @lines);
	
#	$GENE_MARK_HMMP
#	print STDERR "Running GeneMark\n";
#	my $geneMark_out = $genome_file . '.gm';
	#	system("$GENE_MARK $genome_file --prok > GMSN_STDOUT");
	#     -n      Turn OFF partial gene prediction at long substrings of 'N'; default ON
#	system("$GENE_MARK_HMMP -m GeneMark_hmm_combined.mod -n -r $genome_file -o $geneMark_out > GMHMMP_STDOUT");
# 	read gene mark (ask Brian about the weird message at the end?)
#	my @lines = `cat $geneMark_out`;
#	my $gene_mark = read_genemark(@lines);

	# keep all glimmer that don't overlap more than N bp with geneMark
#	my $genes = combine_gene_calls($gene_mark, $glimmer);

	print STDERR "have " . scalar @$genes ." genes\n";

	# cleanup all of the garbadge these programs leave behind
#	system("rm GeneMark*");
#	system("rm GM*STDOUT");
#	system("rm gms.log");
#	system("rm tmp.*");

	return $genes;
}

sub combine_gene_calls {
	my ($gene_mark, $glimmer) = @_;

	print STDERR "combining gene calls\n";

	my @all_genes = @$gene_mark; # keep everything from genemark

	my $allow_overlap = 10;  # allow 10 bp overlap with other genes
	for my $g (@$glimmer) {
		if (is_new_gene($gene_mark, $g, $allow_overlap)) {
			push @all_genes, $g;
		}
	}

	@all_genes =  sort { $a->{start} <=> $b->{start} } @all_genes;

	return \@all_genes;
}

# slow way to do it, but who cares;  assumes $genes is sorted
sub is_new_gene {
	my ($genes, $gene, $overlap_allow) = @_;

	$overlap_allow=0 unless $overlap_allow;

	my $start = $gene->{start};
	my $stop = $gene->{stop};


	my $g1 = { start => $overlap_allow, stop => $overlap_allow };  # start at position zero
	for my $i (0 .. $#$genes) {
		my $g2 = $genes->[$i];
#		print "have $i start $start stop $stop\t\t$g1->{start} $g1->{stop} $g2->{start} $g2->{stop}\n";

#		print "adding gene $start\t$stop\tfrom glimmer\n";
		if ($start > $g1->{stop}-$overlap_allow && $stop < $g2->{start}+$overlap_allow) {
#			print "adding gene $start\t$stop\tfrom glimmer\n";
			return 1;
		}

		last if ($g1->{start} > $start);  # minor speed up

		$g1 = $g2;
	}

	return 0;
}

sub read_genemark {
	my @lines = @_;


	# get rid of the header stuff
	my $i = 0;
	for (@lines) { $i++; last if /#\s+Length/; }
	#die "have i $i: $lines[$i]\n";

	my @genemark_res;
	for my $j ($i .. $#lines) {
		$_ = $lines[$j];
		s/^\s+//;
		s/\s+$//;
		my ($n, $dir, $start, $stop, @stuff) = split /\s+/;
		my $type = 'CDS';
		my $description = 'genemark gene call';		
		my %hash = (
			start => $start,
			stop => $stop,
			dir => $dir,
			description => $description,
			type => $type,
		);
#		print "have $start\t$stop\t$dir\t$description\t$type\n";
		push @genemark_res, \%hash;		
	}

	@genemark_res =  sort { $a->{start} <=> $b->{start} } @genemark_res;

	return \@genemark_res;
}
sub read_prodigal {
	my ($contigs, @lines) = @_;

	my @prodigal_res;

#	my $offset = 0;  # offset due to contigs
	my $header;
	foreach my $l (@lines) {
		
		unless ($l=~/^\#/){
			my @cols = split(/\s+/,$l);
			$header=">" . $cols[0];
			#print "[$header][$contigs->{$header}]\n";
			unless ($contigs->{$header}) {
				next;
				print "warning could not find header $header\n";
			}
			my $dir = $cols[6];
			my $description = $cols[1];	
			my $type = 'CDS';
			my $start=$cols[3];
			my $stop=$cols[4];
			#print "[$start][$stop]\n";
			$start = $start + $contigs->{$header};
			$stop = $stop + $contigs->{$header};
			#print "[$start][$stop]\n";
			my %hash = (
				start => $start,
				stop => $stop,
				dir => $dir,
				description => $description,
				type => $type,
			);
#			print "have $start\t$stop\t$dir\t$description\t$type for $header\n";
			push @prodigal_res, \%hash;	
		}	
	}

	@prodigal_res = sort { $a->{start} <=> $b->{start} } @prodigal_res;
	
	return \@prodigal_res;
}



sub read_glimmer {
	my ($contigs, @lines) = @_;

	my @glimmer_res;

#	my $offset = 0;  # offset due to contigs
	my $header;
	for (@lines) {
		s/^\s+//;
		s/\s+$//;
		if (/>/) {
			$header = $_;
			print "[$header][$contigs->{$header}]\n";
			#print "have $header and $contigs->{$header}\n";
			next;
		}
		unless ($contigs->{$header}) {
			next;
			print "warning could not find header $header\n";
		}
#		else {
#			print "found header $header\n";		
#		}
		
		
		my ($n, $start, $stop, @stuff) = split /\s+/;
		my $dir = '+';
		my $description = 'glimmer gene call';		
		my $type = 'CDS';

		if ($start > $stop) {
			my $temp = $start;
			$start = $stop;
			$stop = $temp;
			$dir = '-';   # tRNA-scan gives direction by start..stop orientation
		}
		print "[$start][$stop]\n";
		$start = $start + $contigs->{$header};
		$stop = $stop + $contigs->{$header};
		print "[$start][$stop]\n";
		my %hash = (
			start => $start,
			stop => $stop,
			dir => $dir,
			description => $description,
			type => $type,
		);
#		print "have $start\t$stop\t$dir\t$description\t$type for $header\n";
		push @glimmer_res, \%hash;		
	}

	@glimmer_res = sort { $a->{start} <=> $b->{start} } @glimmer_res;

	return \@glimmer_res;
}

# concatenate all of the sequences that pass the length threshold, and write a new fasta file
sub concatenate_contigs {
	my ($contig_file, $outdir, $trim_size, $filebase) = @_;
	my $concatamer_file = $outdir . "/$filebase" . ".fna";
	my $contig_outfile = $outdir . "/$filebase" . ".contig";

	# Assign the linker 
	# this is from the Institute for Genome Science at the University of Maryland annotation pipeline, 
	# and should introduce a stop codon in all 6 reading frames.
#	my $linker = "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";
#	my $linker = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
	# the above linker doesn't prevent gene calling
	my $linker = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	my $linkerLength = length($linker);

	my $seqs = read_seq_to_array($contig_file);

	my $concatamer_seq = "";
	my %contigs;
	my @positions; # store contig start, stop and header
	my $count=0;
	for my $s (@$seqs) {
		my ($header, $seq) = @$s;
			
		if (length($seq) >= $trim_size) {
			print STDERR "adding $header\n";
			if ($count > 0) {   
				$concatamer_seq .= $linker;	
			}
			my $contig_start = length($concatamer_seq);
			my $contig_stop  = $contig_start + length($seq);
			push @positions, [$contig_start, $contig_stop, $header];
			$count++;
			$contigs{$header} = $contig_start;
			$concatamer_seq .= $seq;
		}	
		else {
			print STDERR "skipping $header (too short)\n";
		}
	}

	if (scalar @$seqs > 1) {   # if there are more than one contig; assume it is NOT a circle 
				   # so add a linker at the end to keep annotation callers from calling "around" the circle
		$concatamer_seq .= $linker;	
	}
	$concatamer_seq = uc($concatamer_seq);	

	for my $h (keys %contigs) {
		# make a cleaned up header that only contains the first text;  RNAMMER trims the fasta header in its output so this is important
		my $new = $h;
		$new =~ s/^(>\w+)\s+.*$/$1/;
#		print "converting $h to $new\n";
		$contigs{$new} = $contigs{$h};
	}
#	$contigs{$header} = $contig_start;

#	die;
	write_fasta_file($filebase, $concatamer_seq, $concatamer_file);
	write_contig_file(\@positions, $contig_outfile);
	print "total concatenated genome length is " . length($concatamer_seq) . "\n";

	return $concatamer_file, $concatamer_seq, \%contigs;
}

sub write_contig_file {
	my ($positions, $outfile) = @_;

	open (OUT, ">$outfile") or die "Can't open $outfile: $!\n";
	for my $p (@$positions) {
		my $text = join "\t", @$p;
		print OUT "$text\n";
	}
	close OUT;
}

sub write_fasta_file {
	my ($header, $seq, $outfile) = @_;

	open (OUT, ">$outfile") or die "Can't open $outfile: $!\n";
	print OUT toFasta($header, $seq);
	close OUT;
}

sub toFasta {
        my ($seqName,$seq,$len) = @_;

        $len = 80 unless $len;

        my $formatted_seq = ">$seqName\n";
        while (my $chunk = substr($seq, 0, $len, "")) {
                $formatted_seq .= "$chunk\n";
        }

        return $formatted_seq;
}

# would explode for large genomes
# puts in a 2d array with $array[$i][0] = header and $array[$i][1] = sequence
sub read_seq_to_array {
	my $seqfile = shift;

	my @seqs;
	my $seq;
	my $header;
	open (IN, $seqfile) or die "can't open $seqfile: $!\n";
	while (<IN>) {
		chomp;
		next if /^\s$/;
		if (/^>/) {
			if ($seq) {
				push @seqs, [$header, $seq];			
			}		
			$seq = "";
			$header = $_;
		}
		else {
			s/\s+//;  # clean up any whitespace
			$seq .= $_;
		}
	}
	close IN;

	if ($seq) { # add the last one
		push @seqs, [$header, $seq];			
	}

	return  \@seqs;
}

# Function: process_ncbi_genome
# rename protein headers in *.faa and return the filename where new amino acids have been written
# 
# This function is called when the file name passed on the commandline contains NM_*
# 
# Parameters:
# $accession - the RefSeq accession (used to locate the *.ptt and *.faa files)
# $outdir - directory to write the modified protein FASTA file
sub process_ncbi_genome {
	my $accession = shift;
	my $outdir    = shift;

	my $ptt = $accession . ".ptt";
	my $faa = $accession . ".faa";

	# first task is to make an amino acid file containing all genes
	# using a shorter name that we can parse better (for now uses gene symbol)
	my $pid_to_geneSymbol  = parse_ptt($ptt);
	my $outfile = file_from_fullpath($accession);
	my $aafile = $outdir . "/$outfile";
	print_aa_with_symbol($pid_to_geneSymbol, $faa, $aafile);

	return $aafile; # FASTA file of amino acid sequences where each sequence header uses the gene symbol only
}

# Function: run_functional_annotation
# call all of the functions that start functional annotation jobs 
#
# Currently hmmpfam is using the Pfam and TIGRFAM databases (function <run_HMMPFAM>).  While blast is run against
# the STRING database with COG annotations (function <run_COG>).
#
# Parameters:
# $protein_file - FASTA protein file to run the functional annotation jobs against
sub run_functional_annotation {
	my $protein_file=shift;
	my $kingdom = shift;

	# run a COG blast annotation
	run_COG($protein_file);

	# run a KEGG blast annotation
	run_KEGG($protein_file);

	# run against TIGRFAM and PFAM
	run_HMMPFAM($protein_file);

	# run CELLO to predict subcellular localization
	run_CELLO($protein_file, $kingdom) if $kingdom =~ /bacn/ || $kingdom =~ /bacp/;
}

sub run_CELLO {
	my $protein_file = shift;
	my $kingdom = shift;


	my $num = int(rand() * 10000);
	my $temp_dir = "CELLO$num";
	system("mkdir $temp_dir");

	# determine if gram positive or gram negative
	my $gram_stain = "";
	if ($kingdom =~ /bacn/) {
		$gram_stain = "gramn";
	}
	elsif ($kingdom =~ /bacp/)  {
		$gram_stain = "gramp";

	}

	my $outfile = $protein_file;
	$outfile =~ s/\.faa$// if $outfile =~ /\.faa/;
	$outfile .= ".CELLO";

	my $command = "$CELLO $protein_file $gram_stain prot $temp_dir > $outfile";
	#my $header = '#!/bin/sh';


	my $batch_file = "CELLO$num.sh";
	system("echo '$command' > $batch_file");

	print STDERR "starting CELLO job: $command from $batch_file\n";
	
	my $qsub ="sh $batch_file &"; 
	system("$qsub");
}

# returns just the file name, trimming off all of the directory stuff
sub file_from_fullpath {
	my $path = shift;
	
	my (@pieces) = split /\//, $path;

	return $pieces[$#pieces];
}

# create a directory to store things in temporarily
sub get_temp_dir {
	my $id = "ann_" . int(rand() * 100000);

	system("mkdir $id");
	return $id;
}

# Function: run_HMMPFAM
# run hmmpfam against the TIGRFAM and PFAM databases
# 
# Uses split_run_check_combine.pl to run everything in batch and combine
#
# Parameters:
# $hmmpfam_input - FASTA protein file to run hmmpfam against
sub run_HMMPFAM {
	my $hmmpfam_input = shift;
	#my $query = "'$BLAST -p blastp -d $STRING -i INCLUDE_INFILE -m 8 -e 10e-10'";
        #my $program_tigrfam    = "$HMMPFAM --cut_nc --acc --cpu 0 $TIGRFAM INCLUDE_INFILE";
        my $program_tigrfam    = "$HMMSCAN --cut_nc --acc $TIGRFAM14 INCLUDE_INFILE";
        
	#my $program_pfam = "$HMMPFAM --cut_nc --acc --cpu 0 $PFAM INCLUDE_INFILE";
	my $program_pfam = "$HMMSCAN --cut_nc --acc $PFAM27 INCLUDE_INFILE";
	my $num_per_job = $procs;
	
	# run tigrfam
	my $command = "$BATCH -i $hmmpfam_input -p '$program_tigrfam' -P $num_per_job -s TIGRFAM";
	print STDERR "running TIGRFAM: $command on the cluster\n";
	# run it
	system($command);
	
	# run pfam 
	$num_per_job = $procs;  # pfam database is much bigger and requires splitting jobs up further
	$command = "$BATCH -i $hmmpfam_input -p '$program_pfam' -P $num_per_job -s PFAM";
	print STDERR "running PFAM: $command on the cluster\n";
	# run it
	system($command);
}

# Function: run_COG
# run blast against the STRING/COG database
# 
# 
# Uses split_run_check_combine.pl to run everything in batch and combine
#
# Parameters:
# $blast_input - FASTA protein file to run hmmpfam against
sub run_COG {
	my ($blast_input) = @_;
	my $num_per_job = $procs;   # 40 results takes about 15 minutes against the COG database

	# set up the blast query
	#my $query = "'$BLAST -p blastp -d $STRING -i INCLUDE_INFILE -m 8 -e 10e-10'";
	#my $query = "'$BLASTP  -db $STRING -query INCLUDE_INFILE -evalue 10e-10 -outfmt 6 -max_target_seqs 20'";
	my $query = "'$PHMMER --tblout INCLUDE_OUTFILE -E 10e-10 INCLUDE_INFILE $STRING'";
	
	# put it into the format for the batch job processor
	my $command = "$BATCH -i $blast_input -p $query -P $num_per_job -s stringCOG";
	
	# run it
	system($command);

	# help for debugging
	print STDERR "running COG on the cluster: $command\n";
}


# Function: run_KEGG
# run blast against the KEGG database
# 
# 
# Uses split_run_check_combine.pl to run everything in batch and combine
#
# Parameters:
# $blast_input - FASTA protein file to run hmmpfam against
# notes for KEGG
# blastall -p blastp -m 8 -e 1e-10  -i 
# blastall -p blastp -m 8 -e 1e-10 -i NC_DpigerGOR1.faa -d ~/KO_blastdb/KO_blastdb2.faa > NC_DpigerGOR1.KEGG
sub run_KEGG {
	my ($blast_input) = @_;
	my $num_per_job = $procs;   # 40 results takes about 15 minutes against the KEGG database

	# set up the blast query
	#my $query = "'$BLAST -p blastp -d $KEGG -i INCLUDE_INFILE -m 8 -e 1e-10'";
	#my $query = "'$BLASTP  -db $KEGG -query INCLUDE_INFILE -evalue 10e-10 -outfmt 6 -max_target_seqs 20'";
	my $query = "'$PHMMER --tblout INCLUDE_OUTFILE -E 10e-10 INCLUDE_INFILE $KEGG'";
	
	# put it into the format for the batch job processor
	my $command = "$BATCH -i $blast_input -p $query -P $num_per_job -s KEGG";
	
	# run it
	system($command);

	# help for debugging
	print STDERR "running KEGG on the cluster: $command\n";
}

# print the amino acids from the faa file; but change the header to only include the gene symbol (this makes things easier for database entry)
sub print_aa_with_symbol {
	my ($pid_to_symbol, $faa, $outfile) = @_;

	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";

	open (IN, $faa) or die "can't open $faa: $!\n";
	while (<IN>) {
		if (/^>gi\|(\d+)\|/) {
			my $pid=$1;
			print OUT ">$pid_to_symbol->{$pid}\n";
		}
		else {
			print OUT;
		}
	}
	close OUT;
	close IN;

	#die;
}

# parse the ptt file to figure out how to convert from protein-id to gene symbol
sub parse_ptt {
	my $in=shift;
	open (IN, $in) or die "can't open $in: $!\n";
	open (ABV, ">>../abbreviations.txt") or die "can't open abbreviations.txt $!\n";
	my $new_format;
	#$new_format = 1 if $in =~ /^NZ/;### genbank NZ's look fine

	# there are many header lines
	my $header = <IN>;
	$header = <IN>;
	$header = <IN> unless $new_format;
	my %pid_to_symbol;
	my $abbreviation;
	while (<IN>) {
		chomp;
		#next if /^\s*$/;
		next unless /^\d+\.\.\d+/; ##necessary to prevent unnecessary warnings... when ptt files are cat'ed together
		my ($loc, $strand, $len, $pid, $gene, $symbol, @junk, $product, $start, $end, $locus, $blah);

		if ($new_format) {
			my @junk = split /\t/;
			my $text =  join "|", @junk;
	#		print "$text\n";
#			die "here: $text\n";
			
			($product, $start, $end, $strand, $len, $pid, $locus, $blah, $symbol, @junk) = split /\t/;
			
		}
		else {
			($loc, $strand, $len, $pid, $gene, $symbol, @junk) = split /\t/;
		}
		#print "have $pid\t$symbol\t$_\n";
#	die;	
		$symbol =~ /(\S+?)\_/;
		$abbreviation=$1;
		$symbol =~ s/_//g;  # need to remove underscores for consistency with the database
		$pid_to_symbol{$pid}=$symbol;
	}
	print ABV "$in\t$abbreviation\n";
	close IN;
	close ABV;
#die;

	return \%pid_to_symbol;
}

