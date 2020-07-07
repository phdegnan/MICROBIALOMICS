#!/usr/bin/perl -w

# create all of the externally needed files for lwgv.js and microbialomics.js
# these include track files and BLAST database files

use strict;
use lib('standard_modules');
use microbialomics;

my $DB = shift;
my $make_default = shift;

$DB   = "microbialomics_dev" unless $DB;
my $host = "localhost";
my $user = "root";
my $pass = "ptgdb25";

my $microbe = new microbialomics($DB, $host, $user, $pass);
my $BLAST_DB_DIR = "/usr/local/share/blast/";
$BLAST_DB_DIR .= $DB unless $DB eq "microbialomics" || $make_default;


create_blast_DBs();
create_lwgv_files();

sub create_lwgv_files {
	print STDERR "creating lwgv annotation files\n";
	my $dbh=$microbe->{DBH};
	my $sth=$dbh->prepare("SELECT genome_constant_id, genome_name, genome_length, genome_sequence from genome");
	$sth->execute();

	
	while (my ($gid, $gname, $glen, $gseq)  = $sth->fetchrow_array()) {
		$microbe->set_genome($gid);
		print "have $gid\n";

		# create the contig tracks if the genome is not complete
		my $out = "./lwgv/" . $gid . ".contig";
		my $contig = $microbe->genome_to_lwgv_contig();
		if ($contig) {
			open (OUT, ">$out") or die "can't open $out: $!\n";
			print OUT $contig;
			close OUT;
		}

		# create the basic ann file for the species
		create_lwgv_ann($gid, $gname, $glen, $contig);

		# create the sequence file for the species
		$out = "./lwgv/" . $gid . ".seq";
		open (OUT, ">$out") or die "can't open $out: $!\n";
		print OUT $gseq;
		close OUT;

		# create the gene annotation tracks
		$out = "./lwgv/" . $gid . ".genes";
		open (OUT, ">$out") or die "can't open $out: $!\n";
		print OUT $microbe->genome_to_lwgv();
		close OUT;


		# for alejandro
#		my $phage_omics_blast = $microbe->blast_to_lwgv();
#		$out = "./lwgv/" . $gid . ".blast";
#		if ($phage_omics_blast) {
#			open (OUT, ">$out") or die "can't open $out: $!\n";
#			print OUT $phage_omics_blast;
#			close OUT;
#		}
	}
}

sub create_lwgv_ann {
	my ($gid, $gname, $glen, $contig) = @_;
	my $outfile = "./lwgv/$gid.ann";
	open (LWGV_OUT, ">$outfile") or die "can't open $outfile: $!\n";

	# basic configuration stuff
	print LWGV_OUT
	      "#config TMP_PIC_DIR /Library/WebServer/Documents/tmp_img/\n".
	      "#config TMP_PIC_DIR_URL /tmp_img/\n".
	     # "#config IMAGE_WIDTH 30500\n".
	     # "#config COMPRESSION_LIMIT 0\n".
	      "#config GRAPH_HEIGHT 100\n".
	      "#config JSON_DEFAULTS true\n".
	      "#config SEQUENCE_FILE $gid.seq\n";
	    #  "#config SHOW_ZOOM_BAR false\n".
	    #  "#config DONT_QUOTE_LINKS true\n".
	    #  "#config SIDE_NAMES false\n".
	    #  "#config TRACK_SIDE_PADDING 0\n";
	# clean up the genome name into something ok in lwgv language
	$gname =~ s/-//g;
	$gname =~ s/substr.//g;
	$gname =~ s/str.//g;
	$gname =~ s/\s+/_/g;

	print LWGV_OUT "\nbegin genome $gname\n";
	# include contig boundaries if genome not complete
	print LWGV_OUT "\t#include $gid.contig\n" if $contig;
	# include the gene annotation
	print LWGV_OUT "\t#include $gid.genes\n";
	print LWGV_OUT "end genome\n";
	print LWGV_OUT "track contigs setTrackColor(100,100,255)\n" if $contig;
	print LWGV_OUT "setGenomeLength($glen)\n";
	print LWGV_OUT "\nshowGenome($gname)\n";

	close LWGV_OUT;
}

sub create_blast_DBs {
	print STDERR "creating blast databases\n";

	#nucleotides (genome files)
	my $fna_out = "microbialomics.fna";
	write_fasta_nucs($fna_out);
	system("formatdb -i $fna_out -p F"); 

	#amino acids (protein files)
	my $faa_out = "microbialomics.faa";
	write_fasta_aa($faa_out);
	system("formatdb -i $faa_out -p T"); 

	# move the new database to their spot
	system("mv $fna_out* $BLAST_DB_DIR");
	system("mv $faa_out* $BLAST_DB_DIR");
	# remove unneccessary log file
	system("rm formatdb.log");
}

# write nucleotide sequences for making a blast db
sub write_fasta_nucs {
	my $out=shift;

	my $dbh=$microbe->{DBH};
	my $sth=$dbh->prepare("SELECT genome_constant_id, genome_sequence from genome");
	$sth->execute();

	open (OUT, ">$out") or die "can't open $out: $!\n";

	while (my @res = $sth->fetchrow_array()) {
		print OUT $microbe->to_fasta($res[0], $res[1]);
	}
}

# write nucleotide sequences for making a blast db
sub write_fasta_aa {
	my $out=shift;

	my $dbh=$microbe->{DBH};
	my $sth=$dbh->prepare("select genome_constant_id, feature_symbol, protein_sequence FROM genome g, feature f, protein_seq p WHERE p.feature_id=f.feature_id AND f.genome_id=g.genome_id AND primary_feature='Y'");
	$sth->execute();

	open (OUT, ">$out") or die "can't open $out: $!\n";

	while (my @res = $sth->fetchrow_array()) {
		my $header = $res[0] . ":"  . $res[1];
		print OUT $microbe->to_fasta($header, $res[2]);
	}
}
