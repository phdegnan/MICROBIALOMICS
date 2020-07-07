#!/usr/bin/perl -w
##
##	Brian Muegge & J Faith 2008-2009
## 	Patrick Degnan revised 2014
##	load_completed_genome3.pl
##	Upload annotation data to MySQL database
##

use lib('/Users/Shared/Scripts/degnan/build_scripts/');####
use strict;

use microbialomicsDB;
# all the parsers are separate modules to make things modular (i.e. if we abandon a method, we don't have a bunch of dead code hanging around)
use load_pfam3;
use load_COG3;
use load_express;
use load_KO;
use load_CELLO;
use load_KEGG3;
use load_CAZy;
use Getopt::Long;
use GO::Parser;



my $accession;
my $clear_db;
my $dir = "/Users/user/";####

my $GO_obo_file = "/Users/Shared/DB/gene_ontology.obo";###

#my $DB   = "microbialomics";
my $DB   = "XXXXX";
my $host = "localhost";
my $user = "user";
my $pass = "XXXX";
#my $first = "No";
GetOptions ( #"s" => \$start_job,
            "DB=s" => \$DB,
            "clear_db" => \$clear_db,
            "directory=s" => \$dir,
            "accession=s" => \$accession,
#            "first=s" => \$first,
            );  # flag

$dir .= "/" unless $dir =~ /\/$/;

my $microbeDB = new microbialomicsDB($DB, $host, $user, $pass);
$microbeDB->clearOldDatabase() if $clear_db;

die "usage: load_completed_genome.pl -a <accession(s) comma separated> -DB [database] -c [clear database] -dir [annotation file directory] \n" unless $accession;
#-f [First genome Yes? default = No]

# initialize the slow go parser
my $GOparser = new GO::Parser({handler=>'obj'});
# and read the giant ontology file
read_GO($GO_obo_file);

my @accessions = split /,/, $accession;

my $start = time();
my $count = 0;
my $num_genomes = scalar @accessions;
for my $a (@accessions) { 
	load_genome($a, $microbeDB); 

	$count++;
        my $time_so_far = time()-$start;
        my $time_per_round = $time_so_far/$count;
        my $minutes_remaining = ($time_per_round * ($num_genomes-$count))/60;

	printf "finished %d of %d, %d minutes remaining\n", $count, $num_genomes, $minutes_remaining;
}

## END OF MAIN ##

sub load_genome {
	my $accession = shift;

	my $base_file = $dir . $accession;

	my @annotation;

	# GET THE INITIAL GENE/RNA POSITIONS AND SEQUENCES
	# load protein annotation
	my $genome_name = load_annotation($base_file . ".ptt", \@annotation, "CDS");
	# load RNA annotation
	load_annotation($base_file . ".rnt", \@annotation);
	# grab the sequence
	my $genome_seq = load_sequence($base_file . ".fna");
	my $aa_seq = load_aa_sequence($base_file . ".faa");

	# LOAD THEM INTO DB
	my $genomeInfo = load_into_DB($genome_name, $genome_seq, \@annotation, $aa_seq, $accession);

#	return;

	# load the contigs if they genome isn't continuous
	if (file_exists($base_file . ".contig")) {
		load_contigs($base_file . ".contig", $genomeInfo->{genome_id});
	}

	# load M3D expression info
	if (file_exists($base_file . ".express")) {
		print "loading expression\n";
		my $express = new load_express($microbeDB);
		$express->load_express_file($base_file . ".express", $genomeInfo->{genome_id});
	}
	else { print "no expresssion results for $genome_name\n" }
	
	# NOW LOAD WHATEVER FUNCTIONAL ANNOTATIONS ARE AVAILABLE
	#CAZy
	if (file_exists($base_file . ".CAZy")) {
		my $CELLO = new load_CAZy($microbeDB);
		$CELLO->load_CAZy_file($base_file . ".CAZy", $genomeInfo->{genome_id});
	}
	else { print "no CAZy results for $genome_name\n" }
	# NOW LOAD WHATEVER FUNCTIONAL ANNOTATIONS ARE AVAILABLE
	#KEGG
	if (file_exists($base_file . ".KEGG")) {
		print "loading KEGG results\n";
		my $KEGG = new load_KEGG($microbeDB);
		$KEGG->load_KEGG_file($base_file . ".KEGG", $genomeInfo->{genome_id});
	}
	else { print "no KEGG results for $genome_name\n" }

	#CELLO
	if (file_exists($base_file . ".CELLO")) {
		my $CELLO = new load_CELLO($microbeDB);
		$CELLO->load_CELLO_file($base_file . ".CELLO", $genomeInfo->{genome_id});
	}
	else { print "no CELLO results for $genome_name\n" }

	#pfam
	if (file_exists($base_file . ".PFAM")) {
		my $pfam = new load_pfam($microbeDB, $GOparser);
		$pfam->load_pfam_file($base_file . ".PFAM", $genomeInfo->{genome_id});
	}
	else { print "no pfam results for $genome_name\n" }
	
	#tigrfam
	if (file_exists($base_file . ".TIGRFAM")) {
		my $pfam = new load_pfam($microbeDB, $GOparser);
		$pfam->load_pfam_file($base_file . ".TIGRFAM", $genomeInfo->{genome_id});
	}
	else { print "no tigrfam results for $genome_name\n" }


	#stringCOG
	if (file_exists($base_file . ".stringCOG")) {
		my $COG = new load_COG($microbeDB);
		$COG->load_COG_file($base_file . ".stringCOG", $genomeInfo->{genome_id});
	}
	else { print "no stringCOG results for $genome_name\n" }

	# NOW LOAD WHATEVER GENETIC STRAINS ARE AVAILABLE
	# KO libraries from Andy
	if (file_exists($base_file . ".KO")) {
		my $KO = new load_KO($microbeDB);
		$KO->load_KO_file($base_file . ".KO", $genomeInfo->{genome_id});
	}
	else { print "no knockout results for $genome_name\n" }
}

sub read_GO {
	my $obo_file = shift;

	print STDERR "reading go (very slow; but only happens one time)\n";
#	$GOparser->parse($obo_file);
#       die "read go parser\n";
}



sub file_exists {
	my $file=shift;
	return 1 if (-e $file);
	return 0;
}

sub load_into_DB { 
	my ($genome_name, $genome_seq, $annotation, $aa_seq, $accession) = @_;
	my $genome_hash = format_genome_hash($genome_name, $genome_seq, $accession);
	$microbeDB->load_genome_table($genome_hash);
	$annotation = format_annotation($annotation, $genome_hash);
	$microbeDB->load_annotation_table($annotation);
	load_aa_seqs($aa_seq, get_gid_to_feature_id($annotation), get_feature_symbol_to_feature_id($annotation));
	return $genome_hash;
}

sub load_contigs {
	my $in = shift;
	my $genome_id=shift;

	my $first_line=1;
	open (IN, $in) or die "can't open $in: $!\n";
	my @rubbish;
	my $format = 'NCBI';
	if($in !~ /NZ_|NC_/){$format = 'other' }
	while (<IN>) {
		chomp;
		my ($start, $stop, $description);
		#if ($first_line) {
		#	$first_line=0;
		#	if (/whole genome shotgun/) { # contig info downloaded from NCBI
		#		my $track = <IN>;
		#		next;
		#	}
		#	else {
		#		$format = 'other';
		#	}
		#}

		if ($format eq 'NCBI') {
			($description, $start, $stop, @rubbish) = split /\t/;
		}
		else {
			($start, $stop, $description) = split /\t/;
		}

		$description =~ s/^>//;
		my %hash;
		$hash{genome_id} = $genome_id;
		$hash{start_pos} = $start;
		$hash{stop_pos} = $stop;
		$hash{contig_description} = $description;
		$microbeDB->load_contig(\%hash);
	#	print "have $start $stop $description\n";
	}

}

sub load_aa_seqs {
	my ($aa_seq, $gid_to_fid, $symb_to_fid) = @_;

	my @aa_seqs;
	for my $k (keys %$gid_to_fid) {
#		print "here with $k\n";	
		if ($aa_seq->{$k}) { # if we have a sequence for gid store it
#			print "inside with $k\n";	
			my %hash;
			$hash{feature_id} = $gid_to_fid->{$k};
			$hash{protein_sequence} = $aa_seq->{$k};
			push @aa_seqs, \%hash;
		}
	}
	for my $k (keys %$symb_to_fid) {
#		print "sym with $k\n";	
		if ($aa_seq->{$k}) { # if we have a sequence for gid store it
#			print "inside with $k\n";	
			my %hash;
			$hash{feature_id} = $symb_to_fid->{$k};
			$hash{protein_sequence} = $aa_seq->{$k};
			push @aa_seqs, \%hash;
		}
	}
	$microbeDB->load_aa_seqs(\@aa_seqs);

#	printf "have %d aa seqs with match\n", scalar @aa_seqs;
}

sub get_feature_symbol_to_feature_id {
	my $ann=shift;

	my %fsymb_to_fid;
	for my $a (@$ann) { $fsymb_to_fid{$a->{feature_name}} = $a->{feature_id}; }
	return \%fsymb_to_fid;
}

sub get_gid_to_feature_id {
	my $ann=shift;

	my %gid_to_fid;
	for my $a (@$ann) { $gid_to_fid{$a->{genbank_id}} = $a->{feature_id};} #print "have $a->{genbank_id}\t$a->{feature_id}\n"}
	return \%gid_to_fid;
die;
}

# put into the format needed by the DB
sub format_annotation {
	my $annotation=shift;
	my $genome=shift;
	for (@$annotation) {
		$_ = format_annotation_hash( $_, $genome->{genome_id});
	}
	return $annotation;
}

sub format_annotation_hash {
	#my $init=shift;
	#my $feature_id = "";
	#die "$init\t[$feature_id]\n";
	#if($init == 1 && $first =~/[Yy]/){$feature_id = 1; print "$init\t[$feature_id]\n";}
	my $ann=shift;
	my $genome_id=shift;
	my ($start, $stop) = split /\.\./, $ann->{location};
	my $dir = $ann->{strand};
	if ($dir eq "+") { $dir = "F"; }
	elsif ($dir eq "-") { $dir = "R"; }
	my $name = $ann->{gene};
	my $symbol = $ann->{synonym};
	undef $name  if $name eq "-";
	$symbol =~ s/_//g; # clean out the underscores
	$name = $symbol unless $name;
	
	my @keys = keys %$ann;
	my $product = $ann->{product};
	my $type = "";
	$type = $ann->{type} if $ann->{type};
	if ($ann->{product} && !$type) {
		if ($ann->{product} =~ /5S/) { $type = "5S rRNA"; }
		elsif ($ann->{product} =~ /16S/) { $type = "16S rRNA"; }
		elsif ($ann->{product} =~ /23S/) { $type = "23S rRNA"; }
		elsif ($ann->{product} =~ /tRNA/) { $type = "tRNA"; }
		elsif ($ann->{product} =~ /sRNA/) { $type = "sRNA"; }
		elsif ($ann->{synonym} =~ /\d+m/) { $type = "miscRNA"; }
		elsif ($ann->{synonym} =~ /\d+tt/) { $type = "TransTerm"; }
		else { $type = "CDS";} # hope this assumption is ok
	}
	if($product eq "CDS"){$product = "hypothetical protein";}##
#	print "here with $name\t$symbol\t$ann->{product}\t$ann->{pid}\n";
	# "" /NULL/ default fail to work for feature_id and feature_score. Use 0 instead
	my $annHash = {
		feature_id => 0,
		genome_id => $genome_id,
		genbank_id => $ann->{pid},
		feature_name => $name,
		feature_symbol => $symbol,
		feature_description => $product,
		feature_type => $type,
		feature_evidence => "",
		feature_score => 0,
		start_pos => $start,
		stop_pos => $stop,
		direction => $dir,
		primary_feature => "Y",
	};
	return $annHash;
}

sub format_genome_hash {
	my ($genome_name, $genome_seq, $accession) = @_;
	my @pieces = split /\s+/, $genome_name;
	my $genus_species = "$pieces[0] $pieces[1]";
	# make all of the database columns
	# "" /NULL/ default fail to work for genome_id and pmid. Use 0 instead
	my $genomeHash = {
		genome_id => 0,
		genome_constant_id => $accession,
		genome_name => $genome_name,
		genome_description => "",
		genus_species => $genus_species,
		genome_pi => "",
		genome_length => length($genome_seq),
		genome_sequence => $genome_seq,
		sequencing_date => '2014-07-10',
		pmid => 0,
		sequencing_method => "Illumina - NexTera",
		assembly_method => "a5ud_pipeline",
		aerobic_tolerance => "",
		growth_conditions => "",
		public_genome => 1,
		symbiosis => "",
		gram_stain => "",
		isolated_from => "",
		host_phenotype => "",
		host => "",
		host_age_in_years => "",
		host_sex => "",
		host_origin => "",
		literature => "",
		expression_db => "",
	};

	if ($microbeDB->{GENOME_INFO}{$accession}) {
		my $genome_info = $microbeDB->{GENOME_INFO}{$accession};

		for my $k (keys %$genomeHash) {
			if ($genome_info->{$k}) {
				$genomeHash->{$k} = $genome_info->{$k};
				print "assigning $k to $genome_info->{$k} for $accession\n";
			}
		}
	}
#	for my $k
	return $genomeHash;
}

sub load_annotation {
	my $ann_file = shift;
	my $annotation = shift;
	my $type = shift;
	if($ann_file =~/\.rnt/){##
		open (IN, $ann_file) or return(0);##
	}else{###
		open (IN, $ann_file) or die "can't open $ann_file: $!\n";
	}##
	my @pieces = split /\//, $ann_file;
	my $accession = $pieces[$#pieces];
	$accession =~ s/\..*//;
#	die "accession is $accession\t@pieces\n";

	# three header lines
	my $header = <IN>;
	my $genome_name = "";
	if($header =~ /(.*?),/){
		$genome_name = $1; 
	}elsif($header =~ /(.*?)\s+\-\s+1\.\./){
		$genome_name = $1; 
	}else{
		$genome_name = $header; 
	}
	$genome_name=~ s/\s+chromosome//; ### added phd
	$genome_name=~ s/\s+draft genome\.//; ### added phd
	$genome_name=~ s/\s+genome//; ### added phd
	# new genomes have one less line
	#if ($ann_file !~ /NZ_/)	{$header = <IN>; }
	$header = <IN>;
#	print "have $ann_file\n";

	# last header has the column names
	$header = <IN>;
	chomp $header;
	my @headerCol = split /\t/, $header;
	for (@headerCol) { $_=lc($_); } # make case insensitive
	#die "have header @headerCol " .scalar @headerCol. " \n";
	my $count=0;
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		my @cols = split /\t/;
		#Handling descriptions that are longer than 255 characters:
		if(length($cols[8]) > 250){
			#print "$cols[8]\n";
			if($cols[8]=~/(.+?)\;/){
				$cols[8]=$1;
			}else{
				my $sub = substr(0,250,$cols[8]);
				$cols[8]=$sub;
			}
			#print "$cols[8]\n";
		}
		my %hash;
		for my $i (0 .. $#cols) { # put into a hash indexed by the column name
		#	$cols[$i] = "" if $cols[$i] eq '-';
			$hash{$headerCol[$i]} = $cols[$i]; 
		} 
		$hash{type} = $type if $type;
	
		# handle the Refseq FTP format and the table format for newer genomes not in FTP
		if (!$hash{location} && ($hash{start} && $hash{end})) {
			$hash{location} = "$hash{start}..$hash{end}";
		}
		if (!$hash{product} && $hash{'product name'}) {
			$hash{product} = $hash{'product name'};
		}
		if (!$hash{synonym} && $hash{locus_tag}) {
#			print "here\n";
			$hash{synonym} = $hash{locus_tag};
		}
		if (!$hash{pid} && $hash{gi}) {
			$hash{pid} = $hash{gi};
			#print "here $hash{pid}\t$hash{gi}\n";
		}
		elsif (!$hash{pid}) {
			$hash{pid} = $accession;
			$hash{pid} =~ s/\D//g;
#			print "is this ok? $hash{product} $hash{synonym} |@headerCol| have $hash{pid}\t$accession\n";
		}
#		else {
#			$hash{pid
#		}
		if (!$hash{gene}) {$hash{gene}='-'; }

		push @$annotation, \%hash; # add to annotation array
	}
#	my $cols = join "|", @headerCol;
#	die "$cols\n";

	return $genome_name;
}

sub load_sequence {
	my $file = shift;
	open (IN, $file) or die "can't open $file: $!\n";
	my $header = <IN>;
	my $seq;
	while (<IN>) {
		next if /^>/;  # NCBI genomes with multiple contigs are concatenated with NO spacer to maintain coordinants
		s/\s+//g;
		chomp;
		$seq .= $_;
	}
	$seq = uc($seq);
#	print "sequence length is " . length($seq) . "\n";	

	return $seq;
}

sub load_aa_sequence {
	my $file = shift;
	open (IN, $file) or die "can't open $file: $!\n";
	my %aa_seq;
	my $gid="";
	while (<IN>) {
		chomp;
		s/\s+//g;
		if (/gi\|(\d+)\|/) {
#			print "have $gid\t$aa_seq{$gid}\n";
			$gid = $1;
		}
		elsif (/>(.*)/) {
			$gid = $1;
			$gid =~ s/^\s+//;
			$gid =~ s/\s+$//;
			$gid =~ s/_//;  # get rid of the underscore 
		}
		else {
#			print "adding $gid\t$_\n";
			$aa_seq{$gid} .= $_;
		}
	}

	return \%aa_seq;
}
