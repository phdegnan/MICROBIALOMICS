#!/usr/bin/perl -w
##
##	B. Muegge & J. Faith 2008-2009
## 	Patrick Degnan revised 2014
##	load_pfam3.pm
##	Edited to read hmmrscan output
##
=pod
=head1 NAME

microbialomicsDB - Perl module to load pfam results into the microbialomics database

=head1 SYNOPSIS

use load_pfam;

my $pfam = new load_pfam($microbeDB);

$pfam->load_pfam_file($file, $genome_id);

=head1 DESCRIPTION

=head2 METHODS

=over

=cut

## edited to read hmmscan results from hmmer3 Phd

#my $hmminfo_dir = "/Users/faithj/tigrfam/";
# go information http://search.cpan.org/~cmungall/go-perl/go-perl.pod

package load_pfam;

use strict;

# was too slow to parser this over and over again; moved to a global position
#use GO::Parser;
#my $GOparser = new GO::Parser({handler=>'obj'});

my %GO_cache;

sub new {
	my $proto=shift;
	my $microbeDB = shift;
	my $GOparser = shift;
        my $class = ref($proto) || $proto;
        my $self  = {
                VERSION => 0.1,
		MICROBE_DB => $microbeDB,
		HMM_INFO_DIR => "/Users/Shared/DB/",
		GO_PARSER => $GOparser,
		#HMM_INFO_DIR => "pfam/infofiles/",
        };

        bless ($self, $class);
        return $self;
}


sub load_pfam_file {
	my $this=shift;
	my $pfam_file = shift;
	my $genome_id = shift;

	print "loading pfam results into DB\n";
	
	$this->{GENOME_ID} = $genome_id;	
	my $ann = $this->read_pfam_results($pfam_file);
	$this->read_pfam_definitions();
#	$this->read_GO();
#	die "finished reading definitions\n";
	$this->add_pfam_results($ann);
}

sub read_GO {
	my $this=shift;

	my $file = $this->{HMM_INFO_DIR} . "gene_ontology.obo";
	print STDERR "reading go\n";
	#$GOparser->parse($file);
#	die "read go parser\n";
}

sub read_pfam_definitions {
	my $this=shift;
	my $dir = $this->{HMM_INFO_DIR};
	my $file = $dir . "accession_to_description.txt";

	# get the definitions
	open (IN, $file) or die "can't open $file: $!\n";
	my %acc_to_def;
	while (<IN>) {
		chomp;
		my ($accession, $name, $description) = split /\t/;
		$acc_to_def{$accession} = [$name, $description];
	}

	$this->{ACCESSION_TO_DEFINITION} = \%acc_to_def;

	# now get the GO associations for those definitions
	$file=$dir . "HMM_TO_GO.txt";
	my %acc_to_GO;
	open (IN, $file) or die "can't open $file: $!\n";
	while (<IN>) {
		chomp;
		my ($acc, $GO) = split /\t/;
		push @{$acc_to_GO{$acc}}, $GO;
#		print "have $acc\t$GO\n";
	}
	$this->{ACCESSION_TO_GO} = \%acc_to_GO;
}


sub read_pfam_results {
	my ($this, $pfam_file) = @_;
	open (IN, $pfam_file) or die "can't open $pfam_file: $!\n";
	my $in_targets = 0;
	my $query_gene;

	my %gene_hits;
	while (<IN>) {
#		print;	
		chomp;
		if (/Query:\s+(\S+)/) {
#			print "have $1\n";
			$query_gene = $1;
		}
		elsif ($in_targets) {
			if (/^\s*$/ || /No hits/) {
				$in_targets=0;
			}
			else {
				my @temp = split /\s+/;
				my $accession = $temp[9];
				my $eval = $temp[1];
				my $score = $temp[2];
#				print "have $accession, $score, $eval for $query_gene\n";
				$accession =~ s/^(PF\d+)\.\d+/$1/;
#				print "after have $accession, $score, $eval for $query_gene\n\n";
				
				push @{$gene_hits{$query_gene}}, [$accession, $score, $eval];
			}
		}
		elsif (/E-value\s+score\s+bias/) {
			my $trash=<IN>;
			$in_targets=1;
		}
	}
	return \%gene_hits;
}

sub add_pfam_results {
	my $this = shift;
	my $pfam = shift;
	my $microbeDB = $this->{MICROBE_DB};
	my $GENOME_ID = $this->{GENOME_ID};

	my %function_info;
	my @genes = sort {$a cmp $b } keys %$pfam;
#	printf "have pfam hits for %d genes\n", scalar @genes;
	my $num_possible_genes = scalar @genes;
	my $passed_threshold = 0;
	my $domain_hits = 0;
	my $genes_passed = 0;
	my $genes_passed_sum = 0;
	my $accession_to_GO = $this->{ACCESSION_TO_GO};

	for my $g (@genes) {
		$genes_passed = 0;
		my $hits = $pfam->{$g};
		my $gene_id = get_gene_id($microbeDB->{DBH}, $g, $GENOME_ID);
		unless ($gene_id) {
			print STDERR "could not find gene $g in feature table with genome_id $GENOME_ID, skipping...\n";
			next;
		}

		for my $h (@$hits) {
			$domain_hits++;
			my $accession = $h->[0];
#			print "have $g\t$h->[0]\n";
			# !!!!!!!!! ADD IF SCORE > NOISE THRESHOLD DEFINED BY TIGRFAM INFO FILE
			my $dinfo = $this->get_function_info(\%function_info, $accession);

			# keep only hits that are higher than the tigrfam defined noised cutoff
			
#			print "checking $g\t$h->[0] score ($h->[1]) less than noise cutoff $dinfo->{noise_cutoff}\n";
#			if ($h->[1] > $dinfo->{noise_cutoff}) {
			$this->add_function($dinfo);
			$this->add_feature_to_function($dinfo, $gene_id, $h);
			$passed_threshold++;
			$genes_passed = 1;
			if ($accession_to_GO->{$accession}) {
				$this->process_GO($gene_id, $accession_to_GO->{$accession}, \%function_info);
			}
			
#			}
#			else {
#				print "skipping $g\t$h->[0] score ($h->[1]) less than noise cutoff $dinfo->{noise_cutoff}\n";
#			}
#			add_feature_to_function($dinfo);
		}
		$genes_passed_sum++ if $genes_passed;
	}
	print "$passed_threshold of $domain_hits initial matches with scores higher than noise in $genes_passed_sum of the $num_possible_genes genes with initial matches\n";
}

# note that accession_to_GO is no longer necessary
sub process_GO {
	my ($this, $gene_id, $accession_to_GO, $function_info) = @_;

	my %function_info;
	my @hit_info = (0, 0, 0);
	my $indirect = 1;
	my $GOparser = $this->{GO_PARSER};

#	my %already_added_cache;
	for my $GO (@$accession_to_GO) {
		# insert the direct terms first
		my $graph = $GOparser->handler->graph;
		my $term = $graph->get_term($GO);
		next unless $term; # skip terms that no longer map; these are due to term-deprecation

		if ($GO_cache{$gene_id . '_' . $GO}) { # don't need to add the same GO term 3x
	#		print STDERR "inner cache $GO\n"; next;
			next;
		}
#		if ($GO_cache{$gene_id . '_' . $GO_parent_accession}) { # don't need to add the same GO term 3x

#		}
		my $dinfo = $this->get_function_info($function_info, $GO);
		$this->add_function($dinfo);
		$this->add_feature_to_function($dinfo, $gene_id, \@hit_info, $indirect);
		$GO_cache{$gene_id . '_' . $GO} = 1;

		# now do all of the parents of each term
		my $ancestor_terms = $graph->get_recursive_parent_terms($term->acc);
		for my $GO_parent (@$ancestor_terms) {
			my $GO_parent_accession = $GO_parent->acc;
			if ($GO_cache{$gene_id . '_' . $GO_parent_accession}) { # don't need to add the same GO term 3x
#				print STDERR "in cache with $GO_parent_accession\n"; next;
				next;
			}

			my $dinfo = $this->get_function_info($function_info, $GO_parent_accession);
			$this->add_function($dinfo);
			$this->add_feature_to_function($dinfo, $gene_id, \@hit_info, $indirect);

#			$already_added_cache{$GO_parent_accession} = 1;
			$GO_cache{$gene_id . '_' . $GO_parent_accession} = 1;
		}
	}
}

sub add_feature_to_function {
	my ($this, $fun, $feature_id, $hit_array, $indirect) = @_;
	my $microbeDB = $this->{MICROBE_DB};
#	print "function id for $fun->{function_name} is $fun->{function_id}\n";
	my $best="Y";
	my %hash = (
		feature_id => $feature_id,
		function_id => $fun->{function_id},
		score => $hit_array->[1],
		eval => $hit_array->[2],
		pval => 0,
		direct_annotation => 'Y',
		best_hit => $best,
	);
	$hash{direct_annotation} = 'N' if $indirect;

	$microbeDB->load_feature_to_function(\%hash);
}

sub get_gene_id {
	my ($dbh, $gname, $genome_id) = @_;
	$gname =~ s/_//;

	my $sth = $dbh->prepare("SELECT feature_id FROM feature WHERE genome_id=$genome_id AND feature_symbol=?");
	$sth->execute($gname);
	while (my $res = $sth->fetchrow_array()) {
		return $res;
	}
}


sub add_function {
	my $this=shift;
	my $fhash=shift;
	my $microbeDB = $this->{MICROBE_DB};

#	print "adding $fhash->{function_accession}\n";	
	if ($fhash->{function_id}) { # already in the cache
#		print "id in cache\n";
		return;
	}
	
	$fhash->{function_id} = $microbeDB->get_function_id_from_function_accession($fhash->{function_accession});
	if ($fhash->{function_id}) { # already in the database
#		print "id in db\n"; 
		return;
	}

	# otherwise add it to the database
	$microbeDB->load_function($fhash);
#	print "id is now $fhash->{function_id}\n";	
}


sub get_function_info {
	my ($this, $infoCache, $accession) = @_;


	my $hmminfo_dir = $this->{HMM_INFO_DIR};
	my $infoFile = $hmminfo_dir . $accession . ".INFO";
	my $GOparser = $this->{GO_PARSER};
#	print "looking for $infoFile\n";

	my $acc_to_def = $this->{ACCESSION_TO_DEFINITION};

	# go get it if we don't have it already
	if (!$infoCache->{$accession}) {
	#	print "using file\n";
#		if (-e $infoFile) {
		if ($accession =~ /^GO/) {
			my %info;
			$info{function_type} = "GO";
			my $graph = $GOparser->handler->graph;
			$info{function_accession} = $accession;
#			print "have function accession $accession\n";
			my $term = $graph->get_term($accession);
			$info{function_name} = $term->name;
			$info{function_description} = "";
#			$info{function_description} = $term->description;
			$infoCache->{$accession} = \%info;
		}
		elsif ($acc_to_def->{$accession}) {
			my %info;
			$info{function_name} = shift @{$acc_to_def->{$accession}};
			$info{function_description} = shift @{$acc_to_def->{$accession}};
			if ($accession =~ /^PF/) {
				$info{function_type} = "PFAM";
			}
			else {
				$info{function_type} = "TIGRFAM";
			}
			$info{function_accession} = $accession;


#			my $infoData = parse_info_file($infoFile, $accession);
			$infoCache->{$accession} = \%info;
		}
		else {
			#die "unknown domain $accession\n";
			print STDERR "unknown domain $accession\n";
		}
	}
#	else { print "using cache\n";}

	return $infoCache->{$accession};
}

sub parse_info_file {
	my $in = shift;
	my $accession = shift;

	my %info;
	open (IN, $in) or die "can't open $in: $!\n";
	while (<IN>) {
#		print;
		chomp;
		if (/^DE/) {
			s/^DE\s+//;
#			print "DE is $_\n";
			$info{function_name} = $_;
		}
		elsif (/^CC/) {
			s/^CC\s+//;
#			print "CC is $_\n";
			$info{function_description} .= $_;
		}
		elsif (/^NC/) {
			s/^NC\s+//;
			my ($score, @junk) = split /\s+/;
			$info{noise_cutoff} = $score;
#			print "NC is $score\n";
		}
	}
	if ($info{function_description}) {
		$info{function_description} =~ s/^\s+//;
		$info{function_description} =~ s/\s+$//;
	}
	$info{function_type} = "TIGRFAM";
	$info{function_accession} = $accession;
	return \%info;
}

#$microbeDB->load_annotation_table($annotation);

return 1;
