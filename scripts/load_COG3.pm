##
##	B. Muegge & J. Faith 2008-2009
## 	Patrick Degnan revised 2014
##	load_COG3.pm
##	Edited to read phmmer output
##

=pod
=head1 NAME

load_COG - Perl module to load COG annotations into microbialomics

=head1 SYNOPSIS

use load_COG;

my $COG = new load_COG();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package load_COG;

use strict;

#my $COGmapfile = "COG/COG.mappings.v7.1.txt";
my $COGmapfile = "/Users/Shared/DB/COG.mappings.v9.1.txt";
my $COGnamefile = "/Users/Shared/DB/COGs/whog";
my $COGcatfile = "/Users/Shared/DB/COGs/fun.txt";



sub new {
        my $proto=shift;
        my $microbeDB = shift;
        my $class = ref($proto) || $proto;
        my $self  = {
                VERSION => 0.1,
                MICROBE_DB => $microbeDB,
        };

        bless ($self, $class);
        return $self;
}


sub load_COG_file {
	my $this=shift;
	my $blastres = shift;
	my $genome_id = shift;
	$this->{GENOME_ID} = $genome_id;
	print "loading COG results into DB\n";


	# get general descriptions
	my $COGcats = parseCOGcategories($COGcatfile);
	# get specific descriptions
	my $COG = parseCOGnames($COGnamefile);
	# get map of gene to COG accession
	my $COGmap = parse_COGmapfile($COGmapfile);
	# parse the blast results for a species and put in DB
	$this->parse_blast_to_COG($blastres, $COGmap, $COG, $COGcats);
}

# get COG number (e.g. COG5653) and the corresponding description
# return a hash of with COG number as key the points to an array
# of COG category and COG description
sub parseCOGnames {
	my $in=shift;

	my %COG;
	open (IN, $in) or die "can't open $in: $!\n";
	while (<IN>) {
		chomp;
		if (/\[(\w+)\]\s+(COG\d+)\s+(.*)/) {
#		print;
			my $COGname = $2;
			my $COGcat = $1;
			my $COGdescription = $3;
			$COG{$COGname} = [$COGcat, $COGdescription];
	#		print "have $1\t$2\t$3\n";
		}
	}

	return \%COG;
}

# these are the broader category terms;
# these more general terms are useful for things like fisher's exact test
# to determine enrichment
sub parseCOGcategories {
	my $in=shift;

	open (IN, $in) or die "can't open $in: $!\n";
	my %COGcats;
	while (<IN>) {
		chomp;
		if (/\[(\w+)\]\s+(.*)/) {
			#print "$1\t$2\n";
			$COGcats{$1} = $2;
		}
	}
	close IN;

	return \%COGcats;
}


sub parse_blast_to_COG {
	my ($this, $blastres, $COGmap, $COGinfo, $COGcats) = @_;
	open (IN, $blastres) or die "can't open $blastres: $!\n";
	my $microbeDB = $this->{MICROBE_DB};

	my %res;
	my %gene_cache;
	my %function_cache;
	my %COG_cache;
	my $last_query="";##
	while (<IN>) {
		# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
		unless(/^\#/){
		chomp;
		##re vamped for phmmer
		my @temp = split /\s+/;
		my $query=$temp[2];
		my $prot=$temp[0];
		my $score=$temp[5];
		my $eval=$temp[4];
		#print "have $query\t$prot\t$score\t$eval\n";
		if (my $COG = $COGmap->{$prot}) {
			next if $COG =~ /^NOG|^KOG/;  # skip groups with no assigned function or with eukaryotic functions
			unless ($res{$query}{$COG}) {
				my $best ="Y";
				if($query eq $last_query){$best ="N";}##
				my @info = @{$COGinfo->{$COG}};
				$res{$query}{$COG} = [$score, $eval];
				my $gene_id = get_gene_id($microbeDB->{DBH}, $query, $this->{GENOME_ID});
				my ($function_id) = get_function_id($microbeDB, $COG, 0, @{$COGinfo->{$COG}});
				my %hash = (
					feature_id => $gene_id,
					function_id => $function_id,
					score => $score,
					eval => $eval,
					pval => 0,
					direct_annotation => 'Y',
					best_hit => $best,##
				);
				# don't create screen clutter and mysql error messages by loading 2x
				unless ($COG_cache{$gene_id . "_" . $function_id}) { # not sure if this is necessary, I'm writing this many months after the initial coding
					$microbeDB->load_feature_to_function(\%hash);
					load_functional_categories(\%hash, $info[0], $microbeDB, $COGcats, \%COG_cache, $gene_id);
					$COG_cache{$gene_id . "_" . $function_id}=1;
				}

#				unless ($feature_cache{$query}) { update_feature_cache(\%feature_cache, $query, $this->{GENOME_ID}); }
		                unless ($gene_id) {
					print STDERR "could not find gene $query in feature table with genome_id $this->{GENOME_ID}, skipping...\n";
					next;
				}
				$last_query=$query;##
#				print "adding $query $prot $COG $score $eval @info\n";
			}
		}
		}##end unless #
	}
}

# now load the broad functional categories
sub load_functional_categories {
	my ($hash, $function_cat, $microbeDB, $COGcats, $cache, $gene_id) = @_;

	my @cats = split //, $function_cat;
	for my $c (@cats) {
		next if $c eq 'S' or $c eq 'R';  # skip the unknown function categories
		my $accession = "COG_category_" . $c;
		my @info;
		push @info, "";
		push @info, $COGcats->{$c};

		my ($function_id, $function_cat) = get_function_id($microbeDB, $accession, 1, @info);
		#die "have cats @cats => $c to $COGcats->{$c}\n";

		# don't create screen clutter and mysql error messages by loading 2x
		unless ($cache->{$gene_id . '_' . $function_id}) { 
			$hash->{function_id} = $function_id;
			$microbeDB->load_feature_to_function($hash);

			$cache->{$gene_id . '_' . $function_id} = 1;
		}
	}
#	$microbeDB->load_feature_to_function(\%hash);


}

sub get_function_id {
	my ($microbeDB, $COG, $isCat, @info) = @_;
	my $function_id = $microbeDB->get_function_id_from_function_accession($COG);

	return $function_id if ($function_id);
	my $function_type = "COG (STRING v7.1)";
	$function_type = "COG category" if $isCat;
	
	my %fhash = (
		function_id => 0,
		function_accession => $COG,
#		function_name => $info[0],
		function_name => $info[1],
#		function_description => $info[1],
		function_description => "",  # I think it's better to not have description for COG
		function_type => $function_type, 
		
	);

	$microbeDB->load_function(\%fhash);
	
	return $fhash{function_id}; # return the function id and the functional category
}

sub get_gene_id {
        my ($dbh, $gname, $genome_id) = @_;
	$gname =~ s/_//;

	my $query = "SELECT feature_id FROM feature WHERE genome_id=$genome_id AND feature_symbol=?";
        my $sth = $dbh->prepare($query);
#	print "query $query\n";
        $sth->execute($gname);
        while (my $res = $sth->fetchrow_array()) {
                return $res;
        }
}

sub parse_COGmapfile {
	my $in=shift;

	open (IN, $in) or die "can't open $in: $!\n";

	my $count = 0;
	my $head = <IN>;
	my %map;
	print STDERR "reading COG map (will take a while)\n";
	while (<IN>) {
		my ($prot, $start, $end, $COG, $ann) = split /\t/;
		$map{$prot} = $COG;
#		print "have $prot\t$COG\n";

#		die if $count++ > 10;
	}
#	print STDERR "finished reading COG map\n";

	return \%map;

}

return 1;
