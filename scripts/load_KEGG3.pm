##
##	B. Muegge & J. Faith 2008-2009
## 	Patrick Degnan revised 2014
##	load_KEGG3.pm
##	Edited to read phmmer output
##

=pod
=head1 NAME

load_KEGG - Perl module to load KEGG annotations into microbialomics

=head1 SYNOPSIS

use load_KEGG;

my $KEGG = new load_KEGG();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package load_KEGG;

use strict;

my $KEGGfile = "/Users/Shared/DB/ko00001.keg";



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


sub load_KEGG_file {
	my $this=shift;
	my $blastres = shift;
	my $genome_id = shift;
	$this->{GENOME_ID} = $genome_id;
	print "loading KEGG results into DB\n";


	# get general descriptions
	my $KEGGfun = parseKEGGfunction($KEGGfile);
	# parse the blast results for a species and put in DB
	$this->parse_blast_to_KEGG($blastres, $KEGGfun);
}


sub parseKEGGfunction {
	my $in=shift;

	my %KEGG;
	my %KEGGentry;
	my $in_class = 0;
	my $class_data = "";
	my $cat1; my $cat2; my $path;
	open (IN, $in) or die "can't open $in: $!\n";
	while (<IN>) {
		chomp;
		if(/^A.+?\>(.+?)\</){
			$cat1=$1;
		}elsif(/^B.+?\>(.+?)\</){
			$cat2=$1;
		}elsif(/^C\s+\d+\s+(.+?)\s+\[PATH\:(ko\d+)\]/){
			$path=$1;
		}elsif(/^D\s+(K\d+)\s+(.+?)\s*\;\s+(.+)/){
			my @parents; 
			push @parents, [$cat1, $cat2, $path];
			$KEGGentry{'function_name'}=$2;
			$KEGGentry{'function_accession'}=$1;
			$KEGGentry{'function_description'}=$3;
			$KEGGentry{'parents'}=\@parents;
			#foreach my $k (keys(%KEGGentry)){print "[$KEGGentry{$k}]\t";} 
			#foreach my $p (@parents){print "[$p]\t";}print "\n";
			my %KEGGcopy = %KEGGentry;
			$KEGG{$KEGGentry{'function_accession'}} = \%KEGGcopy;
			%KEGGentry = ();	
		}

	}
	close IN;
	
	#print "K11074\n[$KEGG{'K11074'}]\n";
	#my $K1=$KEGG{'K11074'}->{'function_type'};
	#print "$K1\n";
	#die;
	return \%KEGG;
	
}



sub clean_ends {
	my $text = shift;
	$text =~ s/^\s+//;
	$text =~ s/\s+$//;

	return $text;
}

sub parse_blast_to_KEGG {
	my ($this, $blastres, $KEGGfun) = @_;
	open (IN, $blastres) or die "can't open $blastres: $!\n";
	my $microbeDB = $this->{MICROBE_DB};

	my %res;
	my %gene_cache;
	my %function_cache;
	my %COG_cache;
	my %gene_id_cache;
	my %gene_to_function_cache;  # first key gene id/queryname second key KEGGid (since it is VERY common to have multiple hits to geens with the same KO)
	my $last_query="";##
	while (<IN>) {
		# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
		chomp;
		unless(/^\#/){
		##re vamped for phmmer
	
		my @temp = split /\s+/;
		my $query=$temp[2];
		my $prot=$temp[0];
		my $score=$temp[5];
		my $eval=$temp[4];
		my $description;
		foreach my $i (18..$#temp){$description.=$temp[$i] . " ";}
		if($description=~/\;\s+(K\d+)/){ ## ignore matches that don't have K0 numbers
			my $KEGG_id = $1;
			#my ($KEGG_id, $junk) = split /##/, $prot;
			#print "have $KEGG_id\t$query\t$prot\t$score\t$eval\n";
			if (my $KEGGdata = $KEGGfun->{$KEGG_id}) {
				next if $gene_to_function_cache{$query}{$KEGG_id}; # don't try to enter things 2x
				my $gene_id = get_gene_id($microbeDB->{DBH}, $query, $this->{GENOME_ID}, \%gene_id_cache);
				die "Error no gene id found for $query\n" unless $gene_id;
				my $function_id = get_function_id($microbeDB, $KEGG_id, 0, $KEGGdata->{'function_name'}, $KEGGdata->{'function_description'},\%function_cache);
				#print "have [$KEGG_id]\t[$KEGGdata->{'function_name'}]\t[$KEGGdata->{'function_description'}]\n";
#				next if $COG =~ /^NOG|^KOG/;  # skip groups with no assigned function or with eukaryotic functions
#				print "inner data for gene $query id $gene_id function id $function_id\n";
				#die "Error no function id found for $KEGG_id\n" unless $function_id;
				if($function_id){# 
					##$KEGG_id is a brite pathway or updated since 2010
					my $best ="Y";##
					if($query eq $last_query){$best ="N";}##
					load_feature_to_function($gene_id, $function_id, $score, $eval, 'Y', $query, $KEGG_id, $best, \%gene_to_function_cache, $microbeDB);
					$last_query=$query;##
					# load the higher level categories now	
					for my $c (@{$KEGGdata->{parents}}) {
						for my $i (0 .. $#$c) {
							my $function_name = $c->[$i];
							$function_name =~ s/\s+/ /g;
							$function_name =~ s/\s*\[.*?\]//g;
							my $function_accession = get_KEGG_function_class($i+1, 1) . '_' . $function_name;
							$function_accession =~ s/\s+/_/g;
							my $function_description = "";
						
							my $function_id = get_function_id($microbeDB, $function_accession, $i+1, $function_name, $function_description, \%function_cache);
							next if $gene_to_function_cache{$query}{$function_accession}; # don't try to enter things 2x
							load_feature_to_function($gene_id, $function_id, $score, $eval, 'Y', $query, $function_accession, $best,  \%gene_to_function_cache, $microbeDB);
						
#							print "$KEGG_id\t|$function_accession|\t|$function_name|$c->[$i]|\n";
						}
#					$KEGGentry{parents} = parse_KEGG_class($class_data);
					}
				}
#					load_functional_categories(\%hash, $info[0], $microbeDB, $COGcats, \%COG_cache, $gene_id);
			}
		}
		}#end unless #
	}#end while
}

sub load_feature_to_function {
	my ($gene_id, $function_id, $score, $eval, $direct_annotation, $query, $KEGG_id, $best, $cache, $microbeDB) = @_;
#			load_feature_to_function($gene_id, $function_id, $score, $eval, 'Y', $query, $KEGG_id, \%gene_to_function_cache);
	my %hash = (
		feature_id => $gene_id,
		function_id => $function_id,
		score => $score,
		eval => $eval,
		pval => 0,
		direct_annotation => $direct_annotation,
		best_hit => $best,##
	);
	$microbeDB->load_feature_to_function(\%hash);
	$cache->{$query}{$KEGG_id} = 1;
}

sub get_KEGG_function_class {
	my ($class, $short) = @_;

	my $function_type;
	if ($short) {
		$function_type = "KEGGcat$class";
		$function_type = "KEGGpath" if $class == 3;
	}
	else {
		$function_type = "KEGG category $class";
		$function_type = "KEGG pathway" if $class == 3;
	}

	return $function_type;
}

sub get_function_id {
	my ($microbeDB, $KEGG_id, $isClass, $name, $description, $cache) = @_;

	my $function_id;
	if ($cache->{$KEGG_id}) { # move on if it's in the cache
		return $cache->{$KEGG_id};
	}
	$function_id = $microbeDB->get_function_id_from_function_accession($KEGG_id);

	return $function_id if ($function_id);
	my $function_type = "KEGG";
	if ($isClass) {
		$function_type = get_KEGG_function_class($isClass);
#		$function_type = "KEGG category $isClass";
#		$function_type = "KEGG pathway" if $isClass == 3;
	}

	my %fhash = (
		function_id => 0,
		function_accession => $KEGG_id,
#		function_name => $info[0],
		function_name => $name,
#		function_description => $info[1],
		function_description => $description,  # I think it's better to not have description for COG
		function_type => $function_type, 
	);

	$microbeDB->load_function(\%fhash);
	$cache->{$KEGG_id} = $fhash{function_id};
	
	return $fhash{function_id}; # return the function id and the functional category
}

sub get_gene_id {
        my ($dbh, $gname, $genome_id, $cache) = @_;
	$gname =~ s/_//;

	if ($cache && $cache->{$gname}) {
		return $cache->{$gname};
	}

	my $query = "SELECT feature_id FROM feature WHERE genome_id=$genome_id AND feature_symbol=?";
        my $sth = $dbh->prepare($query);
#	print "query $query\n";
        $sth->execute($gname);
        while (my $res = $sth->fetchrow_array()) {
		$cache->{$gname} = $res;
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
__END__

sub parseKEGGfunction {
	my $in=shift;

	my %KEGG;
	my %KEGGentry;
	my $in_class = 0;
	my $class_data = "";
	open (IN, $in) or die "can't open $in: $!\n";
	while (<IN>) {
		chomp;
		if (/^ENTRY/) {
			my ($name, $id, $type) = split /\s+/;
#			print "class: $class_data\n";
			if (%KEGGentry) {
				$KEGGentry{parents} = parse_KEGG_class($class_data);
				my %KEGGcopy = %KEGGentry;
				$KEGG{$KEGGentry{function_accession}} = \%KEGGcopy;
			}

			%KEGGentry = ();
			$class_data = "";
			$KEGGentry{'function_accession'} = $id;
		}
		elsif (/^NAME/) {
			s/^NAME\s+//;
			$KEGGentry{'function_name'} = $_;
		}
		elsif (/^DEFINITION/) {
			s/^DEFINITION\s+//;
#			print "description $_\n";
			$KEGGentry{'function_description'} = $_;
		}
		elsif (/^DBLINKS/ || /^GENES/ || /^\/\/\//) {
			$in_class=0;
		}
		elsif (/^CLASS/ || $in_class) {
			s/^CLASS\s+//;
			$in_class=1;
			$class_data .= $_;
		}
#		print;

	}
	close IN;

	return \%KEGG;
}

sub parse_KEGG_class {
	my $class = shift;
	my @pieces = split /\[PATH:ko\d+\]/, $class;
#	print "starting with $class\n";
#	print "have pieces @pieces\n";
#	($KEGGentry{class});
	my @parents;
	for my $p (@pieces){
#		print "splitting $p\n";
		my ($level1, $level2, $path) = split /;/, $p;
		unless ($level1 && $level2 && $path) {
			die "KEGG parse error for $p class $class\n";
		}
		$level1 = clean_ends($level1);
		$level2 = clean_ends($level2);
		$path = clean_ends($path);
		push @parents, [$level1, $level2, $path];
	}

	return \@parents;
}
