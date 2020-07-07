=pod
=head1 NAME

microbialomicsDB - Perl module to deal with the database stuff for the microbialomics database

=head1 SYNOPSIS

use microbialomicsDB;

my $fe = new microbialomicsDB();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package microbialomicsDB;

use strict;
use DBI;

sub new {
	my $proto = shift;
	my $db_name = shift;
	my $class = ref($proto) || $proto;
	my $host = shift;
	my $user = shift;
	my $pass= shift;
	my $dbh = DBI->connect("DBI:mysql:$db_name;mysql_socket=/var/mysql/mysql.sock",$user,$pass) or die "can't open database: $DBI::errstr\n";

	my $self  = {
		VERSION => 0.1,
		DBH => $dbh,
		USER => $user,
		PASS => $pass,
		DB => $db_name,
		GENOME_INFO => {},
	};

	bless ($self, $class);
	return $self;
}

sub clearOldDatabase {
	my $this=shift;
	my $sql_schema = "create_microbialomics.sql";
	my $query = "mysql -u$this->{USER} -p$this->{PASS} -D$this->{DB} < $sql_schema";
	my $project_url1 = "http://spreadsheets.google.com/pub?key=rQrJjKgf3qBm9wuWDoa8trg&single=true&gid=0&output=txt";
	my $project_url2 = "http://spreadsheets.google.com/pub?key=rQrJjKgf3qBm9wuWDoa8trg&single=true&gid=2&output=txt";


	$this->add_genome_info($project_url1);
	$this->add_genome_info($project_url2);
	

#	print "running $query\n";
	system("$query");
}

sub add_genome_info {
	my $this=shift;
	my $url = shift;

	my $table = "table.tab";
	download_table($url, "table.tab");

	open (IN, $table) or die "can't open $table: $!\n";

	my $header = <IN>;
	chomp $header;
	my @headers = split /\t/, $header;
	while (<IN>) {
		chomp;
		my @values = split /\t/, $_;

		my %data;
		for my $i (0 .. $#values) {
			$data{lc($headers[$i])} = $values[$i];
#			print "$headers[$i]\t$values[$i]\n";
		}
		if ($data{'sequenced'}) {
			$this->{GENOME_INFO}{$data{'sequenced'}} = \%data;
			print "storing $data{sequenced}\n";
		}
	}
#	die;
}

sub download_table {
	my ($url, $outfile) = @_;

	print STDERR "\n\ndownloading $url to $outfile\n";
	my $download_info = `wget '$url' -O $outfile -t 3`;
	sleep(1);

}

sub load_genome_table {
	my ($this, $genome_hash) = @_;

	my $g = $genome_hash; # less typing

	my $dbh = $this->{DBH};
	print STDERR "loading genome $g->{genome_name} into DB\n";
	#my $sth = $dbh->prepare("INSERT INTO genome VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
	my $sth = $dbh->prepare("INSERT INTO genome VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)");## database has a lot fewer columns
	$sth->execute($g->{genome_id}, $g->{genome_constant_id}, $g->{genome_name}, $g->{genome_description}, $g->{genus_species},
	              $g->{genome_length}, $g->{genome_sequence}, $g->{sequencing_date},
	              $g->{pmid}, $g->{sequencing_method}, $g->{assembly_method}, $g->{aerobic_tolerance}, $g->{growth_conditions}, $g->{public_genome}
		);

#	             ,$g->{symbiosis}, $g->{gram_stain}, $g->{isolated_from}, $g->{host_phenotype}, $g->{host}, $g->{host_age_in_years},
#	              $g->{host_sex}, $g->{host_origin}, $g->{literature}, $g->{expression_db}
#		);

	$g->{genome_id} = $sth->{'mysql_insertid'};
}

sub load_contig {
	my ($this, $contig_hash) = @_;

	my $c = $contig_hash;
	my $contig_id = 0;
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("INSERT INTO contig VALUES(?,?,?,?,?)");

	$sth->execute($contig_id, $c->{genome_id}, $c->{contig_description}, $c->{start_pos}, $c->{stop_pos});
}

sub load_genetics {
	my ($this, $g) = @_;

	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("INSERT INTO genetics VALUES(?,?,?,?,?,?,?)");

	$sth->execute($g->{feature_id}, $g->{perturbation_type}, $g->{hit_type}, $g->{strain_id}, $g->{description}, $g->{multiple_matches}, $g->{insertion_pos});	
}

sub load_aa_seqs {
	my ($this, $seqs) = @_;

	my $dbh = $this->{DBH};
	print STDERR "loading protein sequences into DB\n";
	my $sth = $dbh->prepare("INSERT INTO protein_seq VALUES(?,?)");

	for my $s (@$seqs) {
		$sth->execute($s->{feature_id}, $s->{protein_sequence});	
	}
}

sub load_function {
	my ($this, $fun) = @_;
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("INSERT INTO function VALUES(?,?,?,?,?)");
	$sth->execute($fun->{function_id}, $fun->{function_accession}, $fun->{function_name}, $fun->{function_description}, $fun->{function_type});
	$fun->{function_id} = $sth->{'mysql_insertid'};
#	print "added new function id for $fun->{function_name} is $fun->{function_id}\n";
}

sub load_feature_to_function {
	my ($this, $ftf) = @_;
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("INSERT INTO feature_to_function VALUES(?,?,?,?,?,?,?)");
	$sth->execute($ftf->{feature_id}, $ftf->{function_id}, $ftf->{score}, $ftf->{eval}, $ftf->{pval}, $ftf->{direct_annotation},$ftf->{best_hit});
}



sub load_annotation_table {
	my ($this, $ann) = @_;
	my $dbh = $this->{DBH};
	print STDERR "loading feature annotations into DB\n";
	my $sth = $dbh->prepare("INSERT INTO feature VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)");

	for my $a (@$ann) {
		$sth->execute($a->{feature_id}, $a->{genome_id}, $a->{genbank_id}, $a->{feature_name}, $a->{feature_symbol}, $a->{feature_description},
		              $a->{feature_type}, $a->{feature_evidence}, $a->{feature_score}, $a->{start_pos}, $a->{stop_pos},
		              $a->{direction}, $a->{primary_feature});

		$a->{feature_id} = $sth->{'mysql_insertid'};
	}

}

sub load_project_table {
	my ($this, $project_information) = @_;
	
	#$p->{project_id} = $sth->{'mysql_insertid'};
}

sub load_feature_set_table {
	my ($this, $feature_sets, $features) = @_;

	my $dbh = $this->{DBH};
#	print STDERR "loading feature sets into DB\n";
	my $sth = $dbh->prepare("INSERT INTO feature_sets VALUES(?,?,?)");

	my @keys = sort {$a cmp $b} keys %$feature_sets;

	for my $k (@keys) {
		my $feature_array = $feature_sets->{$k};

		for my $f (@$feature_array) {
#			print "loading $k\t$features->{$f->[0]}{feature_id}\t$f->[1]\n";
			$sth->execute($k, $features->{$f->[0]}{feature_id}, $f->[1]);
		}
	}
}

sub load_conditions {
	my ($this, $condition_info, $features, $experiment_id, $condition_skip) = @_;

	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("INSERT INTO conditions VALUES(?,?,?,?)");

	my @keys = sort {$a cmp $b } keys %$condition_info;

	#print "have @keys\n";
	my @k = keys %$condition_skip;

	for my $k (@keys) {
#		print "have $k\t$condition_info->{$k}\n";
		next, print "skip\n" if $condition_skip->{$k}; # skip columns already used for basic experimental table stuff
			
		my $feature_id = $features->{$k}{feature_id};
		# for numeric stuff store the number too
		if ($features->{$k}{type} eq 'integer' || $features->{$k}{type} eq 'real') {
			$sth->execute($experiment_id, $feature_id, $condition_info->{$k}, $condition_info->{$k});
		}
		else { # otherwise just store the text
			my $rv = $sth->execute($experiment_id, $feature_id, $condition_info->{$k}, "");
			print join("|",($experiment_id,$feature_id,$condition_info->{$k})),"\n" if (!$rv);
		}	
	}
}

# BELOW ARE BASIC RETRIEVAL FUNCTIONS
sub get_function_id_from_function_accession {
	my ($this, $accession) = @_;

	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("SELECT function_id from function where function_accession=?");
	$sth->execute($accession);
	while (my $res = $sth->fetchrow_array()) {
		return $res;
	}
}

sub get_function_id_from_function_name {
	my ($this, $name) = @_;

	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("SELECT function_id from function where function_name=?");
	$sth->execute($name);
	while (my $res = $sth->fetchrow_array()) {
		return $res;
	}
}

sub get_genome_id_from_genome_name {
	my ($this, $gname) = @_;

	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("SELECT genome_id from genome where genome_name=?");
	$sth->execute($gname);
	while (my $res = $sth->fetchrow_array()) {
		return $res;
	}
}



sub get_feature_id_from_feature_symbol {
	my ($this, $symbol) = @_;

	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare("SELECT feature_id from feature where feature_symbol=?");
	$sth->execute($symbol);
	while (my $res = $sth->fetchrow_array()) {
		return $res;
	}
}

1;
