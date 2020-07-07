=pod
=head1 NAME

load_COG - Perl module to load functional annotations based on expression into the database

=head1 SYNOPSIS

use load_express;

my $express = new load_express();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package load_express;

use strict;

sub new {
        my $proto=shift;
        my $microbeDB = shift;
        my $class = ref($proto) || $proto;
        my $self  = {
                VERSION => 0.1,
                MICROBE_DB => $microbeDB,
		MIN_PVAL => 0.05,
        };

        bless ($self, $class);
        return $self;
}


sub load_express_file {
	my $this=shift;
	my $KOfile = shift;
	my $genome_id = shift;
	$this->{GENOME_ID} = $genome_id;

	my $KO = $this->parse_express_file($KOfile);
}

sub parse_express_file {
	my $this=shift;
	my $file=shift;

	open (IN, $file) or die "can't open $file: $!\n";
	my $header = <IN>;
	my @genes = split /\t/, $header;
	shift @genes;
	my $gene_to_feature_id = $this->get_gene_to_feature_id(@genes);
#	print "have genes @genes\n";
	
	my $min_pval = $this->{MIN_PVAL};
	my $microbeDB = $this->{MICROBE_DB};
	while (<IN>) {
		chomp;
		my ($condition, @data) = split /\t/;
		for my $i (0 .. $#data) {
			my $d = $data[$i];
			my $g = $genes[$i];
			my ($correlation, $pval) = split /:/, $d;
			if ($pval <= $min_pval) {
				my $gid = $gene_to_feature_id->{$g};
				my $function_id = get_function_id($microbeDB, $condition);
				if ($gid && $function_id) {
					load_annotation($microbeDB, $gid, $function_id, $correlation, $pval);
				}
				
#				print "have $d\t$correlation\t$pval\tfor $g\t";
#				print "have gid $gid\tfunction id $function_id\n";
			}
		}
	}
}

sub load_annotation {
	my ($microbeDB, $gene_id, $function_id, $score, $pval) = @_;

	my %hash = (
		feature_id => $gene_id,
		function_id => $function_id,
		score => $score,
		eval => $pval,  # maybe I should use pval?
		pval => "",
		direct_annotation => 'Y',
	);
	$microbeDB->load_feature_to_function(\%hash);
}

sub get_function_id {
        my ($microbeDB, $function_name) = @_;
        my $function_id = $microbeDB->get_function_id_from_function_name($function_name);

        return $function_id if ($function_id);

        my %fhash = (
                function_id => "",
                function_accession => $function_name,
                function_name => $function_name,
                function_description => "",
                function_type => "expression_modulation",
        );

        $microbeDB->load_function(\%fhash);

        return $fhash{function_id};
}



sub get_gene_to_feature_id {
	my ($this, @genes) = @_;
	my $microbeDB = $this->{MICROBE_DB};

	my %gene_to_feature_id;
	for my $g (@genes) {
		my $probe_name = $g;
#		$g =~ s/_.*?$//g;
		my @pieces = split /_/, $probe_name;

		for my $p (@pieces) {
			next if (length($p) <= 2);
			next if $p !~ /\w\d{4}/;  # use the gene symbols
			my $gid = get_gene_id($microbeDB->{DBH}, $p, $this->{GENOME_ID});
			if ($gid) {
				$gene_to_feature_id{$probe_name} = $gid;
		#		print STDERR "have $g\t$gid\n";
			}
			else {
				print STDERR "can't find $p\n";
			}
		}
	}
	
	return \%gene_to_feature_id;
}

sub get_KO_hash {
	my %hash = (
		hit_type => "direct",
		perturbation_type => "knockout",
	);
}

sub get_gene_id {
        my ($dbh, $gname, $genome_id) = @_;

        my $sth = $dbh->prepare("SELECT feature_id FROM feature WHERE genome_id=$genome_id AND feature_symbol=?");
        $sth->execute($gname);
        while (my $res = $sth->fetchrow_array()) {
                return $res;
        }
}


return 1;
