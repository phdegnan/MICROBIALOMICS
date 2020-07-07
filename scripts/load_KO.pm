=pod
=head1 NAME

load_COG - Perl module to load knockouts into the database

=head1 SYNOPSIS

use load_COG;

my $KO = new load_KO();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package load_KO;

use strict;

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


sub load_KO_file {
	my $this=shift;
	my $KOfile = shift;
	my $genome_id = shift;
	$this->{GENOME_ID} = $genome_id;
	print "loading KO results into DB\n";


	my $KO = $this->parseKOfile($KOfile);
}

sub parseKOfile {
	my $this=shift;
	my $file=shift;

	open (IN, $file) or die "can't open $file: $!\n";
	my $header = <IN>;
	
	my $microbeDB = $this->{MICROBE_DB};
	while (<IN>) {
		chomp;
		my @cols = split /\t/;
		# get the basics
		my ($direct_gene, $KOhash) = get_KO_hash(@cols);
		my @polar_genes;
		@polar_genes = split /::/, $cols[7] if $cols[7];
#		print "polar hits are @polar_genes\n" if @polar_genes;
		# now load direct hit
		load_gene_hit($direct_gene, $KOhash, $this->{GENOME_ID}, $microbeDB->{DBH}, $microbeDB);
		load_polar_hits(\@polar_genes, $KOhash, $this->{GENOME_ID}, $microbeDB->{DBH}, $microbeDB);
		# then load polar hits
	}
}

sub load_polar_hits {
	my ($polar_genes, $KOhash, $genome_id, $dbh, $microbeDB) = @_;
	$KOhash->{hit_type} = "polar";
	for my $g (@$polar_genes) {
		load_gene_hit($g, $KOhash, $genome_id, $dbh, $microbeDB);
	}
}

sub load_gene_hit {
	my ($direct_gene, $KOhash, $genome_id, $dbh, $microbeDB) = @_;
	my $gid = get_gene_id($dbh, $direct_gene, $genome_id);

	unless ($gid) {
		print STDERR "warning couldn't find feature id for $direct_gene\n";
		return;
	}
	$KOhash->{feature_id} = $gid;
	$microbeDB->load_genetics($KOhash);
}

sub get_KO_hash {
	my ($coordinate, $strainID, $multiple_matches, $plate, $well, $mismatches, $direct) = @_;
	my ($direct_g, $percent) = split /;/, $direct;
	
	my $description = "hits $direct_g $percent from start codon; strain is in plate $plate well $well";
	if ($mismatches) { $description .= "; barcode matches with $mismatches mismatches"; }
	my $strain_id = "AG_BTv1_" .  $strainID;

	my %hash = (
		hit_type => "direct",
		perturbation_type => "knockout",
		strain_id => $strain_id, 
		multiple_matches => $multiple_matches,
		insertion_pos => $coordinate,
		description => $description,
	);
#	print "strain $strain_id desc $description multiple matches ($multiple_matches)\n";
#	print "$description\n";

	return $direct_g, \%hash;
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
