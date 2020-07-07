=pod
=head1 NAME

load_CELLO - Perl module to load CELLO annotations into microbialomics

=head1 SYNOPSIS

use load_COG;

my $COG = new load_COG();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package load_CAZy;

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


sub load_CAZy_file {
	my $this=shift;
	my $CAZy_file = shift;
	my $genome_id = shift;
	$this->{GENOME_ID} = $genome_id;
	print "loading CAZy results into DB\n";

	$this->parse_CAZy($CAZy_file);

}

sub parse_CAZy {
	my ($this, $CELLOres) = @_;
	open (IN, $CELLOres) or die "can't open $CELLOres: $!\n";
	my $microbeDB = $this->{MICROBE_DB};

	my $gene_id;
	while (<IN>) {
	#	print STDERR;
		chomp;
		my ($gene_name, $cazy_terms) = split /\t/;
		$gene_id = get_gene_id($microbeDB->{DBH}, $gene_name, $this->{GENOME_ID});
		next unless ($gene_id); # skip if we can't find gene

		my @terms = split /-/, $cazy_terms;
		for my $t (@terms) {
			my $function_id = get_function_id($microbeDB, "CAZy" . '_'  . $t, $t);
			if ($function_id && $gene_id) {
#				print "loading $gene_name\t$t\n";
				my %hash = (
					feature_id => $gene_id,
					function_id => $function_id,
					score => "",
					eval => "",
					pval => "",
					direct_annotation => 'Y',
				);
				$microbeDB->load_feature_to_function(\%hash);
			}
		}

	}
}

sub get_function_id {
	my ($microbeDB, $accession, $function_name) = @_;
	my $function_id = $microbeDB->get_function_id_from_function_accession($accession);

	return $function_id if ($function_id);
	my $function_type = "CAZy";

	my %fhash = (
		function_id => "",
		function_accession => $accession,
#		function_name => $info[0],
		function_name => $function_name,
#		function_description => $info[1],
		function_description => "",  # I think it's better to not have description for COG
		function_type => $function_type, 
	);

	$microbeDB->load_function(\%fhash);
	
	return $fhash{function_id}; # return the function id and the functional category
}

sub get_gene_id {
        my ($dbh, $gname, $genome_id) = @_;

	my $query = "SELECT feature_id FROM feature WHERE genome_id=$genome_id AND feature_symbol=?";
        my $sth = $dbh->prepare($query);
#	print "query $query\n";
        $sth->execute($gname);
        while (my $res = $sth->fetchrow_array()) {
                return $res;
        }
}

return 1;
