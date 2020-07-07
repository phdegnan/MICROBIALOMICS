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

package load_CELLO;

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


sub load_CELLO_file {
	my $this=shift;
	my $CELLO_file = shift;
	my $genome_id = shift;
	$this->{GENOME_ID} = $genome_id;
	print "loading CELLO results into DB\n";

	$this->parse_CELLO($CELLO_file);

}

sub parse_CELLO {
	my ($this, $CELLOres) = @_;
	open (IN, $CELLOres) or die "can't open $CELLOres: $!\n";
	my $microbeDB = $this->{MICROBE_DB};

	my %feature_name_to_accession = (
		'inner membrane' => 'CELLO_im',
		'outer membrane' => 'CELLO_om',
		'cytoplasmic' => 'CELLO_cy',
		'periplasmic' => 'CELLO_pe',
		'extra cellular' => 'CELLO_ex',
		'cell wall'      => 'CELLO_cw',
		'membrane'      => 'CELLO_me',
	);

	my $gene_id;
	while (<IN>) {
#		chomp;
		if (/SeqID:\s+(.*)\b/) {
			$gene_id = get_gene_id($microbeDB->{DBH}, $1, $this->{GENOME_ID});
#			print "have $1\t$gene_id\n";
		}
#		elsif (/(.*)(\d+\.\d+)\*\s*$/) {
		elsif (/(.*)(\d+\.\d+)\s+\*\s*$/) {
			my $score = $2;
			my $feature_name = $1;
			$feature_name =~ s/^\s+//;
			$feature_name =~ s/\s+$//;
			die "Error cannot find CELLO accession for $feature_name\n" unless $feature_name_to_accession{$feature_name};
			my $function_id = get_function_id($microbeDB, $feature_name_to_accession{$feature_name}, $feature_name);
#			print "match $feature_name\t$score\ngid $gene_id fid $function_id\n";
			if ($function_id && $gene_id) {
				my %hash = (
					feature_id => $gene_id,
					function_id => $function_id,
					score => $score,
					eval => "",
					pval => "",
					direct_annotation => 'Y',
				);
				$microbeDB->load_feature_to_function(\%hash);
			}
		}

#                unless ($gene_id) {
#			print STDERR "could not find gene $query in feature table with genome_id $this->{GENOME_ID}, skipping...\n";
#			next;
#		}

	}
}

sub get_function_id {
	my ($microbeDB, $accession, $function_name) = @_;
	my $function_id = $microbeDB->get_function_id_from_function_accession($accession);

	return $function_id if ($function_id);
	my $function_type = "CELLO";

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

	$gname =~ s/_//;

	my $query = "SELECT feature_id FROM feature WHERE genome_id=$genome_id AND feature_symbol=?";
        my $sth = $dbh->prepare($query);
#	print "query $query\n";
        $sth->execute($gname);
        while (my $res = $sth->fetchrow_array()) {
                return $res;
        }
}

return 1;
