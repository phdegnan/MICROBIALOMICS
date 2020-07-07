=pod
=head1 NAME

microbialomicsDB - Perl module to deal with the database stuff for the microbialomics database

=head1 SYNOPSIS

use microbialomics;

my $fe = new microbialomics();

=head1 DESCRIPTION

description here

=head2 METHODS

=over

=cut

package microbialomics;

use strict;
use DBI;

my $TRANS_TABLE = {
        TTT => 'F', TCT => 'S', TAT => 'Y', TGT => 'C', TTC => 'F', TCC => 'S',
        TAC => 'Y', TGC => 'C', TTA => 'L', TCA => 'S', TAA => '*', TGA => '*',
        TTG => 'L', TCG => 'S', TAG => '*', TGG => 'W', CTT => 'L', CCT => 'P',
        CAT => 'H', CGT => 'R', CTC => 'L', CCC => 'P', CAC => 'H', CGC => 'R',
        CTA => 'L', CCA => 'P', CAA => 'Q', CGA => 'R', CTG => 'L', CCG => 'P',
        CAG => 'Q', CGG => 'R', ATT => 'I', ACT => 'T', AAT => 'N', AGT => 'S',
        ATC => 'I', ACC => 'T', AAC => 'N', AGC => 'S', ATA => 'I', ACA => 'T',
        AAA => 'K', AGA => 'R', ATG => 'M', ACG => 'T', AAG => 'K', AGG => 'R',
        GTT => 'V', GCT => 'A', GAT => 'D', GGT => 'G', GTC => 'V', GCC => 'A',
        GAC => 'D', GGC => 'G', GTA => 'V', GCA => 'A', GAA => 'E', GGA => 'G',
        GTG => 'V', GCG => 'A', GAG => 'E', GGG => 'G'
};

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
		GENOME_ID => "",
	};

	bless ($self, $class);
	return $self;
}

sub set_genome {
	my ($this, $genome_name) = @_;

	my $dbh = $this->{DBH};
        my $sth = $dbh->prepare("SELECT genome_id FROM genome where genome_name=? OR genome_constant_id=?");
        $sth->execute($genome_name, $genome_name);

        my %feature_sets;
        while(my @res = $sth->fetchrow_array()) {
		$this->{GENOME_ID} = $res[0];
                return $res[0];
        }

	die "Could not find genome_id for $genome_name\n";

}

sub get_feature_location {
	my ($this, $feature_name) = @_;
	my $dbh = $this->{DBH};
        my $sth = $dbh->prepare("SELECT start_pos, stop_pos, direction FROM feature where genome_id=$this->{GENOME_ID} AND (feature_name=? OR feature_symbol=?)");
        $sth->execute($feature_name,$feature_name);
	
        while(my @res = $sth->fetchrow_array()) {
		return @res;
	}
}

sub get_nucleotide_sequence {
	my ($this, $from, $to, $dir) = @_;
	my $len = ($to-$from) + 1;
	my $dbh = $this->{DBH};
        my $sth = $dbh->prepare("SELECT SUBSTRING(genome_sequence, $from, $len) FROM genome where genome_id=$this->{GENOME_ID}");
	$sth->execute();

        while(my $res = $sth->fetchrow_array()) {
		if ($dir eq 'R') {
			return $this->reverse_complement($res);
		}
		else {
			return $res;
		}
	}
}

sub to_fasta {
        my ($this, $seqName, $seq,$len) = @_;

        $len = 80 unless $len;

        my $formatted_seq = ">$seqName\n";
        while (my $chunk = substr($seq, 0, $len, "")) {
                $formatted_seq .= "$chunk\n";
        }

        return $formatted_seq;
}


sub reverse_complement {
        my $this = shift;
        my $dna = shift;
        my $revcomp = reverse($dna);

        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub translate_coding_sequence {
	my $this = shift;
	my $cds = shift;

	my $aa_seq = "";
        while (my $codon = substr($cds, 0, 3, "")) {
#               print "$codon ";
		unless ($TRANS_TABLE->{$codon}) {
			print STDERR "error couldn't find codon $codon\n";
		}
                $aa_seq .= $TRANS_TABLE->{$codon};
        }
	#print "$aa_seq\n";
	$aa_seq =~ s/^./M/; # replace first character with methionine
	$aa_seq =~ s/\*$//; # replace stop codon

	return $aa_seq;
}
sub genome_to_lwgv_js {
	my $this=shift;
	die "you must set genome" unless $this->{GENOME_ID};
	my $gid = $this->{GENOME_ID};
	my $query = "SELECT feature_name, feature_symbol, start_pos, stop_pos FROM feature where genome_id=$gid AND primary_feature='Y'";
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare($query);
	$sth->execute();

	my $js = "var gene_hash = {";
	while (my $row = $sth->fetchrow_hashref()) {
		my $pos = ($row->{start_pos}+$row->{stop_pos})/2;
		my $line = sprintf("\"%s\":%d,", lc($row->{feature_name}), $pos);
		if ($row->{feature_name} ne $row->{feature_symbol}) {
			$line .= sprintf("\"%s\":%d,", lc($row->{feature_symbol}), $pos);
		}
		$js.=$line;
	}

	chop $js;
	$js .= "};\n";

	return $js;
}

sub get_genome_features {
	my $this=shift;
	my $gid = $this->{GENOME_ID};
	my $query = "SELECT feature_name, feature_id FROM feature where genome_id=$gid AND primary_feature='Y'";
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare($query);
	$sth->execute();
	
	my @features;
	while (my @row = $sth->fetchrow_array()) {
		push @features, \@row;
	}

	return \@features;
}

sub get_feature_functions {
	my $this=shift;
	my $feature_id=shift;
	my $gid = $this->{GENOME_ID};
	my $query = "SELECT function_accession, function_name, function_description, function_type, score, eval, pval from function fun, feature_to_function ff where fun.function_id=ff.function_id AND ff.feature_id=? ORDER BY function_type, eval";
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare($query);
	$sth->execute($feature_id);
	my @fun;
	while (my @row = $sth->fetchrow_array()) {
		push @fun, \@row;
	}
	return \@fun;
}
	


sub genome_to_lwgv {
	my $this=shift;
	die "you must set genome" unless $this->{GENOME_ID};
	my $gid = $this->{GENOME_ID};
	my $query = "SELECT feature_name, start_pos, stop_pos, direction FROM feature where genome_id=$gid AND primary_feature='Y'";
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $forward = "track forward_genes addPairs(";
	my $reverse = "track reverse_genes addPairs(";
	my $have_forward;
	my $have_reverse;

#	my $count=0;
	while (my $row = $sth->fetchrow_hashref()) {
#		my $line = "$row->{start_pos}:$row->{stop_pos}:link(\"$row->{feature_name}\",\"#\\\"hello\"),";
		my $line = "$row->{start_pos}:$row->{stop_pos}:link(\"$row->{feature_name}\",\"`javascript:void(0)` onclick=`parent.g('$row->{feature_name}');`\"),";
		if ($row->{direction} eq "F") {
			$forward .= $line;
			$have_forward = 1;
		}
		elsif ($row->{direction} eq "R") {
			$reverse .= $line;
			$have_reverse = 1;
		}
#		last if $count++ > 10;
	}

	chop $forward;
	chop $reverse;
	$forward .= ")\n";
	$reverse .= ")\n";
	my $tracks = "";
	$tracks .= $reverse if $have_reverse;
	$tracks .= $forward if $have_forward;
	return $tracks;
#	return "$reverse$forward";
#	return "$forward$reverse";
}

# for alejandro's phage_omics
sub blast_to_lwgv {
        my $this=shift;
        die "you must set genome" unless $this->{GENOME_ID};
        my $gid = $this->{GENOME_ID};
        my $query = "SELECT feature_name, start_pos, stop_pos, feature_type FROM feature where genome_id=$gid AND feature_type IN ('GutGenBlast', 'VirBlast')";
        my $dbh = $this->{DBH};
        my $sth = $dbh->prepare($query);
        $sth->execute();
        my $forward = "track microbe_blast addPairs(";
        my $reverse = "track virus_blast addPairs(";
        my $have_forward;
        my $have_reverse;

#       my $count=0;
        while (my $row = $sth->fetchrow_hashref()) {
#               my $line = "$row->{start_pos}:$row->{stop_pos}:link(\"$row->{feature_name}\",\"#\\\"hello\"),";
                my $line = "$row->{start_pos}:$row->{stop_pos}:link(\"$row->{feature_name}\",\"`javascript:void(0)` onclick=`parent.g('$row->{feature_name}');`\"),";
                if ($row->{feature_type} eq "GutGenBlast") {
                        $forward .= $line;
                        $have_forward = 1;
                }
                elsif ($row->{feature_type} eq "VirBlast") {
                        $reverse .= $line;
                        $have_reverse = 1;
                }
#               last if $count++ > 10;
        }

        chop $forward;
        chop $reverse;
        $forward .= ")\n";
        $reverse .= ")\n";
        my $tracks = "";
        $tracks .= $reverse if $have_reverse;
        $tracks .= $forward if $have_forward;

        return $tracks;
#       return "$reverse$forward";
#       return "$forward$reverse";
}


sub genome_to_lwgv_contig {
	my $this=shift;
	die "you must set genome" unless $this->{GENOME_ID};
	my $gid = $this->{GENOME_ID};
	my $query = "SELECT contig_description, start_pos, stop_pos FROM contig where genome_id=$gid";
	my $dbh = $this->{DBH};
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $track = "track contigs addPairs(";

	my $count=0;
	while (my $row = $sth->fetchrow_hashref()) {
		$track .= "$row->{start_pos}:$row->{stop_pos},";
		$count++;
	}
	chop $track;

	$track .= ")\n";

	if ($count) {
		return $track;
	}

	return;
}

1;
