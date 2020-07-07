#!/usr/bin/perl 

use strict;
use DBI;
use lib('standard_modules');
use microbialomics;

my $db_name = shift;

die "Usage: perl db_check.pl <db_name>\n" unless $db_name;


#my $db_name = "microbialomics_practice";
my $host = "localhost";
my $user = "root";
my $pass = "ptgdb25";
my $dbh = DBI->connect("DBI:mysql:$db_name;mysql_socket=/var/mysql/mysql.sock",$user,$pass) or die "can't open database: $DBI::errstr\n";
my $microbe = new microbialomics($db_name, $host, $user, $pass);


# optimize tables
optimize_tables();
check_for_best_hit();
optimize_tables();
check_for_unique_feature_symbols();
check_reading_frames();

sub check_for_best_hit {
	print STDERR "adding best hits\n";
	my $table_update = "UPDATE feature_to_function SET best_hit='Y' WHERE feature_id=? AND function_id=?";
	my $sthUpdate = $dbh->prepare($table_update);

	# PFAM
	my $PFAMquery = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='PFAM' GROUP BY feature_id";
	my $sth = $dbh->prepare($PFAMquery);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# TIGRFAM
	my $TIGRFAMquery = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='TIGRFAM' GROUP BY feature_id";
	$sth = $dbh->prepare($TIGRFAMquery);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# COG
	my $COGquery = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='COG (STRING v7.1)' GROUP BY feature_id";
	$sth = $dbh->prepare($COGquery);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# COG category
	my $COGcategory = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='COG category' GROUP BY feature_id";
	$sth = $dbh->prepare($COGcategory);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# KEGG
	my $KEGG = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='KEGG' GROUP BY feature_id";
	$sth = $dbh->prepare($KEGG);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# KEGG category 1
	my $KEGGcat1 = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='KEGG category 1' GROUP BY feature_id";
	$sth = $dbh->prepare($KEGGcat1);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# KEGG category 2
	my $KEGGcat2 = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='KEGG category 2' GROUP BY feature_id";
	$sth = $dbh->prepare($KEGGcat2);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }

	# KEGG pathway
	my $KEGGpath = "select feature_id, ff.function_id, MAX(score) FROM feature_to_function ff, function fu WHERE fu.function_id=ff.function_id AND function_type='KEGG pathway' GROUP BY feature_id";
	$sth = $dbh->prepare($KEGGpath);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) { $sthUpdate->execute($res[0], $res[1]); }
}

sub check_reading_frames {
	my $query = "SELECT genome_id, genome_sequence, genome_constant_id, genome_name FROM genome";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	#my $query1 = "SELECT feature_symbol, feature_type, feature_description, start_pos, stop_pos, direction FROM feature f WHERE primary_feature='Y' AND feature_type='CDS' AND f.genome_id=? ORDER BY start_pos";
#my $sth2 = $dbh->prepare($query1);

	my $query2 = "SELECT feature_symbol, feature_description, start_pos, stop_pos, direction, protein_sequence FROM feature f, protein_seq p WHERE primary_feature='Y' AND feature_type='CDS' AND f.genome_id=? AND p.feature_id=f.feature_id";
	my $sth2 = $dbh->prepare($query2);
	while (my @res = $sth->fetchrow_array()) {
		my $error_genes = 0;


		print STDERR "checking reading frames of $res[2]\t$res[3]\t$res[0]\n";
		$sth2->execute($res[0]);
		my $genome_seq = \$res[1];

		my $num = $sth2->rows;
		print STDERR "genes is $num\n";
		while (my @res2 = $sth2->fetchrow_array()) {
			my $seq = get_sequence($genome_seq, $res2[2], $res2[3], $res2[4]);
			my $aa = $microbe->translate_coding_sequence($seq);
			# ignor first amino acid
			$aa = substr($aa,1);
			my $aa_NCBI = substr($res2[5],1);
			if ($aa ne $aa_NCBI) {
				print "Warning: incorrect translation $res[2]\t$res[3]\t$res2[0]\t$res2[4]\n$seq\n$aa\n$res2[5]\n\n";
				$error_genes++;
			}
		}

		if ($error_genes/$num > 0.25) {
			printf  STDERR "$res[2]\t$res[3]\tmissed: $error_genes (%f)\n", $error_genes/$num;
		}
		printf "$res[2]\t$res[3]\tmissed: $error_genes (%f)\n", $error_genes/$num;

#		return;
	}
#	$sth2->execute($res[0]);
}

sub get_sequence {
	my ($genome_seq, $start, $stop,$dir) = @_;
	my $len = ($stop-$start) + 1;

	my $s = $start-1;
	my $len = $len;
	my $res = substr($$genome_seq, $start-1, $len);

#	print "dir is $dir\n";
	if ($dir eq 'R') {
		$res = reverse_complement($res);
	}

	return $res;
}

sub reverse_complement {
	my $dna = shift;
	my $revcomp = reverse($dna);

	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

sub optimize_tables {
	my $queryA = "show tables";
	my $sth = $dbh->prepare($queryA);
	$sth->execute();
	
	while (my @res = $sth->fetchrow_array()) {
		print STDERR "optimizing table $res[0]\n";
		my $optimize_query = "optimize table $res[0]";
		$dbh->do($optimize_query);
	}

}

sub check_for_unique_feature_symbols {
	# to check if gene names are unique
# give warning for all  non-CDS and die/error for CDS
	my $query = "select genome_name, genome_constant_id, feature_symbol, feature_type, count(feature_id) as feature_count " . 
		    "FROM feature f, genome g WHERE g.genome_id=f.genome_id GROUP BY feature_symbol ORDER BY feature_count DESC LIMIT 20";

	my $sth = $dbh->prepare($query);
	$sth->execute();
	while (my @res = $sth->fetchrow_array()) {
		if ($res[4] > 1) { # non-unique
			if ($res[3] eq 'CDS') {
				#die "Error: multiple gene_symbols not allowed for CDS (@res)\n";
				print "Warning: multiple gene_symbols for: @res\n";
			}
			else {
				print "Warning: multiple gene_symbols for: @res\n";
			}
		}

	}
}


sub to_fasta {
        my ($seqName, $seq, $len) = @_;

        $len = 80 unless $len;

        my $formatted_seq = ">$seqName\n";
        while (my $chunk = substr($seq, 0, $len, "")) {
                $formatted_seq .= "$chunk\n";
        }

        return $formatted_seq;
}

sub reverse_complement {
        my $dna = shift;
        my $revcomp = reverse($dna);

        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

