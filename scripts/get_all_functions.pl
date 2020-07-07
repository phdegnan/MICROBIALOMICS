#!/usr/bin/perl -w

# make a function list for import into excel

use lib('standard_modules');
use strict;
use microbialomics;

my $DB   = "microbialomics_dev";
my $host = "localhost";
my $user = "root";
my $pass = "ptgdb25";

my $microbe = new microbialomics($DB, $host, $user, $pass);
#my $genomes = get_genomes($microbe->{DBH});

my $genome = "Ruminococcus hydrogenotrophicus";
$microbe->set_genome($genome);
my $features = $microbe->get_genome_features();

for my $f (@$features) {
	my $fun = $microbe->get_feature_functions($f->[1]);
	my $text = $f->[0];
	for my $f (@$fun) {
		my $fline = join "|", @$f;
		$text .= "\t$fline";
	}
	$text .= "\n";
	print $text;

}
