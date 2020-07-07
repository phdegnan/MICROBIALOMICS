#!/usr/bin/perl -w

use lib('standard_modules');
use strict;
use microbialomics;

my $DB   = "microbialomics";
my $host = "localhost";
my $user = "root";
my $pass = "ptgdb25";

my $microbe = new microbialomics($DB, $host, $user, $pass);
$microbe->set_genome("Bacteroides thetaiotaomicron VPI-5482");
$microbe->genome_to_lwgv_js();
