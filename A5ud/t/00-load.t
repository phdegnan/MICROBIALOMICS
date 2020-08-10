#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'a5ud_pipeline' ) || print "Bail out!\n";
}

diag( "Testing a5ud_pipeline $a5ud_pipeline::VERSION, Perl $], $^X" );
