use 5.006;
use strict;
use warnings FATAL => 'all';
use Module::Build;

opendir IN, "./bin/" or die $!;
my @binfiles = readdir IN;
my @binscripts = grep(/.pl$/, @binfiles);
@binscripts = map{$_ = "./bin/" . $_} @binscripts;


my $builder = Module::Build->new(
    module_name         => 'a5ud_pipeline',
    license             => 'perl',
    dist_author         => q{Nick Youngblut <nyoungb2@illinois.edu>},
    dist_abstract      => 'a5_pipeline with idba_ud',
    dist_version_from   => 'lib/a5ud_pipeline.pm',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => 0,
    },
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        #'ABC'              => 1.6,
        #'Foo::Bar::Module' => 5.0401,
    },
    add_to_cleanup     => [ 'a5ud_pipeline-*' ],
    create_makefile_pl => 'traditional',
    scripts => \@binscripts
);

$builder->create_build_script();

