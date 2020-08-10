a5ud_pipeline
=============
This package was an update to original software developed by 
Andrew Tritt and Aaron Darling for the original A5 pipeline.
This was further modified by Nick Youngblut but is no longer available on GitHub
https://github.com/nyoungb2/a5ud_pipeline


The pipeline requires multiple dependencies which are in the bin/ folder.
Add the bin/ folder to your PATH to install.

Check any additional file paths in use in the a5ud_pipeline.pl script, make
them match you system.


SUPPORT AND DOCUMENTATION


• Install using standard install to global or user directory for Perl Modules

$ perl ./Build.PL --prefix /home/user/PM
$ ./Build
$ ./Build test
$ ./Build install


• If using local directory make sure it is added to the PERL Path:

export PERL5LIB=/home/usr/PM/share/perl/5.14.2/


• IDBA_UD is not distributed with this package (only idba). 

Install independently https://github.com/loneknightpy/idba



• After installing, you can find documentation using the -h option


$ a5ud_pipeline.pl -h

/data/software/PM/bin/a5ud_pipeline.pl: No files given.
Usage:
    Method 1) Run using paired-end reads in separate files:
        a5ud_pipeline.pl [options] -read1 read1.fastq -read2 read2.fasta -out my_assembly

    Method 2) Run using paired-end reads interleaved in the same file:
        a5ud_pipeline.pl [options] -inter read1-2.fasta -out my_assembly

    Method 3) Run using a library file:
        a5ud_pipeline.pl [options] -lib library_file

  Options:
    -begin <int>
        Step in the pipeline to begin at (Steps 1-5). [1]

    -end <int>
        Step in the pipeline to end at (Steps 1-5). [5]

    -threads <int>
        Number of threads to use with IDBA-UD. [1]

    -min <int>
        Minimum contig size to report with IDBA-UD. [1]

    -preprocessed <bool>
        Use if starting after Step 2.

    -h This help message



Making a Library file is easy: 


[LIB]
p1=Sample_R1.fastq
p2=Sample_R2.fastq
ins=500

**Note: If you write your own file make sure there are no extra/blank lines in the file or else you will get an error.**

Now you can run a5ud pipeline on the data:

$ a5ud_pipeline.pl -threads 1 -lib Sample_lib.txt -out Sample_run1 -min 500
[a5] Begin: 16:03 on 1-8-2016 
[a5] Found the following libraries:
     raw1:
      id=raw1
      p1=Sample_R1.fastq
      p2=Sample_R2.fastq
      ins=500
[a5] Found 1 libraries
[a5] Starting pipeline at step 1
[a5] Cleaning reads with SGA
[a5] Cleaning reads with SGA
...
[a5] Finish:  16:18 1-8-2016
[a5] Took 0:0:0 to completed

Output folders & files:

Sample_run1.s1			      		// Quality trimming report
Sample_run1.s2			      		// Assembly report
Sample_run1.s3 		            	// Initial scaffolding report
Sample_run1.s4    					// Final scaffold quality check
Sample_run1.ec.fastq.gz             // Error corrected reads
Sample_run1.contigs.fasta           // Contigs
Sample_run1.crude.scaffolds.fasta   // Crude scaffolds:  Not checked for misassemblies
Sample_run1.broken.scaffolds.fasta  // Broken scaffolds: Checked for misassemblies, but not rescaffolded
Sample_run1.final.scaffolds.fasta   // Final scaffolds:  Checked for misassemblies, and rescaffolded




LICENSE AND COPYRIGHT

Copyright (C) 2014-2020 Patrick Degnan

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See <http://dev.perl.org/licenses/> for more information.
