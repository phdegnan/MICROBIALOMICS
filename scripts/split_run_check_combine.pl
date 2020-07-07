#!/usr/bin/perl -w

#    Script: split_run_check_combine.pl
#    ___________________________________________________________________________
#
#    Version 1
#
#    Copyright (C) 2008-2009 Brian Muegge and Jeremiah Faith
#	Additional changes to code made by Patrick Degnan 2011-2012
#	Original code commented out, changes noted w/ intials
#
#    http://hamlet.wustl.edu/microbialomics_dev/docs_db/
#
#    About: License
#
#       Licensed under the GNU General Public License
#
#        This program is free software; you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation; either version 2 of the License, or
#        (at your option) any later version.
#
#        This program is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with this program; if not, visit http://www.gnu.org/licenses/gpl.txt
#        or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
#        Boston, MA  02111-1307  USA.
#
# Topic: Getting started
# 
#
# (start code) 
#  perl split_run_check_combine.pl -i <input (in fasta or scarf format)> -p <program and parameters (put in single quotes; include INCLUDE_INFILE)> -n <num_lines_per_job>
#
#  # Example (will blastp proteins.faa against the string database using 40 sequences per batch job)
#  perl split_run_check_combine.pl -i proteins.faa -p '/srv/cgs/local/i386/ncbi-blast/latest/bin/blastall -p blastp -d /srv/cgs/data/STRING_v8.0/protein.sequences.v8.0.fa -i INCLUDE_INFILE -m 8 -e 10e-10' -n 40
# (end)
#
# This script aims to provide a robust/generalizable framework for starting and checking jobs on the 
# cluster. It splits a FASTA or scarf file into many smaller files, runs a user-defined program on each small file, checks to make sure that each of the programs
# run on the smaller files completed successfully, and if they did complete successfully, it combines all of the outputs into a single file.
# 
# You specify where to place the infile name for each subfile in your program using the INCLUDE_INFILE.  You should try to set -n so that jobs on the cluster take around 15 minutes.
#
# Topic: How it works
# This script splits a scarf file or a fasta file into many smaller files.  It run a particular program 
# creating a successful-completion file *.complete when finished and a *.start file before running (in case the job doesn't finish)
# run a job that waits for the first job to finish and checks to make sure everything is ok (otherwise it should re-run the failed jobs, but this isn't implemented yet)
# if all of the jobs are ok, it should combine all of the outputs and remove the temp files from the split
#
# IMPORTANT: in your program program parameters put the word INCLUDE_INFILE whereever you want the infile to be placed
# this script assumes your program prints the output to STDOUT
#

use strict;
#use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case);
use FindBin qw($Bin);

my $infile;
my $num_lines_per_job;
my $program;
my $mode = 'split';
my $code = "c" . int(rand() * 100000); # used to track things
my $file_suffix;
my $infile_suffix;
my $this_program = "$Bin/split_run_check_combine.pl";
my $outfile;
my $memory;
my $use_long_queue;
my $node_arch;
my $JOB_ID;
my $processors=1; 	##phd
my $time="24:00:00";	##phd
my $CWD=`pwd`; 		##phd
my $queue="mdi";	##phd
chomp($CWD);

GetOptions ( #"s" => \$start_job,
            "infile=s"   => \$infile,
            "num_lines_per_job=i"   => \$num_lines_per_job,      # integer
            "code=s"   => \$code,      # integer
	    "params=s"     => \$program,
	    "suffix=s" => \$file_suffix, # used after combining the files to 
	    "Suffix_infile=s" => \$infile_suffix, 
	    "mode=s" => \$mode,
	    "Memory=i" => \$memory,
	    "LongQueue" => \$use_long_queue,
	    "arch_node=s" => \$node_arch,						#Mod to define the node arch
	    "outfile=s" => \$outfile,
	    "jobid=i" => \$JOB_ID,
	    "Proc=i" => \$processors, 	##phd
	    "time=s"=> \$time, 		##phd
	    "queue=s"=> \$queue, 	##phd
            );  # flag


#print "have memory $memory long queue $use_long_queue\n";

#$code = $code .'_'. $file_suffix if $file_suffix;


if ($mode eq 'split') { # in split mode, split the file, start the job and start a job to check the jobs 
	usage() unless $infile && $program && $processors;
#	die "usage: perl $this_program -i <input (in fasta or scarf format)> " .
#	"-p <program and parameters (put in single quotes; include INCLUDE_INFILE)> -n <num_lines_per_job> -s [file suffix for joined results]\n" unless $infile && $program && $num_lines_per_job;

	print STDERR "code is $code\n";
	##mkdir "SGE";
	
	# split the big sequence file into a bunch of little ones
	my ($infile_suffix, @split_files) = split_jobs($infile, $processors, $code);  #replace number of lines per job to number of processors
	my @job_ids = run_jobs(\@split_files, $program, $code, $memory, $use_long_queue, $node_arch);

	my $num_jobs = scalar @split_files;
	append_job_check(\@job_ids, $code, $num_jobs, $infile, $file_suffix, $infile_suffix);
}
elsif ($mode eq 'recheck') { # in split mode, split the file, start the job and start a job to check the jobs 
	# blah
	my ($num_jobs, $infile, $file_suffix, $infile_suffix, $lines_jobs, $jobs_to_rerun) = get_recheck_jobs($code, 1);

	# remove all of the complete and start jobs so they can be properly recreated in "run" mode
	system("rm $code.*.complete $code.*.start");

#	print "have num jobs $num_jobs\nhave rerun jobs:\n$jobs_to_rerun\n";

	# rerun the failed jobs and check for completion
	my @job_ids = run_jobs("", "", $code, $memory, $use_long_queue, $node_arch, $jobs_to_rerun);
	append_job_check(\@job_ids, $code, $num_jobs, $infile, $file_suffix, $infile_suffix);
}
elsif ($mode eq 'run') { # run mode is called by the jobs in split mode; it simply runs each job of the user's program; creating a *.start file and a *.complete file which are used by the check mode 
	my $hostname = my $host = `hostname`;
#	print "in run on host $hostname\n";

	# create a start file with the stuff to be executed in case something fails and the job needs to be restarted

	my $sub_id = parse_sub_id($program, $code);

	# recommended by Brian (thanks to the heads up from Alejandro)
	# this allows you to use just the program name (and NOT the full path) on the different nodes
	# and it imports the correct environment for the node (i.e. for 32 bit it looks in the 32-bit path)
	# The line is commented since apparently is no longer needed on the CGS cluster. Each node should be able to identify the adequate path
	# eval('cgsenv -f -b');

	system("echo '$program' > $code.$sub_id.start");
	# run the program with the parameters passed by the split mode
	my $status = system("$program");

	# check return status
	if ( $? == -1 ) {  # program never executed
		my $error_msg  = sprintf "Error: command ($program) failed to execute: $!\n";
		system("echo '$error_msg' > $code.$sub_id.failed");
	}
	elsif ($status != 0) { # program executed but failed for some reason
		my $exit_value = $? >> 8;
		my $error_msg = sprintf "command failed (did not exit with status zero); exited with value %d", $exit_value;
		system("echo '$error_msg' > $code.$sub_id.failed");
	}
	else { # if the program completed successfully write a complete file
		my $exit_value = $? >> 8;
		printf "command exited with value %d", $exit_value;
		system("touch $code.$sub_id.complete");
	}
}
elsif ($mode eq 'check') { # in check mode, we just make sure everyfile has a *.start and *.complete and concatenate the outputs
	my $num_files = $num_lines_per_job;
	my ($failed, $files_to_delete) = check_files($code, $num_files, $JOB_ID);
	print "to delete @$files_to_delete\n";
	if (scalar @$failed) {
		create_recheck_files($failed, $code, $files_to_delete, $num_files, $infile, $file_suffix, $infile_suffix);
		unlink(@$files_to_delete);
		die "Error missing some files: @$failed\n";
	}
	else {
		combine_files($code, $num_files, $infile, $outfile, $file_suffix, $files_to_delete, $infile_suffix)
	}
}

sub get_recheck_jobs { 
	my $code = shift;
	my $get_rerun_jobs = shift;

	my $recheck_file = "$code.recheck";

	open (IN, $recheck_file) or die "can't open $recheck_file: $!\n";
	my $num_lines = <IN>;
	chomp $num_lines;
	my $infile = <IN>;
	chomp $infile;
	my $file_suffix = <IN>;
	chomp $file_suffix;
	my $infile_suffix = <IN>;
	chomp $infile_suffix;

	my $text = "";	
	my @lines;
	while (<IN>) {
		chomp;
		s/\s+//g;
		$text .= $_;
		@lines = split /,/, $text;
	}
	close IN;

	my @lines_jobs;
	for my $l (@lines) {
		my ($line, $subjob) = split /:/, $l;
		push @lines_jobs, [$line,$subjob];
	}

	# make things easy to check for using a hash
	my %need_to_recheck;
	for my $l (@lines_jobs) {
		$need_to_recheck{$l->[0]} = 1;
	}


	my $recheck_jobs = "";
	my $old_submission = "$code.sh";


	if ($get_rerun_jobs) {
		# go into the old submission file and grab out the jobs that failed
		open (IN, $old_submission) or die "can't open $old_submission: $!\n";
		my $count = 0;
		while (<IN>) {
			if ($need_to_recheck{$count}) {
				$recheck_jobs .= $_;
			}
			$count++;
		}
		close IN;
	}
#	print "lines is @lines from $text from file $recheck_file jobs are:\n$recheck_jobs\n";

	return $num_lines, $infile, $file_suffix, $infile_suffix, \@lines_jobs, $recheck_jobs;
}

sub create_recheck_files { 
	my ($failed, $code, $files_to_delete, $num_files, $infile, $file_suffix, $infile_suffix) = @_;


	my $recheck_file = "$code.recheck";
	open (OUT, ">$recheck_file") or die "can't open $recheck_file: $!\n";
	my $to_check = join ",", @$failed;
	
	# first line is the number of files in the original query
	print OUT "$num_files\n";
	# second line is the name of the original infile
	print OUT "$infile\n";
	# third line is the file_suffix to be used by the outfile
	print OUT "$file_suffix\n";
	# fourth line is the file_suffix to be used by the outfile
	print OUT "$infile_suffix\n";
	# second line is id of each subfile for the resubmission
	print OUT "$to_check\n";
}

# concatenate all of the results in order
#sub combine_files { 
#	my ($code, $num_files, $infile, $outfile, $files_to_delete) = @_;
#
#	my $outfiles = "";
#	for my $i (0 .. $num_files-1) {
#		push @$files_to_delete, "$infile.$code.$i";
##		system("rm $infile.$code.$i"); # remove the sequence subset file
#		$outfiles .= " $infile.$code.$i.out.gz";
#		my $gz_outfile = "$infile.$code.$i.out.gz";
#		system("zcat $gz_outfile >> $outfile"); # concatenate all of the files into outfile
#		push @$files_to_delete, $gz_outfile;
##		system("rm $gz_outfile"); # remove the sequence subset file
#	}
#
	#my $recheck_file = "$code.recheck";
	# second line is id of each subfile for the resubmission
	#print OUT "$to_check\n";
#}#

# concatenate all of the results in order
sub combine_files { 
	my ($code, $num_files, $infile, $outfile, $file_suffix, $files_to_delete, $infile_suffix) = @_;
	#	combine_files($code, $num_files, $infile, $outfile, $file_suffix, $files_to_delete)

	my $outfiles = "";
	for my $i (0 .. $num_files-1) {
#		system("rm $infile.$code.$i"); # remove the sequence subset file
		my $gz_outfile;
		if ($infile_suffix) {
			push @$files_to_delete, "$infile.$code.$i.$infile_suffix";
			$outfiles .= " $infile.$code.$i.$infile_suffix.out.gz";
			$gz_outfile = "$infile.$code.$i.$infile_suffix.out.gz";
		}
		else {
			push @$files_to_delete, "$infile.$code.$i";
			$outfiles .= " $infile.$code.$i.out.gz";
			$gz_outfile = "$infile.$code.$i.out.gz";
		}

		if ($i == 0) {
			system("zcat $gz_outfile > $outfile"); # concatenate all of the files into outfile (first file creates a new file)
		}
		else {
			system("zcat $gz_outfile >> $outfile"); # concatenate all of the files into outfile
		}
		push @$files_to_delete, $gz_outfile;
#		system("rm $gz_outfile"); # remove the sequence subset file
	}

	my $recheck_file = "$code.recheck";
	if (-e $recheck_file) { # don't need the recheck file anymore
		push @$files_to_delete, $recheck_file;
	}

	push @$files_to_delete, "$code.check.sh";
	push @$files_to_delete, "$code.sh";
	push @$files_to_delete, "$code.check_status";

	print "deleting @$files_to_delete\n";

	# now that everything is joined, delete all of the files
	unlink(@$files_to_delete);
	#foreach my $f (@$files_to_delete){ unlink($f) or die "Can't delete $f : $!";} 
	###
	#system("mv *$code* SQ_Files_$JOB_ID.omega-rocks.local"); 		##phd
	#system("mv SQ_Files_$JOB_ID.omega-rocks.local SQ_Files_$JOB_ID");	##phd
	#system("tar -cvmf SQ_Files_$JOB_ID.tar SQ_Files_$JOB_ID.bulldogi-rocks.hpc.yale.edu");
	#system("gzip SQ_Files_$JOB_ID.tar");
	#system("rm -r SQ_Files_$JOB_ID.bulldogi-rocks.hpc.yale.edu");
}

# make sure that all of the files are there otherwise return the missing ones so they can be rerun if desired
sub check_files {
	my ($code, $num, $JOB_ID) = @_;

	my $out = "$code.check_status";
	open (OUT, ">$out") or die "can't open $out: $!\n";
	my @failed_jobs;
	my $num_failed = 0;
	my @files_to_delete;
	my $sge_err;
	my $infile;


	# CHECK FROM HERE, NEED TO MAKE SURE THIS WORKS

	my $num_to_check = $num-1;

	# not going through the normal way, need to only scan the recheck files
	my ($num_jobs, $lines_jobs);
	if (-e "$code.recheck") {
	        ($num_jobs, $infile, $file_suffix, $infile_suffix, $lines_jobs) = get_recheck_jobs($code);

		print "using recheck\n";
		$num_to_check = $#$lines_jobs;	
	}

	# CHECK FOR start and complete files
	for my $i (0 .. $num_to_check) {
		my $line = $i;
		my $subjob = $i;

		# if working from a recheck, need to use a different coordinant system
		if ($lines_jobs) {
#			$line  = $lines_jobs->[$i][0];
			$subjob = $lines_jobs->[$i][1];
		}

		
		my ($start, $complete) = (0,0);##for louise, symbiont
		my $start_file = "$code.$subjob.start";
		if (-e $start_file) { 
			#print OUT "found\n";  # printing success cluttered everything up
			$start=1;
			#print "$start_file\n";
		}else{ 
			print OUT "$start_file not found\n"; 
		}

		my $complete_file = "$code.$subjob.complete";
		if (-e $complete_file) { 
	#		print OUT "found\n"; # printing success cluttered everything up
			$complete=1;
			#print "$complete_file\n";
		}else{ 
			print OUT "$complete_file not found\n"; 
		}
		


		#print "[$start][$complete][$sge_err_OK]\n";
		if ($start && $complete ) {  # remove the check files if we have all three
			push @files_to_delete, $start_file;
			push @files_to_delete, $complete_file;
			print "adding $start_file $complete\n";
		}else{
			push @failed_jobs, "$line:$subjob";  # save the jobs that failed
		}
	}

	if (scalar @failed_jobs) {
		print OUT "\n\nSummary: failed jobs @failed_jobs\n";
	}else{
		print OUT "\n\nSummary: no failed jobs\n";
	}
	close OUT;

	return \@failed_jobs, \@files_to_delete;
}



# start a job to check that all of the previous jobs completed; this job starts after the first job-array completes
sub append_job_check {
	my ($job_ids, $code, $num_jobs, $infile, $file_suffix, $infile_suffix) = @_;
	#print "In append_job_check\n";

	my $outfile = $infile . ".$code.combined";
	if ($file_suffix) {
		$outfile = $infile . ".$file_suffix";
	}
	
	# have the output checked
	my $program = "$this_program -c $code -n $num_jobs -m check -i $infile -o $outfile -j $$job_ids[0]";
	$program .= " -s $file_suffix" if $file_suffix;
	$program .= " -S $infile_suffix" if $infile_suffix;

	system("echo '$program' > $code.check.sh");
	#my $qsub = 'qsub -e $HOME/.sge -o $HOME/.sge -V -cwd -r y -hold_jid' .  " $job_id";
#	my $qsub = 'qsub -e $PWD/SGE -o $PWD/SGE -cwd -l arch=lx24-amd64 -r y -hold_jid' .  " $job_id";
#	my $qsub = 'qsub -e $PWD/SGE -o $PWD/SGE -cwd -r y -hold_jid' .  " $job_id";

## difference between PBS Torque and PBS for Sun microsystems
##	my $qsub = 'qsub -e $PWD/SGE -o $PWD/SGE -d $PWD -r y -W depend=afterany:' .  "$job_id";##
	#my $qsub = 'qsub -q $queue -e $PWD -o $PWD -d $PWD -r y -W depend=afterany:' .  "$job_id";####phd
	#print "End append_job_check\n[$qsub $code.check.sh]\n";##
	#system("$qsub $code.check.sh");
	
	my $command="/home/pdegnan/Scripts/wait.pl $code.check.sh";
	foreach my $p (@$job_ids){
		$command.=" $p";
	}
	print "$command\n";
	system("$command 2>> temp1.log&");
	
}

sub parse_sub_id {
	my $text = shift;
	my $code = shift;

	if ($text =~ /\.$code\.(\d+)/) {
		return $1;
	}

	return;
}


#$outfile = "$scarf_file.$suffix";
#print SH_OUT "$program -q $outfile -D $db $params | gzip > $outfile.out.gz\n";

#close SH_OUT;

# start the job on the cluster
#if ($start_job) {
#	system("nq $sh_out | qsub");
#}

# set up a batch file to run all of the jobs with nq;  rather than run the user's job directly, 
# we're going to all this script to run the job, so we can create a *.start and *.complete file
# to report on the progress/failure of each job
sub run_jobs {
	my ($files, $program, $code, $memory, $use_long_queue, $node_arch, $recheck_jobs) = @_;

	# make a batch file
	my $batch_file = $code . ".sh";
#	print "opening $batch_file: $!\n";
#	open (OUT, ">$batch_file") or die "can't open: $batch_file: $!\n";
	my @PROCESSES=();
	my $pid="";
	my $command="";
	
	if ($recheck_jobs) {
		print "opening $batch_file: $!\n";
		open (OUT, ">$batch_file") or die "can't open: $batch_file: $!\n";
		print OUT "$recheck_jobs";
		close OUT;
	## Submit job & recover pid
		$command="sh $batch_file >> temp1.log 2> pid &"; ##STDERR holds the pid
		print "$command\n";
		system($command);
		sleep(3);
		open(PID,"<pid") or die "Cannot open pid\n";
		$pid=<PID>; chomp($pid);
		close(PID);
		`rm pid`;
		print "[$pid]\n";
		push(@PROCESSES, $pid);
	
	
	
	
	}else {
		for my $f (@$files) {
			print "opening $batch_file: $!\n";
			open (OUT, ">$batch_file") or die "can't open: $batch_file: $!\n";

			my $command = $program;	
			$command =~ s/INCLUDE_INFILE/$f/g;
			if ($command =~ /INCLUDE_OUTFILE/){							#Modified to add the possibility of an output file on the command line of the program
				$command =~ s/INCLUDE_OUTFILE/$f.out/g;
				print OUT "echo \$\$ \>\&2\ndate\ncd $CWD\n$this_program -c $code -m run -p '$command'\ngzip $f.out\ndate\n";##phd
			}else{	
				print OUT "echo \$\$ \>\&2\ndate\ncd $CWD\n$this_program -c $code -m run -p '$command | gzip > $f.out.gz'\ndate\n";##phd
			}
			close OUT;
			## Submit job & recover pid
			$command="sh $batch_file >> temp1.log 2> pid &"; ##STDERR holds the pid
			print "$command\n";
			system($command);
			sleep(2);
			open(PID,"<pid") or die "Cannot open pid\n";
			$pid=<PID>; chomp($pid);
			close(PID);
			`rm pid`;
			print "[$pid]\n";
			push(@PROCESSES, $pid);
			
			
		}
	}

	#close OUT;

#	my $params = '-e $PWD/SGE -o $PWD/SGE ';

	

#	print "using params $params\n";
#	system("sqPBS.py -q $queue -n $processors -w $time -N $code $batch_file | qsub");##phd

#	my $id;
#	while (!($id = get_job_id($code))) {
#		sleep(10);
#	}
	print "Returning from run_jobs [$PROCESSES[0]]\n";
	return @PROCESSES;
}

sub get_job_id {
	my $code = shift;

	print STDERR "getting job id\n";
#	sleep(10);

	my @jobs = `qstat`;

#	my $job_id;
#	my @stuff;
	my $jobid;
	for my $j (@jobs) {
		$j =~ s/^\s+//;
		$j =~ s/\s+$//;
		if ($j =~ /$code/) {
			my ($jid, @stuff) = split /\s+/, $j;	
			$jobid = $jid;
			$jobid=~s/\.omega\-rocks//;##
#			print "have jobid $jobid from @jobs\n";
			print STDERR "found job id $jobid\n";
			return $jobid;
		}
	}
	print STDERR "found job id $jobid\n";

	return;	
}

# split the file into smaller files with $num_lines_per_job in each file
sub split_jobs {
	my ($infile, $num_jobs, $code) = @_; ## numb lines per job switched to number of processors allowed

	my $file_type = get_file_type($infile);
	my $numlines = get_num_seqs($file_type, $infile);
	my $num_lines_per_job=0;
	if(($numlines%$num_jobs)==0){
		$num_lines_per_job=int($numlines / $num_jobs);
	}else{$num_lines_per_job=1+int($numlines / ($num_jobs));}
	printf "splitting $numlines sequence $infile into $num_jobs jobs w/ $num_lines_per_job sequences\n";
	my @files; # of format infile: outfile
	my $file_suffix;


	if ($file_type eq 'scarf') {
		@files = split_scarf($infile, $num_lines_per_job, $code);
	}
	elsif ($file_type eq 'fasta') {
		($file_suffix, @files) = split_fasta($infile, $num_lines_per_job, $code);
	}

	return $file_suffix, @files;
}


# split a scarf file into smaller files with $num_lines_per_job in each file
sub split_scarf {
	my ($scarf_file, $num_lines_per_job, $code) = @_;
	my @files;
	my $suffix = 0;
	my $line_num=0;

	open (IN, $scarf_file) or die "can't open $scarf_file: $!\n";

	my $prefix = "$scarf_file.$code";
	my $outfile = "$prefix.$suffix";
	# open the first file
	open (OUT, ">$outfile") or die "can't open $scarf_file: $!\n";
	push @files, $outfile;
	while(<IN>) {
		if (0 == ($line_num++ % $num_lines_per_job)) {  # reached the desired number for a file
			unless ($line_num == 1){
				close OUT;	
				$suffix++;
				$outfile = "$prefix.$suffix";
				open (OUT, ">$outfile") or die "can't open $outfile: $!\n"; # open a new file
				push @files, $outfile; # store the file name
			}
		}
		print OUT $_;
	}
	close OUT;

	return @files;
}

sub get_file_suffix {
	my $file = shift;

	my @pieces = split /\./, $file;
	if (scalar @pieces > 1) {
		return $pieces[$#pieces];
	}

	return;
}

# split a fasta file into smaller files with $num_lines_per_job in each file
sub split_fasta {
	my ($fasta_file, $num_lines_per_job, $code) = @_;
	my @files;
	my $suffix = 0;
	my $line_num=0;
	my $file_suffix = get_file_suffix($fasta_file);

	open (IN, $fasta_file) or die "can't open $fasta_file: $!\n";

	my $prefix = "$fasta_file.$code";
	my $outfile = "$prefix.$suffix";
	$outfile .= ".$file_suffix" if ($file_suffix); # retain the user's file suffix
	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
	push @files, $outfile; # store the file name
	while(<IN>) {
		#print "[$line_num][$num_lines_per_job]\n";
		if (/^>/ && 0 == ($line_num++ % $num_lines_per_job) ) {
			unless($line_num ==1){
				close OUT;
				$suffix++;
				$outfile = "$prefix.$suffix";					
				$outfile .= ".$file_suffix" if ($file_suffix); # retain the user's file suffix
				open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
				push @files, $outfile; # store the file name
			}
		}
        	print OUT $_;
        }
	close OUT;

	return $file_suffix, @files;
}



# determine if file is fasta or scarf (i.e. solexa)
sub get_file_type {
	my $in=shift;	
	open (IN, $in) or die "can't open $in: $!\n";

	while (<IN>) {
		next if /^\s*$/; # skip blank lines

		if (/^>/) {
			close IN;
			return 'fasta';
		}
		elsif (/:\d+:\d+:\d+:/) {
			close IN;
			return 'scarf';
		}
		else {
			close IN;
			die "Could not recognize file format for: $in\n";
		}
	}

}


sub get_num_seqs {
	my ($type, $infile) = @_;

	if ($type eq 'fasta') {
		return get_num_seqs_fasta($infile);
	}
	elsif ($type eq 'scarf') {
		return get_num_seqs_scarf($infile);
	}
}

sub get_num_seqs_scarf {
	my $in=shift;

	open (IN, $in) or die "can't open $in: $!\n";
	
	my $line_count=0;
	while (<IN>) { $line_count++; }
	close IN;

	return $line_count;
}

sub get_num_seqs_fasta {
	my $in=shift;

	open (IN, $in) or die "can't open $in: $!\n";
	
	my $line_count=0;
	while (<IN>) { $line_count++ if /^>/; }
	close IN;

	return $line_count;
}


sub usage {
	print "\n\tUsage: perl split_run_check_combine.pl\n".
	"\t\t-i <input (in fasta or scarf format)>\n" .
	"\t\t-p <program and parameters (put in single quotes; include INCLUDE_INFILE)>\n".
#	"\t\t-n <num_lines_per_job>\n".
	"\t\t-s [file suffix for joined results]\n".
#	"\t\t-M [minimum memory requirements in GB : NOT TESTED ON OMEGA]\n".
#	"\t\t-L [use long queue : NOT TESTED ON OMEGA]\n".
#	"\t\t-a [architecture of the node needed. i.e.: -a arch=lx24-amd64 : NOT TESTED ON OMEGA]\n".
	"\t\t-m [recheck (to recheck, only provide -c <code> and optionally -M and -L)]\n".
	"\t\t-P [Number of nodes to request from sqPBS.py default = 1 node 8 processors]\n\n".
#	"\t\t-t [Time for job to run in format 'HH:MM:SS', default = '24:00:00']\n".
#	"\t\t-q [Which OMEGA queue, default = mdi]\n\n".
	"\tIf the program requires an output file just write INCLUDE_OUTFILE where the outfile should go. The joined output will be still Input.suffix\n\n".
	"\tExamples:\n\n".
	"\t# split test.faa info 5 sequences per file and blast against proteinDB, joining the results into a file called test.faa.blastp\n".
	"\tsplit_run_check_combine.pl -P 1 -i test.faa -p '/home/mdi/goodman/shared/scripts/blast-2.2.25+/bin/blastp -db /home/mdi/goodman/pd275/DB/micro_v3_aa -query INCLUDE_INFILE -out INCLUDE_OUTFILE -evalue 1e-4 -outfmt 6' -n 5 -s M36.br\n\n".
	"\t# try to recover failed jobs from c4789\n".
	"\tsplit_run_check_combine.pl -c c4789 -m recheck \n".
	"\n";
	exit(1);
}
