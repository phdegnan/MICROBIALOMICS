#!/usr/bin/perl
##
## 	Patrick Degnan
##	pre_process_genomes.pl
##	Decompress necessary files and that all required files are present
##	Make necessary files if missing.
##
unless(defined(@ARGV)){die "usage: pre_process_genomes.pl <list of folder names>\n";}

open(IN,"<$ARGV[0]") or die "Cannot open $ARGV[0]\n";
while($f=<IN>){
	chomp($f);
	print "[$f]\n";
	if($f =~ /^NZ\_/){
		@extensions=(".ptt",".faa",".rnt",".fna",".gbk");
## gunzip *
		print "Unzipping and compressing raw files\n";
		@files=glob("$f/*z");	# getting list of files from local directory that match extension
		foreach $compressed (@files){
			$uncompressed=$compressed;
			if($uncompressed=~/\.tgz/){$uncompressed=~s/\.tgz/\.tar/;}
			else{$uncompressed=~s/\.gz//;}
			$command="gunzip < $compressed > $uncompressed";
			##command accounts for impromperly compressed files which are occasionaly problematic e.g.,
			## gzip: NZ_ACAC00000000.contig.gbk.tgz: unexpected end of file
			#print "$command\n";
			system($command);
		}
		#$unzip ="gunzip $f/*gz";
		#print "$unzip\n";
		#system($unzip);
	
## tar -xvmf each

		@files=glob("$f/*tar");	# getting list of files from local directory that match extension
		foreach $compressed (@files){
			$command="tar -xmf $compressed -C $f";
			#print "$command\n";
			system($command);
		}

	}else{
		@extensions=(".ptt",".faa",".rnt",".fna");
	}
	$mainprefix=$f;
	
## check if *ptt for each *faa
## if not use gbk2ptt.pl
	print "Checking for missing ptt files\n";
	@files=glob("$f/*faa");
	foreach $fasta (@files){
		$ptt=$fasta;
		$ptt=~s/\.faa/\.ptt/;
		unless(-e $ptt){
			$gbk=$fasta;
			$gbk=~s/\.faa/\.gbk/;
			$command="gbk2ptt.pl $gbk";
			print "$command\n";
			system($command);
		}
	
	}
	if(length($files[1]) > 34){
		$files[1]=~/(N\S+\/N\S{8})\S+\.faa/;
		$contigprefix=$1;
	}elsif($files[1]=~/^NC/){###
		$files[1]=~/(N\S+\/N\S{7})\S+\.faa/;
		$contigprefix=$1;
	}else{
		$files[1]=~/(N\S+\/N\S{6})\S+\.faa/;
		$contigprefix=$1;
	}
	
## cat *ptt and each *faa using parent prefix
## && remove uncompress files 
	print "Concatenating and removing individual files\n";
	foreach $e (@extensions){
		$command="cat $contigprefix*$e > $mainprefix/$mainprefix$e"; 
		print "$command\n";
		system($command);
		$command="rm $contigprefix*$e";
		print "$command\n";
		system($command);

	}

## gzip *
	print "Re-compressing original files\n";
	$gzip ="gzip $f/*tar";
	#print "$gzip\n";
	system($gzip);

}##end of foreach
