#!/usr/bin/perl
##
## 	Patrick Degnan
##	post_process_genomes.pl
##	Check for presence of required files and compress unnecessary files.
##

unless(defined(@ARGV)){die "usage: post_process_genomes.pl <list of folder names>\n";}

open(LIST,">>processed_genomes_list.txt") or die "Cannot open processed_genomes_list.txt\n";

#open(IN,"<$ARGV[0]") or die "Cannot open $ARGV[0]\n";
#while($f=<IN>){
foreach $f (@ARGV){
	chomp($f);
	print "[$f]\n";
	if($f =~ /^NZ\_/){
		@extensions=(".ptt",".faa",".rnt",".fna",".gbk");
## gunzip *
		print "Unzipping and compressing raw files\n";
		@files=glob("$f/*.gz");	# getting list of files from local directory that match extension
		foreach $compressed (@files){
			$uncompressed=$compressed;
			$uncompressed=~s/\.gz//;
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
			#$command="tar -xmf $compressed";
			#print "$command\n";
			system($command);
		}
		if( -d "$f/$f"){
			$command="mv $f/$f/* $f/";
			system($command);
			$command="rm -r $f/$f";
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
		$files[1]=~/(N\S+\/N\S{5})\S+\.faa/;
		$contigprefix=$1;
	}
	#print "CONTIG PREFIX:[$contigprefix]\n";
	unless(defined($contigprefix)){die "CONTIG PREFIX:[$contigprefix]\n";}

## Open *contig file to get accessions and coordinate lengths 
	print "Getting contig coordinates and accessions\n";
	$contig_file="$f/$f.contig";
	open(IN,"<$contig_file") or die "Cannot open $contig_file\n";
	while($l=<IN>){
		@cols=split(/\t/,$l);
		
		if(defined($SIZE{$cols[0]})){die "\nWarning 2 contigs have same id:[$cols[0]]\n";}
		else{$SIZE{$cols[0]}=$cols[1];}
		$total=$cols[2];
	}
	close(IN);
	
	$definition=`head -1 $f/$f.faa`;
	$definition=~/\[(.+?)\]/;
	$name=$1;
	print "Organism = [$name]\n";
## Re-make rnt and ptt files and copy original concatenated versions to *_cat.ptt
	print "Remaking ptt and rnt files\n";
	if(-e "$f/$f.ptt"){# if file exists
		$abbv=&REMAKE($contigprefix,"ptt",$f,$total,$name,%SIZE);
	}
	if(-e "$f/$f.rnt"){# if file exists
		$blank=&REMAKE($contigprefix,"rnt",$f,$total,$name,%SIZE);
	}
	
## tar new files overwriting old ones
	print "Re-tar with added files\n";
	foreach $e (@extensions){
		($compressed,@rest)=glob("$f/*$e*.tar");
		$command="tar -cvmf $compressed $contigprefix*$e"; 
		print "$command\n";
		system($command);
		$command="rm $contigprefix*$e";
		print "$command\n";
		system($command);

	}

## gzip new files
	print "Re-compressing new tar files\n";
	$gzip ="gzip -f $f/*tar";
	print "$gzip\n";
	system($gzip);
	print LIST "$mainprefix\t$abbv\t\t$name\n";
}##end of foreach
close(LIST);
##subroutine


sub REMAKE {
	#&REMAKE($contigprefix,"ptt",$f,$total,%SIZE);
	my ($star,$ext,$folder,$genome,$name,%CONTIGS)=@_;
	my $prefix="";my $number="";
	my @files=glob("$star*$ext");
	my $first="Y"; my $origin; my $out; my $ptt_file; my $l; my $count=0; my $head="";
	open(OUT,">$folder/temp") or die "Cannot open $folder/temp\n";
	foreach $ptt_file (@files){
		open(IN,"<$ptt_file") or die "Cannot open $ptt_file\n";
		#$head="";$out="";
		$ptt_file=~/\/(.+)\.$ext/;
		$origin=$CONTIGS{$1};
		while ($l=<IN>){
		
			if($l !~ /^\d+/){	
				if($first eq "Y"){
					if($l=~/\d+\.\.\d+/){
						#$l=~s/\d+\.\.\d+/(1\.\.$genome)/;
						$range="1..$genome";
					}
					$head="$name, whole genome shotgun - $range\n";
				}
			}elsif($l=~ /(^\d+)\.\.(\d+)/){
				if($first eq "Y"){$first="N";}
				
				$old_start=$1;
				$old_stop=$2;
			
				$new_start=$origin + $old_start - 1;
				$new_stop=$origin + $old_stop - 1;
				
				my @cols=split(/\t/,$l);
				if($cols[5]=~/\_/){
					($prefix,$number)=split(/\_/,$cols[5]);
				}elsif($cols[5]=~/(^.+?)\d+$/){
					$prefix=$1;
				}

				$l=~s/^\d+\.\.\d+/$new_start\.\.$new_stop/;
				$out=$out . $l;
				$count++;
	
			}
		}##end while
		#print OUT $out;
	}##end foreach
	print OUT $head;
	if($ext eq "ptt"){print OUT "$count proteins\n";}
	else{print OUT "$count RNAs\n";}
	print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
	print OUT $out;
	
	close(OUT);
	my $command="mv $folder/$folder.$ext $folder/$folder\_cat.$ext";
	print "$command\n";
	system($command);
	$command="mv $folder/temp $folder/$folder.$ext";
	print "$command\n";
	system($command);

	return($prefix);

}
