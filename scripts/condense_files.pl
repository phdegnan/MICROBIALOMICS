#!/usr/bin/perl
##
## 	Patrick Degnan
##	condense_files.pl
##	Copy raw annotation files from Microbialomics subdirectories to new folders
##	Typically for transfer to different server and upload to MySQL
##



unless(@ARGV){die "usage: condense_files.pl <file of folder names>\n";}

@ext1=(".ptt",".contig",".fna",".faa",".gbk",".rnt");
@ext2=(".CELLO",".KEGG",".PFAM",".TIGRFAM",".stringCOG");
open(IN,"<$ARGV[0]") or die "Cannot open $ARGV[0]\n";

##NC vs. NZ don't expect *contig

while($l=<IN>){
#starting with list of files
	chomp($l);
	print "$l\n";
	mkdir($l);
	$file="";
	foreach $e (@ext1){		
		if(-e "/home/pd275/proc_genomes/$l/$l$e"){# if file exists
			$command="cp /home/pd275/proc_genomes/$l/$l$e $l/";
			system($command);
		}else{print "File: /home/pd275/proc_genomes/$l/$l$e does not exist!\n";}
	}
	foreach $e (@ext2){
		($file,@rest)=glob("/home/pd275/proc_genomes/$l/ann_*/$l$e");
		if(defined($file)){# if file exists
			$command="cp /home/pd275/proc_genomes/$l/ann_*/$l$e $l/";
			system($command);
		}else{print "File: /home/pd275/proc_genomes/$l/ann_*/$l$e does not exist!\n";}

	}
	print "Done\n";
}
