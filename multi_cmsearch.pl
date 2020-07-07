#!/usr/bin/perl
##
## 	Patrick Degnan
##	multi_cmsearch.pl
##	Use a list of RFAM prefixes to search selected models only
##


use Getopt::Long;
unless(@ARGV){
        print "\nUsage: multi_cmsearch.pl -r rfam_list.txt -g genome.fna\n";
        print "\t-r\tlist file of RFAM *cm sequences to use\n";
        print "\t-g\tgenome sequence(s) in fasta format\n\n";
        print "\t-o\toutput rfam results\n\n";
        exit;
}
GetOptions("out=s"=>\$output,"rfam=s"=>\$rfam,"genome=s"=>\$genome);

unless($output){
	$output=$genome;
	$output=~s/\..+/\.rfam/;
}
open(IN,"<$rfam") or die "cannot open $rfam\n";
while($l=<IN>){
        chomp($l);
        if($l =~/.cm/){$cm="/data/DB/RFAM_11.0/" . $l;}
        else{$cm="/data/DB/RFAM_11.0/" . $l . ".cm";}
        
        $command="cmsearch --cut_tc --acc --tblout table.txt $cm $genome >> rfam_alignments.txt";
        print "$command\n";
        system($command);
        $command="cat table.txt >> $output";
        print "$command\n";
        system($command);	
	

}close(IN);
