#!/usr/bin/perl
##
## 	Patrick Degnan
##	gbk2ptt.pl
##	Convert gbk or gbff file into ptt files
##


unless(defined($ARGV[0])){die "Enter GenBank formatted file\n";}
$file=$ARGV[0];

$/="     gene            ";
@types=("CDS","rRNA","RNA","tRNA","tmRNA","misc_RNA","ncRNA");
$out=$file;
$out=~s/\.gbk/\.ptt/;
open(OUT,">$out") or die "Cannot open $out\n";
	
open(IN,"<$file") or die "Cannot open $file\n";
$first=<IN>;
if($first=~ /^LOCUS\s+(\S+)\s+(\d+)\s+bp/){
	print OUT "- 1..$2\n";
	$cds=`grep -c "     CDS             " $file`;
	$cds=~/^(\d+)/;
	print OUT "$1 proteins\n";
}else{
	print OUT "- 1..X\nX proteins\n";
}
print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";

while($line=<IN>){

	if($line=~ /\/locus\_tag\=\"(.*?)\"/){
 		$name=$1;
 		#$entire= ">gnl|NAM|$name ";
 		#print "$name\t";
 	}else{$name=""; }#print "\t";
 	
 	if($line=~ /gene\=\"(.*)\"/){
		$gene=$1;
		#print "$gene\t";
		#$entire= $entire . "$gene ";
	}else{$gene="-";}#print "\t";
	
	if($line=~ /product\=\"(.*?)\"\s+.+?/s){
		$product=$1;
		$product=~ s/\n\s+/ /g;
		#$entire= $entire . "$product ";
		#print "$product\t";
	}elsif($line=~ /note\=\"(.*?)\"\s+.+?/s){
		$product=$1;
		$product=~ s/\n\s+/ /g;
		#print "$note\t";
	}else{$product="";}#print "\t";
	
	if($line=~/\/db\_xref\=\"GI\:(\d+)\"/){
		$gi=$1;
	}else{$gi="";}

	if($line=~/\/db\_xref\=\"COG\:(COG\d+)\"/){
		$cog=$1;
	}else{$cog="-";}	

	if($line=~/EcoGene\:(EG\d+)\"/){
		$eg=$1;
	}else{$eg="";}
		
	if($line=~ /\/pseudo\n/){
		if($line=~ /(\d+)\.\.(\d+)/){
			$start=$1;
			$stop=$2;
			$strand="+";
		}elsif($line=~ /complement\((\d+)\.\.(\d+)\)/){
			$start=$1;
			$stop=$2;
			$strand="-";	
		}
		$type="pseudogene";	
	}else{		
		foreach $cat (@types){
			if($line=~ /\s+$cat\s+<*(\d+)\.\.>*(\d+)/){
				$start=$1;
				$stop=$2;
				$strand="+";
				$type=$cat;
				last;
			}elsif($line=~ /\s+$cat\s+complement\(<*(\d+)\.\.>*(\d+)\)/){
				$start=$1;
				$stop=$2;
				$strand="-";
				$type=$cat;	
				last;
			}
		}
	}

	if($type eq "CDS"){
		$aa=((($stop-$start+1)/3) - 1);
		$rounded = sprintf("%.0f", $aa);
		print OUT "$start..$stop\t$strand\t$rounded\t$gi\t$gene\t$name\t-\t$cog\t$product\n";
	}
	
 
}