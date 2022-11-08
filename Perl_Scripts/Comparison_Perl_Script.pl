$input1=$ARGV[0];	#Islander Common Allele file
$input2=$ARGV[1];	#Europeans Rare Allele File
$output=$ARGV[2];	#Output File


open(I1,"< $input1")|| print "Unable to read the file\n";


while(<I1>){

	$e=$_;

	chop($e);

	@b=split(' ',$e);
	
	$arr[$b[0]][$b[1]]=$e; #storing the value in an array
	
}

close(I1);


open(I2,"< $input2")|| print "Unable to read the file\n";

open(O,"> $output")|| print "Unable to wrte the file\n";

$i=0;

$u=0;

$highfrequency=0;


while(<I2>){

	$e=$_;
	
	chomp($e);
	
	@r=split(' ',$e);	

	#print "hi\n";
	
		if($arr[$r[0]][$r[1]] ne "") { #if the European genome coordinate matches the  Islanders genomes
	
				$i++;
	
	  			@p=split(' ',$arr[$r[0]][$r[1]]); 
	  	
	  			if($p[2] eq $r[2]){ # if alleles from both File matches
	  			
	  				$u++;
	  				
	  				#$w=$p[3]%$r[3];
	  				
	  				#if($w>=2){
	  				
	  					#$highfrequency++;
	  				
	  				#}
	  				
	  				#else{}
	  	
	  				print O $arr[$r[0]][$r[1]],"\t\t\t", $e,"\n";
			
				}
		}
	
		else{}
		
}

$i--;

print "Total number of Deleterious variants associated with Epilepsy genes: ",$i,"\n";

print $u," are rare in Europeans genomes from 1000 genomes but common in the Islander VCF file\n";

#print $task1,"\n";

#print $highfrequency,"\n";

close(I2);

close(O);

