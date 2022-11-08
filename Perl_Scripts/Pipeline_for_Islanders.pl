#usr/local/bin/perl


#the program is to obtain 
$input=$ARGV[0];	#IslandersEpilepsyWithHg37Coordinates.vcf
$output1=$ARGV[1];	#MinorOutputFile
$output2=$ARGV[2]; #MajorOutputFile
#$vcf=$ARGV[3]; #VCF

open(I1,"< $input") or die "Unable to read the file\n";
open(O1,"> $output1") or die "Unable to read the file\n";
open(O2,"> $output2") or die "Unable to read the file\n";

#open(V,"> $vcf") or die "Unable to read the file\n";

$rare=0;
$common=0;
$total=0;

print O1 "Chromosome\tPosition\tCommon_Allele\tRareAlleleFrequency\n";

print O2 "Chromosome\tPosition\tRare_Allele\tCommonAlleleFrequency\n";

#To read the 1000G file

while(<I1>){

	$e=$_;

	chop($e);
	
	if($e=~ /\#\#/){
	
		print V $e,"\n";
	
	}
	
	elsif($e=~/\#/){
	
		@header=genomes($e); #total genomes and CScore column array index is returned
		
	#	print $header[2];
		
		#print V $e,"\n";
		
	}
	
	else{
	
		@t=split(" ",$e);
		
		$t[0]=~ s/chr//g;
		
		$chr=int($t[0]);
		
		#print $t[117],"\n";
		
		$each_site=site($e, $header[0],$header[1]);
		
		if($each_site eq "yes"){

			@allele_count=count($e);
			
			#print V $chr,"\t", $t[117], "\t",$t[2],"\t",$t[3], "\t",$t[4] ,"\t",$t[5], "\t",$t[6],  "\t",$t[7],  "\t",$t[8], "\t", $allele_count[3],"\n";
		
			@allele_freq=alleleFreq($allele_count[0],$allele_count[1],$allele_count[2]); #ReferenceHomo, AlternateHomo, Reference Hetero
			
			@major_allele=major($t[3],$allele_freq[0],$t[4],$allele_freq[1]);
			
			if($major_allele[1] ne ""){
			
				print O2 $chr,"\t",$t[117],"\t",$major_allele[0],"\t",$major_allele[1],"\n"; #chromosome, position, common allele, allele frequency
			
				$common++;
			
			}

			
			@rare_allele=minor($t[3],$allele_freq[0],$t[4],$allele_freq[1]);
			
			if($rare_allele[1] ne ""){
			
				print O1 $chr,"\t",$t[117],"\t",$rare_allele[0],"\t",$rare_allele[1],"\n"; #chromosome, position, rare allele, rare allele frequency
		
				$rare++;
			
				
			
			}
			
			$total++
			
	}
	
	else{}
	}
}

print "total number of alleles\t", $total,"\n";

print "rare alleles\t",$rare,"\n";

print "common alleles\t",$common,"\n";

close(I1);

close(O);

close(V);


	

sub count
{
	my $a=$_; #pass line
	
	@t=split(" ",$a);
	
	my $homo11=0;
	my $homo00=0;
	my $hetero=0;
	
	my $str="";
	
		
	for($x=9;$x<=$genomes_till;$x++){
			@u=split(/\:/,$t[$x]);
			
			$str=$str."\t".$u[0];
			
		#	print $str;
			
			$u[0]=~ s/\|//g;
			$u[0]=~ s/\///g;
	
			if($u[0] eq "00"){	
			
			#print $u[0],"\t";
				$homo00++;
			}
			elsif($u[0] eq "11"){
			
				$homo11++;
			}
			elsif(($u[0] eq "01")||($u[0] eq "10")){
			
				$hetero++;
			}	
	}
	
	#print $homo00,"\t",$homo11,"\t",$hetero,"\n";
	
	return $homo00,$homo11,$hetero, $str;
}

sub alleleFreq
{
    # Initializing Variables a, b & c
  	 my  $a = $_[0]; #RefHomo
  	 my  $b = $_[1]; #AltHomo
   	 my  $c = $_[2]; #Hetero
      
    # Performing the operation 
    
 	$total_samples=$a+$b+$c;
 	$freqP=0;
 	$freqQ=0;
 	$freqPq=0;
 	
 	if($a!=0){
 		if($c!=0){
			$freqP=(($a/$total_samples)+(($c/$total_samples)/2)); #Calculates the allele frequency for P(refHomo)
		}
		else{
			$freqP=($a/$total_samples);
		}
	}
	else{
		if($c!=0){
			$freqP=($c/$total_samples)/2;
		}
	}
	if($b!=0){
		if($c!=0){
    		$freqQ=(($b/$total_samples)+(($c/$total_samples)/2)); #Calculates the allele frequency for Q(altHomo)
    	}
    	else{
			$freqQ=($b/$total_samples);
		}
    }
    else{
   		if($c!=0){
     		$freqQ=($c/$total_samples)/2;
     	}
    }
    
    $p=sprintf("%.4f",$freqP);
    $q=sprintf("%.4f",$freqQ);
    $check=$freqP+$freqQ;
   
  # 	print  $a,"\t", $b,"\t", $c, "\t", $freqP, "\t", $freqQ, "\t",$p,"\t",$q,"\t", $check,"\n";
    
    # Function to return the value back to the function
    
    return $p, $q;
    
    #$freqP,$freqQ;
    
}


sub major
{

	my $ref_allele=$_[0]; #Ref Allele
	my $freqP=$_[1]; #Ref Allele Frequency
	my $alt_allele=$_[2]; #Alt Allele
	my $freqQ=$_[3]; #Alt Allele Frequency
	
	if(($freqP>=0.05)and($freqP!=0)){

	    	return $ref_allele,$freqP,"P";  	

    }
   	 elsif(($freqQ>=0.05)and($freqQ!=0)){
  		 return $alt_allele,$freqQ,"P";
  	}
  	else{}
}

sub minor
{

	my $ref_allele=$_[0]; #Ref Allele
	my $freqP=$_[1]; #Ref Allele Frequency
	my $alt_allele=$_[2]; #Alt Allele
	my $freqQ=$_[3]; #Alt Allele Frequency
	
	if(($freqP<=0.01)and ($freqP!=0)){
			#print O $ref_allele, "\t", $freqP,"\n";
	    	return $ref_allele,$freqP,"\tP";  	

    }
   	 elsif(($freqQ<=0.01)and ($freqQ!=0)){
   		 
   		# print O $alt_allele,"\t",$freqQ ,"\tQ\n";#"\tQ is the rare Allele\n";
  		 return $alt_allele,$freqQ,"Q" ;
  	}
  	else{}
}

sub site{

    # Initializing Variables a, b & c
  	 my  $a = $_[0]; #line
  	 my  $b = $_[1]; #total genomes
   	 my  $c = $_[2]; #Cscore
   	 
   	 @t=split(" ",$a);
   	 
   	# print $t[$c], "\n";
   	
   #	if($t[$c]>=20){
   	
    	if(($t[3] eq "A" || $t[3] eq "T" || $t[3] eq "G" || $t[3] eq "C") && ($t[4] eq "A" || $t[4] eq "T" || $t[4] eq "G" || $t[4] eq "C")){
   	 	 
   	 		return "yes";
   	 	
   	 	 }
   	 	 else{
   	 	 	return "no";
   	 	 }
   #	 }

}

sub genomes{

	my $a=$_; #pass line
	
	 $genomes_till=0;
	
	my @t=split(" ",$a);

	my $ii=0;
	
	my $coordinate=0;

	while($t[$ii] ne ""){

		$header[$ii]=$t[$ii];
		
		if(($t[$ii] eq "ANC" )||($t[$ii] eq "Anc")){
		
			$anc=$ii;
			
			$genomes_till=$anc-1;
			
			$c_score=$anc+1;
			
			#print $ii,"\t", "ANC\n";
		
		}
		else{
		
			if($t[$ii] eq "CScore" ){
			
				$c_score=$ii;
				
				$genomes_till=$c_score-1;
			
			}
			elsif($t[$ii] eq "HG37Coordinate"){
			
				$genomes_till=$ii-1;

			}
			
			else{}
		}
		if($genomes_till==0){
		
			$genomes_till=$ii-1;
			
			$c_score=0;
			
		}
		
		$ii++;
	
	}
	
	return $genomes_till,$c_score;
	
}






  