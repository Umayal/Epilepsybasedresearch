#usr/local/bin/perl

$rare=0;
$common=0;
$total=0;

$minoutput="/RDS/Q1233/Uma/Uma_Work/coordinates/Minor_1000G.txt";	#Minor OutputFile

$majoutput="/RDS/Q1233/Uma/Uma_Work/coordinates/Major_1000G.txt";	#Major OutputFile

$output="/RDS/Q1233/Uma/Uma_Work/coordinates/All_Chromosomes_GBR.vcf"; #VCF File with high deleterious Variants

open(O1,"> $minoutput") or die "Unable to read the file\n";

open(O2,"> $majoutput") or die "Unable to read the file\n";

open(S,"> $output") or die "Unable to read the file\n";

for($l=1;$l<=22;$l++){

	$input="/RDS/Q1233/Uma/Humans/Delt/GBR/GBR_Pops_1000G_".$l."_with_Anc_and_Phred.vcf";	#vcf File

	print $input,"\n";

	open(I1,"< $input") or die "Unable to read the file\n";

$flag=0;
	#To read the 1000G file

	while(<I1>){

		$e=$_;

		chop($e);
	
		if($e=~ /\#\#/){
		
			if($flag==0){
			
			print S $e,"\n";
			
			}
		
		}
	
		elsif($e=~/\#/){
	
			@header=genomes($e); #total genomes and CScore column array index is returned
		
			if($flag==0){
			
			print S $e,"\n";
			
			$flag=1;
			
			}
		
		}
	
		else{
	
			@t=split(" ",$e);
			
			
		
			$each_site=site($e, $header[0],$header[1]);
		
			if($each_site eq "yes"){
			
				print S $e,"\n";

				@allele_count=count($e);
				
				#print $allele_count[0],"\t",$allele_count[1],"\t",$allele_count[2],"\n";
		
				@allele_freq=alleleFreq($allele_count[0],$allele_count[1],$allele_count[2]); #ReferenceHomo, AlternateHomo, Reference Hetero
			
				@major_allele=major($t[3],$allele_freq[0],$t[4],$allele_freq[1]);
			
				if($major_allele[1] ne ""){
			
			
					print O2 $t[0],"\t",$t[1],"\t",$major_allele[0],"\t",$major_allele[1],"\n"; #chromosome, position, common allele, allele frequency
			
					$common++;
			
				}
				else{}

				@rare_allele=minor($t[3],$allele_freq[0],$t[4],$allele_freq[1]);
			
				if($rare_allele[1] ne ""){
			
					
					print O1 $t[0],"\t",$t[1],"\t", $rare_allele[0],"\t", $rare_allele[1],"\n";
			
					$rare++;
			
				}
				else{}
				
				
			
				$total++;	
			}
	
			else{}
		}
		
}

close(I1);

}

close(O1);
close(O2);

print "total number of alleles\t", $total,"\n";

print "rare alleles\t",$rare,"\n";

print "common alleles\t",$common,"\n";



sub count
{
	my $a=$_; #pass line
	
	@t=split(" ",$a);
	
	my $homo11=0;
	my $homo00=0;
	my $hetero=0;
	my $r=9;
	
	
	while($t[$r] ne ""){
		
	
	
	#print $t[$r],"\t";
	
			#@u=split(/\:/,$t[$x]);
			
			$t[$r]=~ s/\|//g;
			$t[$r]=~ s/\///g;
	
			if($t[$r] eq "00"){	
			
			#print $u[0],"\t";
				$homo00++;
			}
			elsif($t[$r] eq "11"){
			
				$homo11++;
			}
			elsif(($t[$r] eq "01")||($t[$r] eq "10")){
			
				$hetero++;
			}	
			
			#print $homo00,"\t";
			$r++;
	}
	
	return $homo00,$homo11,$hetero;
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
   
 	 # print  $a,"\t", $b,"\t", $c, "\t", $freqP, "\t", $freqQ, "\t",$p,"\t",$q,"\t", $check,"\n";
    
    # Function to return the value back to the function
    
    return $p, $q;

    
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
	    	return $ref_allele,$freqP,"\tP";  	
    }
   	 elsif(($freqQ<=0.01)and ($freqQ!=0)){
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
   	
   	if($t[$c]>=15){
   	
    	if(($t[3] eq "A" || $t[3] eq "T" || $t[3] eq "G" || $t[3] eq "C") && ($t[4] eq "A" || $t[4] eq "T" || $t[4] eq "G" || $t[4] eq "C")){

   	 		return "yes";
   	 	
   	 	 }
   	 	 else{
   	 	 	return "no";
   	 	 }
   	 }

}

sub genomes{

	my $a=$_; #pass line
	
	my $genomes_till=0;
	
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






  