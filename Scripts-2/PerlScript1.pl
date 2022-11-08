$input1=$ARGV[0];	#genelist
$input2=$ARGV[1];	#Annotated_File
$output=$ARGV[2];	#OutputFile


open(I1,"< $input1")|| print "Unable to read the file\n";

#To read the Epilepsy gene list first

$str="";

while($f=<I1>){

	chop($f);

	@bb =split "\n",$f;
	
	$str=$str.$bb[0].";"; 

}

close(I1);


#********************************************************************************
# To compare the read Epilepsy gene list with the UCSC complete annotated gene list
#********************************************************************************

open(I2,"< $input2")|| print "Unable to read the file\n";

open(O,"> $output")|| print "Unable to write the file\n";

while(<I2>){

	$e=$_;
	
	chomp($e);
	
	@b=split(' ',$e);	
	
	if($e !~ /\#/){
	
		if($b[2]=="gene"){
		
			@t=split(/\;/,$b[8]);
			
			@gene=split(/\=/,$t[3]);
			
			@k=split(/\;/,$str);
			
			$i=0;
			
			while($k[$i] ne ""){
			
				if($gene[1] eq $k[$i]){	
					
					print O $gene[1],"\t",$e,"\n";
					
					$x[$i]=$gene[1]."\t".$e."\n";
									
				}
				
				$i++;
				
			}
		}
	}
}


close(I2);



