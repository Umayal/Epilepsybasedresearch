#usr/local/bin/perl

$input= $ARGV[0];
$output= $ARGV[1];

@chr=qw(chr1	chr2	chr2	chr3	chr5	chr5	chr5	chr5	chr8	chr9	chr9	chr9	chr9	chr10	chr10	chr11	chr11	chr12	chr14	chr15	chr16	chr16	chr16	chr17	chr17	chr19	chr20),

@pos1= qw(42925353	165194993	165984641	10992186	45254948	88717117	126531200	161847063	132120859	127579370	128191655	128552558	135702185	87862638	93757840	790475	72189558	13437942	28766787	92900189	9753404	29811382	78099400	6684719	63682336	13206442	63400208);

@pos2= qw(42958893	165392310	166149214	11039247	45696498	88904257	126595362	161899981	132481095	127696027	128255248	128633662	135795508	87971930	93806272	798281	72196323	13982002	28770277	93027996	10182928	29815892	79212667	6713377	63741986	13633025	63472677);

$length=scalar @chr;
 
open(I,"< $input")|| print "Unable to read the file\n";

open(O,"> $output")|| print "Unable to write the file\n";

while(<I>){

		$e=$_;

		chop($e);
		
		@b=split(' ',$e);
		
		if($e !~ /\#/){
		
			if((length($b[3]) eq 1) and (length($b[3]) eq 1)) {
			
			for($i=0;$i<=$length;$i++){
			
				if($chr[$i] eq $b[0]){
			
					if(($b[1]>=$pos1[$i]) && ($b[1]<=$pos2[$i])) {
			
						print O $e, "\n";
			
					}
			
				}
			
			}

			}
		

				
		}
		
		else{
		
			print O $e,"\n";
			
		}
		
}
			
close(I);
		
close(O);