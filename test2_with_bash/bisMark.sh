#!/bin/sh -l 

#SBATCH -J bismark.Test
#SBATCH -n 4
#SBATCH --mem-per-cpu=16000

# load the appropriate compute environment
#source activate BS-seq


# check if all required arguments are specified 
if [ -z "$1" ] || [ -z "$2" ]; then
	echo ""
	echo "	Please specify all required arguments;"
	echo "	For Example:"
	echo "		bash bisMark.sh <SE or PE> <WGBS,RRBS or Targeted>"
	echo ""
	exit 1
fi

# initiate workflow for PE data sets 
if [ "$1" = "PE" ] && [ ! -e Reads.list ]  
	then			
			ls -1 *R1* > .Read1.list
			ls -1 *R2* > .Read2.list
			paste -d "," .Read1.list .Read2.list > Reads.list

			if [ "$2" = "WGBS" ]  

				then

					fastqs=( `cat "Reads.list" `) # for Mac + Ubuntu
					#readarray fastqs < Reads.list # for Ubuntu only 

					for i in "${fastqs[@]}" 

						do
							a=`echo "$i" | cut -d ',' -f1`
							b=`echo "$i" | cut -d ',' -f2`
							#echo $i
							bismark --multicore 4 --genome /network/rit/lab/ahmedlab/genomes/hg38.bs.ucsc/ -1 $a -2 $b  # al
								
						done

						for i in *.bam

							do
								deduplicate_bismark --bam $i 		# only recommended for WGBS data-sets 
							
							done 
						
						for i in *.deduplicated.bam

			                do
			                   bismark_methylation_extractor --multicore 4  --bedGraph $i  	
			                
			                done

					

			elif [ "$2" = "RRBS" ] || [ "$2" = "Targeted" ] 

				then 
					fastqs=( `cat "Reads.list" `) # for Mac
			                #readarray fastqs < Reads.list # for Ubuntu

			                for i in "${fastqs[@]}" 

			                        do
			                                a=`echo "$i" | cut -d ',' -f1`
			                                b=`echo "$i" | cut -d ',' -f2`
			                                #echo $i
			                                bismark --multicore 4 --genome /network/rit/lab/ahmedlab/genomes/hg38.bs.ucsc/ -1 $a -2 $b

			                        done

                            for i in *.bam

			                        	do
			                        		bismark_methylation_extractor --multicore 4 --bedGraph $i
			                        	
			                        	done
                    mkdir BAMS COVG
                    mv *.bam BAMS
                    mv *.txt COVG


			else 
				echo " please indicate the correct library type ("WGBS", "Targeted" or "RRBS")"
				exit  

			fi

elif [ "$1" = "SE" ] && [ ! -e Reads.list ] 
	then
	
		ls -1 *R1* > Reads.list

			if [ "$2" = "WGBS" ]  

							then

								fastqs=( `cat "Reads.list" `) # for Mac + Ubuntu
								#readarray fastqs < Reads.list # for Ubuntu only 

								for i in "${fastqs[@]}" 

									do
										a=`echo "$i" | cut -d ',' -f1`
										bismark --multicore 4 --genome /network/rit/lab/ahmedlab/genomes/hg38.bs.ucsc/ --se $a
											
									done

									for i in *.bam

										do
											deduplicate_bismark --bam $i 		# only recommended for WGBS data-sets 
										
										done 
									
									for i in *.deduplicated.bam

						                do
						                   bismark_methylation_extractor --multicore 4 --bedGraph $i  	
						                
						                done

								

			elif [ "$2" = "RRBS" ] || [ "$2" = "Targeted" ] 

				then 
					fastqs=( `cat "Reads.list" `) # for Mac
			                #readarray fastqs < Reads.list # for Ubuntu

			                for i in "${fastqs[@]}" 

			                        do
			                                a=`echo "$i" | cut -d ',' -f1`
			                                bismark --multicore 4 --genome /network/rit/lab/ahmedlab/genomes/hg38.bs.ucsc/ -se $a 

			                        done

			                        for i in *.bam

			                        	do
			                        		bismark_methylation_extractor --multicore 4 --bedGraph $i
			                        	
			                        	done

			else 
				echo " please indicate the correct library type ("WGBS", "Targeted" or "RRBS")"
				exit 1

			fi

elif [ -e Reads.list ]
	then

		echo ""
		echo "	Reads.list file already exists; "
		echo "	Please remove the Reads.list file to re-run this workflow!"
		echo ""
		exit 1
		
else
	echo "please indicate whether a single-end or paired-end experiment (SE or PE)"
	exit 1
fi






#source deactivate


