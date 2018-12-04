###### Copy all scripts and star count files i.e. files ending with *.ReadsPerGene.out.tab into a new directory and execute the following steps in order:

  ### Step 1 : Execute the processCounts.sh script using the following syntax:
           $ bash processCounts.sh <1> or <2> ( 1 for first strand, 2 for second strand)
           
  ### Step 2: Create a phenoData.csv file. Template is available in the Utilites folder.
          Make sure col 2 has file names generated in Step 1
  
  ### Step 3: Execute the countMatrix.sh script using the following syntax:
          $ bash countMatrix.sh
          
  ### Step 4: Use the following command to run the final DESeq2.R script:
          $ RScript DESeq2.R <Experiment.Name> <Numerator> <Denominator> # (OS X)
          $ Rscript DESeq2.R <Experiment.Name> <Numerator> <Denominator> # (Linux)
          
          
  
           
           
