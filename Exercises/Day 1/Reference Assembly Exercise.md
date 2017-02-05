#dDocent Reference Assembly Exercise
Designed by Jon Puritz

**NOTE: You can download the Ref.Ex file from the depository and run this tutorial from the command line**

#GOALS
1.	To set up test data for our exercise
2.	To demultiplex samples with process_radtags and rename samples 
3.	To use the methods of dDocent (via rainbow) to assemble reference contigs
4.	To learn how optimize a de novo reference assembly
5.	To learn how to utilize pyRAD to assemble loci

#Tutorial
*dDocent must be properly installed for this tutorial to work*

Welcome to the first exercise of B@G 2017!!

Let's get started.  First let's create a working directory for yourself for the Day 1 workshop
```bash
mkdir D1W
```
Let's change into that directory and load the dDocent2.2.7 module
```
cd D1W
module load dDocent/v2.2.7
```
Now let's get the test data I created for the course.


No surprises here from our simulated data, butI highly recommend familiarizing yourself with grep, awk, and regular expressions to help evaluate de novo references.

#Bonus Section

Here, I am going to let you in on an experimental script I have been using to help optimize reference assemblies.

```bash
curl -L -O https://raw.githubusercontent.com/jpuritz/WinterSchool.2016/master/RefMapOpt.sh
```

This script assembles references across cutoff values and then maps 20 random samples and evaluates mappings to the reference, along with number of contigs and coverage.  
It takes a long time to run, but here's a sample command and output

```bash
#RefMapOpt.sh 4 8 4 8 0.9 64 PE
```

This would loop across cutoffs of 4-8 using a similarity of 90% for clustering, parellized across 64 processors, using PE assembly technique.

The output is stored in a file called `mapping.results`

```bash
curl -L -o mapping.results https://www.dropbox.com/s/x7p7j1xn1hjltzv/mapping.results?dl=0
cat mapping.results
```

```
Cov		Non0Cov	Contigs	MeanContigsMapped	K1	K2	SUM Mapped	SUM Properly	Mean Mapped	Mean Properly	MisMatched
37.3382	39.6684	1000	942.25				4	4	747510		747342			37375.5		37367.1			0
37.4003	39.7343	1000	942.25				4	5	748753		748546			37437.7		37427.3			0
37.4625	39.7919	1000	942.45				4	6	749999		749874			37499.9		37493.7			0
37.4967	39.8282	1000	942.45				4	7	750685		750541			37534.2		37527.1			0
37.486	39.8169	1000	942.45				4	8	750469		750205			37523.4		37510.2			0
37.3517	39.6785	1000	942.35				5	4	747780		747612			37389		37380.6			0
37.4147	39.7454	1000	942.35				5	5	749042		748835			37452.1		37441.8			0
37.4701	39.7999	1000	942.45				5	6	750151		750009			37507.6		37500.4			0
37.4852	39.8161	1000	942.45				5	7	750453		750210			37522.7		37510.5			0
37.4551	39.7824	999		941.55				5	8	749102		748837			37455.1		37441.8			0
37.3561	39.6833	1000	942.35				6	4	747870		747731			37393.5		37386.6			0
37.453	39.7776	1000	942.55				6	5	749809		749734			37490.4		37486.7			0
37.4923	39.8193	1000	942.55				6	6	750595		750376			37529.8		37518.8			0
37.4784	39.8089	1000	942.45				6	7	750318		750075			37515.9		37503.8			0
37.4437	39.766	999		941.65				6	8	748874		748616			37443.7		37430.8			0
37.4013	39.7312	1000	942.35				7	4	748774		748698			37438.7		37434.9			0
37.4592	39.7907	1000	942.4				7	5	749934		749835			37496.7		37491.8			0
37.4682	39.7981	1000	942.45				7	6	750114		749897			37505.7		37494.8			0
37.4239	39.7468	1000	942.55				7	7	749227		748993			37461.3		37449.7			0
37.417	39.736	998		940.75				7	8	747591		747320			37379.6		37366			0
37.4413	39.761	1000	942.65				8	4	749575		749499			37478.8		37474.9			0
37.4492	39.7843	1000	942.3				8	5	749733		749562			37486.7		37478.1			0
37.4441	39.7711	998		940.6				8	6	748133		747888			37406.7		37394.4			0
37.4274	39.7517	997		939.7				8	7	747052		746779			37352.6		37338.9			0
37.5014	39.8269	989		932.25				8	8	742528		742279			37126.4		37113.9			0
```

I have added extra tabs for readability.  The output contains the average coverage per contig, the average coverage per contig not counting zero coverage contigs, the number of contigs, the mean number of contigs mapped, the two cutoff values used, the sum of all mapped reads, the sum of all properly mapped reads, the mean number of mapped reads, the mean number of properly mapped reads, and the number of reads that are mapped to mismatching contigs.
Here, we are looking to values that maximize properly mapped reads, the mean number of contigs mapped, and the coverage.  In this example, it's easy.  Values 4,7 produce the highes number of properly mapped reads, coverage, and contigs.  
Real data will involve a judgement call.  Again, I haven't finished vetting this script, so use at your own risk.

# pyRAD assembly tutorial

Now, let's take a look at another way to assemble RAD data from the software package pyRAD.  Please note that many of these steps have been altered from Deren Eaton's tutorial
See http://nbviewer.ipython.org/gist/dereneaton/dc6241083c912519064e/tutorial_pairddRAD_3.0.ipynb for more details
First let's make a new directory and move into it

```bash 
mkdir pyrad
cd pyrad
module load pyRAD/3.0.6
```
Next, let's make a symoblic link to the original fastq.gz data files and barcodes

```bash 
ln -s ../SimRAD.barcodes .
ln -s ../SimRAD_R1.fastq.gz SimRAD_R1_.fastq.gz
ln -s ../SimRAD_R2.fastq.gz SimRAD_R2_.fastq.gz
```
THE FORMATTING HERE IS CRITICAL.  The pyRAD package requires the files to have the _R1_ and _R2_ in the names.  
In case you weren't aware, this makes a virtual link (like one on a desktop) to a file and saves disk space by not recopying the files

Now, let's create a parameters file

```bash 
pyrad -n
```	
	
This creates a file (params.txt) that we can edit to adjust the settings of pyRAD

We will need to edit a few values.  You can do this in a text editor like nano or emacs, but for this exercise it is easier just to use sed
First let's change the restriction sites to match our data

```bash 
sed -i '/## 6. /c\AATT,CCG                 ## 6. cutsites... ' ./params.txt
```	
Next, let's change the number of processors to use in parallel to 3

```bash 
sed -i '/## 7. /c\8\t                 ## 7. N processors... ' ./params.txt
```
Change the datatype to paired ddRAD

```bash 
sed -i '/## 11. /c\pairddrad                 ## 11. datatype... ' ./params.txt
```
Now, we are ready to proceed with the pyRAD pipeline.  First step is demultiplexing files

```bash 
pyrad -p params.txt -s 1
```
You should see:
```
  ------------------------------------------------------------
   pyRAD : RADseq for phylogenetics & introgression analyses
  ------------------------------------------------------------


  step 1: sorting reads by barcode
	 .
```
This now created demultiplexed files with the proper pyRAD naming convention.  There are stats about the demultiplexing in the ./stats directory

The next step is quality filtering.  This is basic filtering, removing any reads with Illumina adapters in them, and replacing low quality bases with Ns.


```bash 
pyrad -p params.txt -s 2
```
You should see
```
  ------------------------------------------------------------
   pyRAD : RADseq for phylogenetics & introgression analyses
  ------------------------------------------------------------


  step 2: quality filtering 
  ........................................
```
Stats about the filtering can be found in the ./stats directory.  For the simulated data set, no reads are filtered.
Unlike the previous method, pyRAD first clusters reads together within individuals for assembly
```bash 
pyrad -p params.txt -s 3
```
The output should look like:
```
   ------------------------------------------------------------
    pyRAD : RADseq for phylogenetics & introgression analyses
   ------------------------------------------------------------


	de-replicating files for clustering...

	step 3: within-sample clustering of 40 samples at 
        '.88' similarity. Running 8 parallel jobs
 	with up to 6 threads per job. If needed, 
	adjust to avoid CPU and MEM limits

	sample PopB_10 finished, 847 loci
	sample PopB_08 finished, 863 loci
	sample PopA_20 finished, 863 loci
	sample PopA_08 finished, 865 loci
	sample PopB_12 finished, 862 loci
```
This will continue through all 40 samples
When the clustering step completes we can examine the results by looking at the file s3.clusters.txt in the /stats directory
```bash 
head -6 ./stats/s3.clusters.txt
```
```
taxa	total	dpt.me	dpt.sd	d>5.tot	d>5.me	d>5.sd	badpairs
PopA_01	869		19.824	9.084	823		20.747	8.424	78
PopA_02	867		19.722	8.822	824		20.576	8.193	95
PopA_03	870		20.322	9.59	828		21.176	9.026	88
PopA_04	899	1	9.339	9.081	860		20.076	8.582	75
```
This output shows us the total number of clusters for each individual, along with some information about mean depth and standard deviation of depth.

It also shows us the number of bad pairs, or mismatched 1st and 2nd reads.  In this example, we are seeing a large number ~10% of mismatched forward and reverse reads.

Considering this simulated data does NOT have any paralogs in it, there should be a very low percentage of mismatched reads.
Let's examine some good and bad clusters
The clusters are in the `./clust88 directory`.  Let's look at a bad one first.
```bash 
zcat ./clust.88/PopA_01.badpairs.gz | head -12
```
```
>PopA_01_9392_pair;size=9;
AATTTGTGGGTTTCTCCTTAAAAGATTACCAAATTCTAGTATCAATCATCCTCCTCCCAATGCATGGAGACTGGCAACACCGTGCAGTAGCCT---nnnnTCTCGGCGGATTTGTTTACCCGCGAAGTCGTAA-CTA--CCACCACTCGACCCAACCGGTCCTAGATGACTGCTGTCATACAAT-GTCGTACCGATGA-AGA---CGG
>PopA_01_9402_pair;size=6;+
AATTTGTGGGTTTCTCCT--AAAGATTACCAAATTCTAGTATCAATCATCCTCCTCCCAATGCATGGAGA-TGGCAACACCGTGCGGTAGCCTAGAnnnn-------CGATTTGTTTACCC-CGAAGTCGTAAGCTGACCAACCACTCTACCCAACCGGTCCTAGATGACTGGTGTCATACAATCGTCGTACCGATGATAGACTGCGG
>PopA_01_9401_pair;size=1;+
AATTTGTGGGTTTCTCCT--AAAGATTACCAAATTCTAGTATCAATCATCCTCCTCCCAATGCATGGAGA-TGGCAACACCGTGCGGTAGCCTAGAnnnn-------CGATTTGTTTACCC-CGAAGTCATAAGCTGACCAACCACTCTACCCAACCGGTCCTAGATGACTGGTGTCATACAATCGTCGTACCGATGATAGACTGCGG
>PopA_01_9409_pair;size=1;+
AATTTGTGGGTTTCTCCT--AAAGATTACTAAATTCTAGTATCAATCATCCTCCTCCCAATGCATGGAGA-TGGCAACACCGTGCGGTAGCCTAGAnnnn-------CGATTTGTTTACCC-CGAAGTCGTAAGCTGACCAACCACTCTACCCAACCGGTCCTAGATGACTGGTGTCATACAATCGTCGTACCGATGATAGACTGCGG
>PopA_01_9408_pair;size=1;+
AATTTGTGGGTTTCTCCT--AAAGATTACCAAATTCTAGTATCAATCATCCTCCTCCCAATGCATGGAGA-TGGCCACACCGTGCGGTAGCCTAGAnnnn-------CGATTTGTTTACCC-CGAAGTCGTAAGCTGACCAACCACTCTACCCAACCGGTCCTAGATGACTGGTGTCATACAATCGTCGTACCGATGATAGACTGCGG
```
This cluster has 5 different unique sequences in it.  Three of them are only one copy (shown by the size=1 flag in the header).

The first two sequences are the only ones with any high numbers.  With the current settings, pyRAD is treating this as a paralog because the PE reads have 7 gaps in the alignment.  The default setting is to only allow 3 indels.  To improve this assembly, we will likely need to increase the setting.  Let's change it to 10.

```bash 
sed -i '/## 27./c\10,99               ## 27. maxIndels: within-clust,across-clust (def. 3,99) ' ./params.txt
```
Now, let's delete all the initial cluster files and redo this step

```bash 
rm ./clust.88/* && mv ./stats/s3.clusters.txt ./stats/s3.clusters.txt.old
pyrad -p params.txt -s 3
```
Let's check the results
```bash 
head -6 ./stats/s3.clusters.txt
```
```
taxa	total	dpt.me	dpt.sd	d>5.tot	d>5.me	d>5.sd	badpairs
PopA_01	901		19.829	9.083	852		20.782	8.396	46
PopA_02	901		19.91	8.81	860		20.702	8.214	61
PopA_03	903		20.474	9.506	862		21.283	8.956	55
PopA_04	925		19.564	9.085	887		20.268	8.599	49
```
This looks better, but still not ideal.  I leave it to you to experiment further.  With real data, you will again have to make a judgement call.  Keeping looking at the alignments in the clust88 directory and let them be your guide.
You can also alter the percentage of similarity parameter to cluster by as well.  It's option  in the params.txt file.  Another option to consider is the minimum number of read pairs to form a cluster.  The default is 6. and controlled by option #8 in the params.txt.  For the rest of this example, I am going to use a minimum coverage of 3 and a gap limit of 20.
```bash 
sed -i '/## 27./c\20,99                ## 27. maxIndels: within-clust,across-clust (def. 3,99) ' ./params.txt
sed -i '/## 8./c\3                     ## 8. Mindepth: min coverage for a cluster ' ./params.txt
```	
The next step of the pyRAD assembly calls the consensus sequence for each within-individual cluster.  It also applies filters aiming to remove potential paralogs.

It does this by estimating the error rate and level of heterozygosity in the data set and filters clusters that have too many heterozygous sits, more than 2 haplotypes, and too many low quality bases
```bash 
pyrad -p params.txt -s 45
```
The output should be like this:
```
     ------------------------------------------------------------
      pyRAD : RADseq for phylogenetics & introgression analyses
     ------------------------------------------------------------


        step 4: estimating error rate and heterozygosity
        ........................................
        step 5: created consensus seqs for 40 samples, using H=0.00732 E=0.00100
        ........................................
```
The next step is to cluster between samples
```bash 
pyrad -p params.txt -s 6
```
The output on the screen should look like:
```
     ------------------------------------------------------------
      pyRAD : RADseq for phylogenetics & introgression analyses
     ------------------------------------------------------------


        step 6: clustering across 40 samples at '.88' similarity

vsearch v1.1.1_linux_x86_64, 883.4GB RAM, 160 cores
https://github.com/torognes/vsearch

Reading file /gdc_home4/jpuritz/test/D1W/pyrad/clust.88/cat.firsts_ 100%  
3115443 nt in 33487 seqs, min 91, max 99, avg 93
Indexing sequences 100%  
Counting unique k-mers 100%  
Clustering 100%  
Writing clusters 100%  
Clusters: 1050 Size min 1, max 40, avg 31.9
Singletons: 19, 0.1% of seqs, 1.8% of clusters

	finished clustering
[```
We can see that pyRAD (via the program vsearch) found 1049 different shared reference sequences

Next we call the last step of pyRAD to produce usable outputs of all the data
```bash 
pyrad -p params.txt -s 7
```
The screen should look like:
```
  	------------------------------------------------------------
    	pyRAD : RADseq for phylogenetics & introgression analyses
   	------------------------------------------------------------

	ingroup PopA_01,PopA_02,PopA_03,PopA_04,PopA_05,PopA_06,PopA_07,PopA_08,PopA_09,PopA_10,PopA_11,PopA_12,PopA_13,PopA_14,PopA_15,PopA_16,PopA_17,PopA_18,PopA_19,PopA_20,PopB_00,PopB_01,PopB_02,PopB_03,PopB_04,PopB_05,PopB_06,PopB_07,PopB_08,PopB_09,PopB_10,PopB_11,PopB_12,PopB_13,PopB_14,PopB_15,PopB_16,PopB_17,PopB_18,PopB_19
	addon 
	exclude 
	................................................................
	final stats written to:
	 /gdc_home4/jpuritz/test/D1W/pyrad/stats/c88d6m4p3.stats
	output files being written to:
	 /gdc_home4/jpuritz/test/D1W/pyrad/outfiles/ directory
```
Let's take a look at the stats.

```bash 
head ./stats/c88d6m4p3.stats 
```
```
1018        ## loci with > minsp containing data
108         ## loci with > minsp containing data & paralogs removed
108         ## loci with > minsp containing data & paralogs removed & final filtering

## number of loci recovered in final data set for each taxon.
taxon	nloci
PopA_01	68
PopA_02	67
```
What the heck happened to all our data?  We went from 1018 RAD fragments to 106???????
It looks like pyRAD is inferring that almost all of the loci are paralogs.
Remember, pyRAD is designed to generate phylogenetic data sets and is not default configured to deal with highly polymorphic populations.
Setting number 13 sets the maximum number of individuals with a shared heterozygous site.  The default configuration is only 3.  
In a population we expect that heterozygosity maxes out at 50%.  In this simulated data, we have two populations of 20 individuals each, and with little genetic structure between them.  Let's try setting this to 20 and rerunning step 7.
```bash 
sed -i '/## 13./c\20                     ## 13. MaxSH: max inds with shared hetero site ' ./params.txt
rm ./outfiles/* && pyrad -p params.txt -s 7
```	
Let's see if that helped.
```bash
head ./stats/c88d6m4p3.stats
```
```
1018        ## loci with > minsp containing data
970         ## loci with > minsp containing data & paralogs removed
970         ## loci with > minsp containing data & paralogs removed & final filtering

## number of loci recovered in final data set for each taxon.
taxon	nloci
PopA_01	806
PopA_02	802
```	
That looks much better! 970 is very close to the actual value!
Now that you know how to manipulate the different parameters in pyRAD, experiment on your own to see if you can find the right settings to get to the correct number of loci!

##Bonus
Want to play with PyRAD more?  Try adding more outputs via line # in the params.txt file
Check out the general use tutorial and paired ddRAD tutorial here http://dereneaton.com/software/pyrad/
