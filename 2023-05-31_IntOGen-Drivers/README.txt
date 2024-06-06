# IntOGen RELEASE 20230531
# Headers of the driver files 
# For further information about methods please visit the documentation website: https://intogen.readthedocs.io/en/latest/index.html
# To access the source code of intogen please visit: https://bitbucket.org/intogen/intogen-plus/src/master/


Compendium_Cancer_Genes.tsv

1   SYMBOL: Hugo Symbol of the gene. 
2   TRANSCRIPT: ENSEMBL TRANSCRIPT identifier. 
3   COHORT: Name of the cohort where the gene has been detected as driver. 
4   CANCER_TYPE: Tumor type of the cohort. 
5   METHODS: Methods that detect this gene as driver (q-value <0.1). When the method is combination it means that the combination of p-values renders the gene as a driver. 
6   MUTATIONS: Number of somatic mutations (including short indels) of the gene in the cohort. 
7   SAMPLES: Number of mutated samples of the gene in the cohort. 
8   QVALUE_COMBINATION: Significance of the combined output. 
9   ROLE: Consensus mode of action of the gene. Derived from intOGen and from literature.  LoF (Loss of Function), Act (Activating) or amb (Ambiguous).  
10  CGC_GENE: Gene annotated in the Cancer Gene Census. 
11	CGC_CANCER_GENE: Gene annotated in the Cancer Gene Census as driver in the CANCER_TYPE of the COHORT. 	
12  DOMAINS: Comma-separated domains that are considered significant by smRegions (q-value <0.1).  The format is PFAM_ID:START_AA:END_AA. 
13  2D_CLUSTERS: Comma-separated linear clusters that are considered significant by OncodriveCLUSTL (p-value <0.05).  The format is START_AA:END_AA. 
14  3D_CLUSTERS: Comma-separated 3D clusters that are considered significant by HOTMAPS (q-value <0.05).  The format is AA_1, AA_2, AA_3, etc.
15  EXCESS_MIS: dNdScv missense excess rate.  
16  EXCESS_NON:  dNdScv nonsense excess rate. 
17  EXCESS_SPL:  dNdScv splicing excess rate.


Unfiltered_driver_results_05.tsv

1	SYMBOL: Hugo Symbol of the gene. 
2	TRANSCRIPT: Ensembl canonical transcript ID
3	COHORT: Name of the cohort where the gene has been detected as driver. 
4	CANCER_TYPE: Tumor type of the cohort. 
5	MUTATIONS: Number of somatic mutations (including short indels) of the gene in the cohort. 
6	SAMPLES_COHORT: Number of mutated samples of the gene in the cohort.
7	ALL_METHODS: Methods that detect this gene as driver.
8	SIG_METHODS: Methods that detect this gene as driver after the combination (q-value <0.1). When the method is combination it means that the combination of p-values renders the gene as a driver. 
9	QVALUE_COMBINATION: Significance of the combined output for a non-cgc gene. 
10	QVALUE_CGC_COMBINATION: Significance of the combined output for a cgc gene. 
11	RANKING: Ranking position at the combination. 
12	TIER: intOGen TIER classification. See documentation for information about the tiers classification. 
13	ROLE: Raw mode of action of the gene derived from intOGen. LoF (Loss of Function), Act (Activating) or amb (Ambiguous).  
14	CGC_GENE: Whether the gene is reported in Cancer Gene Census (CGC). 
15	TIER_CGC: Tier from Cancer Gene Census. Not Applicable if not reported in CGC. 
16	CGC_CANCER_GENE: Whether the gene has been reported as driver in the tumor type. See documentation of intOGen about how did we map tumor types. 
17	SIGNATURE9: Percentage of mutations attributed to Signature9. See documentation about detecting hypermutating events with signatures. 
18	SIGNATURE10: Percentage of mutations attributed to Signature10. See documentation about detecting hypermutating events with signatures. 
19	WARNING_EXPRESSION: Percentage of mutations attributed to Signature9. See documentation about 
20	WARNING_GERMLINE: True if the gene is labelled as polymorphic (see documentation for information about methods). False otherwise. 
21	SAMPLES_3MUTS: Number of samples with more than 2 mutations in the gene. 
22	OR_WARNING: Olfactory receptor.  
23	WARNING_ARTIFACT: True if the gene is labelled as likely an artifact (see artifacts.txt file from the repository to see the detailed list of genes included in the file). 
24	KNOWN_ARTIFACT: True if the gene is an artifact reported in the literature. 
25	NUM_PAPERS: Number of papers citing this gene as driver according to CancerMine. 
26	WARNING_ENSEMBL_TRANSCRIPTS: True if more than one Ensembl canonical transcript is assigned to the gene SYMBOL
27	DRIVER: True if passing the driver filters. See documentation about the driver filters applied. 
28	FILTER: Category of the according to their filters: 
        {"PASS": The gene is considered as driver, 
        "Germline Warning": The gene is filtered due to germline warning, 
        "Known artifact": The gene is filtered because is considered as an artifact, 
        "Warning artifact": The gene is filtered because is suspected to be an artifact, 
        "Lack of literature evidence": The gene is discarded because there is lack of literature evidence and lack of support from the methods, 
        "Samples with more than 2 mutations": The gene is filtered because there are samples with more than two mutations and it is not a CGC gene, 
        "Warning expression": The gene is filtered because is likely non-expressed in that tissue, 
        "Warning Signature9": The gene is filtered because is likely target of somatic hypermutation by AID, 
        "No driver": The gene is not considered driver by the driver filters (False in DRIVER column)} 