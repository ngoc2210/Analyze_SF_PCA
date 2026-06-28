## Goals:
- I want to verify if the PCA phase 2 of SF_GWAS is correct by implementing a simple python script mimic the same logic of phase 2 of SF_GWAS. All main information about phase 2 PCA
of SF_GWAS you should find in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\sfgwas\gwas\gwas.go, C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\sfgwas\gwas\pca.go, (you may need the other important information for understanding the logics in other files in the same folder C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\sfgwas)
This phase should take C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\out\geno_pca_party1.tsv and C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\out\geno_pca_party2.tsv as inputs , which are encoded from the true genotype data after Quality Control filtering and LD filtering, and output something like C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\cache\party1\Qpc1.txt and
C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\cache\party2\Qpc2.txt, The out put should be placed in out C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\out
## How to apply these goals:
- Your code in python and simple as in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\.github\copilot-instructions.md description
- Your code should be placed in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\eigenvec_sf.py, you should use the available functions and variables for example computing mean, stdInv and matmul_lazynorm and avoid fixing the available things , except really needed.
- You should use the exact logic and parameters from taking geno_pca_party1.tsv as encoded from geno_pca.bin of party 1 and geno_pca_party2.tsv as encoded from geno_pca.bin of party 2 to output Qpc1.txt and Qpc2.txt from phase 2 SF_GWAS, the only diffenrence is you
dont use any MHE or MPC mechanism, just do the operations in cleartext.
- If the the output files of SF_GWAS is not in the same domain as in the cleartext what you compute, then please could you transform it the the cleartext domain and 
compare this with the one you get in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\eigenvec_sf.py by using angle metric in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\Analyze_SF_PCA\comparison.py

copilot --resume=f9487f9f-091f-47b5-af2e-60c06f83c556
