# SepsisData

This repository holds all raw data that have been analyzed in:

**"Plasma mass spectrometry fingerprints induced by bacterial endotoxins in murine models as markers of pathophysiological evolution in sepsis."**

Authors: Martín Ledesma, María Florencia Todero, Lautaro Maceira, Monica Prieto, Carlos Vay, Marcelo Galas, Beariz López, Noemí Yokobori, Bárbara Rearte.  

This repository contains two folders with raw MALDI-TOF-MS (Folder.1 and Folder.2) data of plasma mice samples under different LPS stimuli that emulate different human sepsis phases. It also includes 2 CSV folders, one that contains the metadata (Metadata.csv) of the former folders and the other that covers all the bioinformatic assigned proteins (Assigned_proteins.csv). 

* Metadata.csv:

This table shows sample names (Names), the sepsis phase condition of each sample (Class), the number of spectra replicates regarding each sample (N.spectra), the folder in which those samples were present (Archive), and finally, the number code of each sample in the Folder raw data (Code). This last Code should be interpreted as follow: the first sample (Cnt1.Exp1) present one replicate, it appears in Folder.1, and the Code is 10; this means that folder 10 in Folder.1 belongs to sample Cnt1.Exp1. The second sample (Cnt2.Exp1) present two replicates, it appears in Folder.1, and the Code is 11; this means that folder 11 and 12 in Folder.1 belongs to sample Cnt2.Exp1. 

* Assigned_proteins.csv: 

This table presents the experimental mass of each peak (PM.peak), the condition where each peak is characteristic (Group
.expressed), the ranking within each algorithm (Rankingn.BDA and Rankingd.Rf), the bioinformatic PM assigned (PM.assigned), and the other columns show uniprot features of each hit, including links to the assignment. 
