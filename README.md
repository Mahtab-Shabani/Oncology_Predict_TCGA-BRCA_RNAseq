# Project Report – Classification of Breast Cancer Samples (TCGA-BRCA)
## 1. Dataset Description
This project uses the TCGA-BRCA dataset from The Cancer Genome Atlas (TCGA), one of the largest genomic repositories available.
The dataset provides comprehensive molecular and clinical information for breast cancer patients.

### Data Files Used
We worked with two core files downloaded from the TCGA data portal:

#### 1) Gene Expression Matrix (RNA-seq – HiSeqV2)
•	File name: HiSeqV2
•	Data type: RNA-seq gene expression
•	Size: 20,530 genes × 1,218 samples
•	Rows represent genes; columns represent patient samples.
•	Sample types included:
o	Primary Tumor (code 01)
o	Solid Tissue Normal (code 11)
This dataset is extremely high-dimensional, which poses significant modeling challenges.
[Link for download](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)<br />

#### 2) Clinical Matrix
•	File name: TCGA.BRCA.sampleMap_BRCA_clinicalMatrix
•	Contains:
o	Sample IDs
o	Patient metadata
o	Basic clinical information
We mainly used this file to match sample identifiers and filter valid cases.
[Link for download](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)<br />


## 2. Initial Challenges
**Challenge 1 — Strong Class Imbalance**
Among 1,211 usable samples:
•	1,097 Tumor
•	114 Normal
This ~1:10 imbalance causes typical ML models to:
•	Predict "Tumor" for nearly everything
•	Achieve artificially high accuracy
•	Fail to detect Normal samples (very low recall)
This was the most critical challenge in the dataset.


**Challenge 2 — High Dimensionality (20,530 genes)**
Only ~1,200 samples but >20k genes.
This leads to:
•	Overfitting
•	High noise
•	Overly expensive computations
Dimensionality reduction was essential.


**Challenge 3 — Mixed Sample Types**
Some samples had codes other than 01/11 (e.g., metastatic, recurrent).
These were removed.


## 3. Preprocessing and Solutions
**Solution 1 — Feature Selection**
We computed the variance of each gene and selected the top 300 most variable genes.
20,530 genes → 300 features
This:
•	Reduced noise
•	Improved training speed
•	Increased model stability


**Solution 2 — Class Balancing (Undersampling)**
Training set before balancing:
•	Tumor = 878
•	Normal = 91
To balance:
•	Randomly undersampled 91 Tumor samples
•	Paired with 91 Normal samples
Train set after balancing:
•	182 samples (91 Tumor, 91 Normal)
This significantly improved the classifier’s ability to detect both classes.


**Solution 3 — Training Three Classification Models**
We built three separate models to compare performance:
1.	SVM with RBF kernel
2.	PCA (50 components) + SVM
3.	Random Forest (TreeBagger, 120 trees)

## 4. Final Results
**Model A — SVM (RBF Kernel)**
•	Accuracy: 97.52%
•	Precision: 1.000
•	Recall: 0.973
•	F1-score: 0.986
<img width="500" alt="TCGA_compare_SVM_counts" src="https://github.com/user-attachments/assets/b75210a8-d051-4d84-9c1b-0ae108822acb" /><img width="500" alt="TCGA_compare_SVM_norm" src="https://github.com/user-attachments/assets/9f6cf8b2-08a6-4059-9498-5d9e7724b704" />


**Model B — PCA (50 components) + SVM**
Performance was essentially identical to Model A:
•	Accuracy: 97.52%
•	F1-score: 0.986
<img width="500" alt="TCGA_compare_PCA_SVM_counts" src="https://github.com/user-attachments/assets/200ab70c-8db7-47d6-923a-b5ffcf4252b1" /><img width="500" alt="TCGA_compare_PCA_SVM_norm" src="https://github.com/user-attachments/assets/83d7feb5-da98-4247-bebf-9f67980d2733" />


**Model C — Random Forest (TreeBagger)**
Best performing model in the project
•	Accuracy: 98.35%
•	Precision: 1.000
•	Recall: 0.982
•	F1-score: 0.991  <br />
This model misclassified only a few tumor samples and never predicted Normal incorrectly, making it the strongest model.
<img width="500" alt="TCGA_compare_RF_counts" src="https://github.com/user-attachments/assets/9d45a2cc-19fe-41c4-8cd1-8cc750b161f8" /><img width="500" alt="TCGA_compare_RF_norm" src="https://github.com/user-attachments/assets/88200f8a-e201-432a-a713-d50ee0c9493c" />



## 5. Summary

In this project: <br />
✔️ We processed and analyzed the TCGA-BRCA gene expression dataset <br />
✔️ Addressed challenges of high dimensionality and extreme class imbalance <br />
✔️ Performed feature selection (top 300 variable genes) <br />
✔️ Balanced the training set using undersampling <br />
✔️ Built and compared three ML models  <br />
✔️ Achieved a high classification performance, with Random Forest reaching 98.35% accuracy <br />
These results demonstrate that gene-expression–based classification of Tumor vs Normal samples can be highly accurate with appropriate preprocessing and model selection.
