# Project Report – Classification of Breast Cancer Samples (TCGA-BRCA)
## 1. Dataset Description
This project uses the TCGA-BRCA dataset from The Cancer Genome Atlas (TCGA), one of the largest genomic repositories available.
The dataset provides comprehensive molecular and clinical information for breast cancer patients.

### Data Files Used
We worked with two core files downloaded from the TCGA data portal:

#### 1) Gene Expression Matrix (RNA-seq – HiSeqV2)
* File name: ``` HiSeqV2 ``` <br />
* Data type: RNA-seq gene expression <br />
* Size: 20,530 genes × 1,218 samples <br />
* Rows represent genes; columns represent patient samples. <br />
* Sample types included: <br />
•	Primary Tumor (code 01) <br />
•	Solid Tissue Normal (code 11) <br />

This dataset is extremely high-dimensional, which poses significant modeling challenges.
[Link for download](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)<br />

#### 2) Clinical Matrix
* File name: ``` TCGA.BRCA.sampleMap_BRCA_clinicalMatrix ``` <br />
* Contains: <br />
•	Sample IDs <br />
•	Patient metadata <br />
•	Basic clinical information <br />

We mainly used this file to match sample identifiers and filter valid cases. 
[Link for download](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)<br />


## 2. Initial Challenges
### Challenge 1 — Strong Class Imbalance <br />
Among 1,211 usable samples: <br />
•	`1,097 Tumor` <br />
•	`114 Normal` <br />

This ~1:10 imbalance causes typical ML models to: <br />
* Predict "Tumor" for nearly everything <br />
* Achieve artificially high accuracy <br />
* Fail to detect Normal samples (very low recall) <br />

This was the most critical challenge in the dataset. <br />



### Challenge 2 — High Dimensionality (20,530 genes) <br />
Only ~1,200 samples but >20k genes. <br />
This leads to: <br />
* Overfitting <br />
* High noise <br />
* Overly expensive computations <br />

Dimensionality reduction was essential. <br />


### Challenge 3 — Mixed Sample Types <br />
Some samples had codes other than 01/11 (e.g., metastatic, recurrent). <br />
These were removed. <br />


## 3. Preprocessing and Solutions
### Solution 1 — Feature Selection <br />
We computed the variance of each gene and selected the top 300 most variable genes. <br />
`20,530 genes → 300 features` <br />
This: <br />
*	Reduced noise <br />
*	Improved training speed <br />
*	Increased model stability <br />


### Solution 2 — Class Balancing (Undersampling) <br />
Training set before balancing: <br />
•	`Tumor = 878` <br />
•	`Normal = 91` <br />

To balance: <br />
*	Randomly undersampled 91 Tumor samples <br />
*	Paired with 91 Normal samples <br />

Train set after balancing: <br />
•	`182 samples (91 Tumor, 91 Normal)` <br />

This significantly improved the classifier’s ability to detect both classes. <br />


### Solution 3 — Training Three Classification Models <br />
We built three separate models to compare performance: <br />
1.	SVM with RBF kernel <br />
2.	PCA (50 components) + SVM <br />
3.	Random Forest (TreeBagger, 120 trees) <br />

## 4. Final Results
### Model A — SVM (RBF Kernel) <br />
•	`Accuracy: 97.52%` <br />
•	`Precision: 1.000` <br />
•	`Recall: 0.973` <br />
•	`F1-score: 0.986` <br />
<img width="500" alt="TCGA_compare_SVM_counts" src="https://github.com/user-attachments/assets/b75210a8-d051-4d84-9c1b-0ae108822acb" /><img width="500" alt="TCGA_compare_SVM_norm" src="https://github.com/user-attachments/assets/9f6cf8b2-08a6-4059-9498-5d9e7724b704" />


### Model B — PCA (50 components) + SVM <br />
Performance was essentially identical to Model A: <br />
•	`Accuracy: 97.52%` <br />
•	`F1-score: 0.986` <br />


<img width="500" alt="TCGA_compare_PCA_SVM_counts" src="https://github.com/user-attachments/assets/200ab70c-8db7-47d6-923a-b5ffcf4252b1" /><img width="500" alt="TCGA_compare_PCA_SVM_norm" src="https://github.com/user-attachments/assets/83d7feb5-da98-4247-bebf-9f67980d2733" />


### Model C — Random Forest (TreeBagger) <br />
Best performing model in the project <br />
•	`Accuracy: 98.35%` <br />
•	`Precision: 1.000` <br />
•	`Recall: 0.982` <br />
•	`F1-score: 0.991`  <br />
This model misclassified only a few tumor samples and never predicted Normal incorrectly, making it the strongest model. <br />
<img width="500" alt="TCGA_compare_RF_counts" src="https://github.com/user-attachments/assets/9d45a2cc-19fe-41c4-8cd1-8cc750b161f8" /><img width="500" alt="TCGA_compare_RF_norm" src="https://github.com/user-attachments/assets/88200f8a-e201-432a-a713-d50ee0c9493c" />



## 5. Summary

In this project: <br />
✔️ We processed and analyzed the TCGA-BRCA gene expression dataset <br />
✔️ Addressed challenges of high dimensionality and extreme class imbalance <br />
✔️ Performed feature selection (top 300 variable genes) <br />
✔️ Balanced the training set using undersampling <br />
✔️ Built and compared three ML models  <br />
✔️ Achieved a high classification performance, with Random Forest reaching `98.35% accuracy` <br />
These results demonstrate that gene-expression–based classification of `Tumor vs Normal` samples can be highly accurate with appropriate preprocessing and model selection.


## 6. How to Run the Code

* Run the `demo_TCGA_BRCA_compare_models.m` script (Results and figures will be generated automatically).

Dependency
------------
This code is implemented in MATLAB 2014a and doesn't depend on any other toolbox.
