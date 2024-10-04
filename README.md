# QSAR_Dengue_Fever
Developing a Quantitative Structure-activity Relationship (QSAR) Model for Discovery of Drugs for Dengue_Fever.
Machine Learning Model Development Steps Summary

Steps followed to find the drug candidates:

1-Find the protein target from the ChEMBL Database. The following protein is used as a target:
CHEMBL5980 (Dengue virus type 2 NS3 protein).
2-Find the ligands associated with the NS3 protein and filter them based on the availability of their IC50.
<img width="644" alt="image" src="https://github.com/user-attachments/assets/7e2d2c43-ca42-4342-be04-e8693b9cf6e0">

3- Assigin activity levels to each of the ligands based on their IC50. The associated activity levels are called Active (i.e., IC50 less than 1000 nM), Intermediate (i.e., IC50 between 1000 to 10000 nM) and Inactive (i.e., IC50 greater than 10000 nM).

4-In the next step a featurizer function is created based on the features available in RDKit library. This function was used to generate features based on the SMILES structure of the molecules. In the following, you can see the resulting dataframe that will be used to train the associated models.

5-After 80/20 split of training and testing datasets, the following metrics were obtained for the trained RandomForestClassifier.

There are some important points needed to be noted. It seams that there is no Active ligands for our target. Therefore, the process of drug discovery will be limited to find intermedaited compounds. Also, since the dataset is imbalanced, we needed to make the database more balanced. I used the Synthetic Minority Oversampling Technique (SMOTE) to balance our database for the sake of finding the intermdiate ligands more accurately.


I deployed the model into the aws lambda but given the size of the python library I had to use a dockerized version of the lambda function in Amazon ECr and deploy the corresponding lambda functions using the image from ECR. This lambda function return the output to the model by only receiving the smiles structure of the drug candidates. An API has been created using flask python framework so that we can eventually use from our function in our front-end streamlit application. In the following you can see a summary of the architecture that I have used. 

<img width="435" alt="image" src="https://github.com/user-attachments/assets/5c65cdfd-ca24-436a-ac6f-74b8b9a43ea8">
