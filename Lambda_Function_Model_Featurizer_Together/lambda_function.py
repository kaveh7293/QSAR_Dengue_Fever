import json
import boto3
import pickle
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


# Initialize AWS clients
s3 = boto3.client('s3')
lambda_client = boto3.client('lambda')

# S3 model information
bucket_name = "drug-discovery-files";#os.environ['bucket_name']
model_key = "random_forest_class.pickle";# os.environ['model_key']
Featurizer_function="arn:aws:lambda:us-east-1:277707143091:function:Featurizer_function";#os.environ['featurizer']


def featurize_molecule(smiles, descriptor_names=None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    features = {}
    
    if descriptor_names:
        for desc_name in descriptor_names:
            try:
                desc_func = getattr(Descriptors, desc_name)
                features[desc_name] = desc_func(mol)
            except AttributeError:
                print(f"Descriptor {desc_name} not found in  Descriptors.")
    
    return pd.DataFrame([features])



# def get_featurized_smile(smile_string):
#     """Invokes the featurizer Lambda function."""
#     response = lambda_client.invoke(
#         FunctionName=Featurizer_function,  # Replace with your featurizer Lambda function name
#         InvocationType='RequestResponse',
#         Payload=json.dumps({'smile_string': smile_string})
#     )
#     response_payload = json.loads(response['Payload'].read())
#     return json.loads(response_payload['body'])  # Adjust based on how the featurizer returns data

class load_prediction:
    def __init__(self):
        obj = s3.get_object(Bucket=bucket_name, Key=model_key)
        model = pickle.loads(obj['Body'].read())
        self.model = model
    def predict(self,smiles,features):
        df=featurize_molecule(smiles,features)
        return self.model.predict(df)
    


def lambda_handler(event, context):
    # 1. Extract SMILES string from the event
    smile_string = event['smile_string']
    descriptor_names = event.get('descriptor_names', None)
    # 2. Get the features by invoking the featurizer Lambda
    load_from_s3_and_predict = load_prediction()
    descriptor_names = [desc_name for desc_name, _ in Descriptors.descList]

    prediction=load_from_s3_and_predict.predict(smile_string,descriptor_names)

    
    # 3. Load the model from S3
    
    # 4. Predict using the model
    
    # 5. Return the prediction result
    return {
        'statusCode': 200,
        'body': json.dumps({'prediction': prediction.tolist()})
    }


