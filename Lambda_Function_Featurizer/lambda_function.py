import json
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

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
                print(f"Descriptor {desc_name} not found in RDKit Descriptors.")
    
    return pd.DataFrame([features])

def lambda_handler(event, context):
    try:
        smiles = event['smiles']
        descriptor_names = event.get('descriptor_names', None)
        
        features_df = featurize_molecule(smiles, descriptor_names)
        features_dict = features_df.to_dict(orient='records')[0]
        
        response = {
            'statusCode': 200,
            'body': json.dumps(features_dict, default=str)
        }
    except Exception as e:
        response = {
            'statusCode': 400,
            'body': json.dumps({'error': str(e)})
        }
    
    return response


event={
  "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  "descriptor_names": ["MolWt", "MolLogP", "NumHDonors", "NumHAcceptors"]
}
print(lambda_handler(event, None))