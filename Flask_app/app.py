from flask import Flask, request, jsonify
import boto3
import json

app = Flask(__name__)

# Initialize a Lambda client
lambda_client = boto3.client('lambda')

@app.route('/predict', methods=['POST'])
def predict():
    # Get the SMILES string from the request
    data = request.json
    smiles = data.get('smile_string')

    # Invoke the Lambda function
    response = lambda_client.invoke(
        FunctionName='model_prediction_lambda',  # Replace with your Lambda function name
        InvocationType='RequestResponse',
        Payload=json.dumps({'smile_string': smiles,'descriptor_names':None})
    )

    # Read the response from the Lambda function
    response_payload = json.loads(response['Payload'].read())
    
    # Return the response as JSON
    return jsonify(response_payload)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
