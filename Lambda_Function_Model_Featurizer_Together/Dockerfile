# Use an image that supports Python 3.12.5
FROM public.ecr.aws/lambda/python:3.12

# Copy function code
COPY lambda_function.py ${LAMBDA_TASK_ROOT}

# Copy requirements file
COPY requirements.txt .

# Install the function's dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Set the CMD to your handler
CMD [ "lambda_function.lambda_handler" ]
