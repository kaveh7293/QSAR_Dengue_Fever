FROM public.ecr.aws/lambda/python:3.8

# Install system dependencies
RUN yum update -y && \
    yum install -y gcc gcc-c++ make && \
    yum clean all

# Copy function code
COPY lambda_function.py ${LAMBDA_TASK_ROOT}

# Install Python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Set the CMD to your handler
CMD [ "lambda_function.lambda_handler" ]
