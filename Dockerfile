FROM public.ecr.aws/docker/library/python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libgl1-mesa-glx \
    gcc \
    g++ \
    libxcb-xinerama0 \
    libxcb1 \
    libxcb-render0 \
    libxcb-render-util0 \
    libxcb-shape0 \
    libxcb-randr0 \
    libxcb-xfixes0 \
    && rm -rf /var/lib/apt/lists/*

# Set environment variable to use system-installed libraries
ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH

WORKDIR /app

COPY requirements.txt .

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "Streamlit_app.py", "--server.address", "0.0.0.0"]
