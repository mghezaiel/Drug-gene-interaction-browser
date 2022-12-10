
FROM r-base:latest

# Author
LABEL maintainer.author.name = Morad Ghezaiel 
LABEL maintainer.email = ghezaiel.morad@gmail.com 

# Requirements and folders
COPY python-requirements.txt /tmp/python-requirements.txt
COPY classes /app/classes 
COPY data /app/data 
COPY cache /app/cache

# Dependencies for devtools (R)
RUN apt-get update && apt-get install -y libfontconfig1-dev libgit2-dev libharfbuzz-dev libfribidi-dev libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

# Install R dependencies 
RUN Rscript /app/classes/requirements.R

# Install pip 
RUN apt-get install -y python3-pip

# Upgrade pip
RUN pip install --upgrade pip

# Install python
RUN apt-get install -y python3-dev 

# Install psutil
RUN pip3 install psutil==5.4.7

# Install python requirements 
RUN pip install -r /tmp/python-requirements.txt

# Working directory 
WORKDIR /app/classes

# Commands to run the app
EXPOSE 8501
ENTRYPOINT ["streamlit", "run", "main.py", "--server.port=8501","--server.address=0.0.0.0"]
