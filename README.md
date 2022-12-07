# Drug gene interaction browser 

A web application to browse drug gene interactions, link genes to known protein/protein interactions and inspect the annotations. 

The tool uses the Drug gene interaction database (DGIDB) to identify genes interacting with a target drug. 
Interacting genes can be filtered by interaction type, interaction score, and linked to protein/protein interactions from databases such as MINT. 
It also allows to query uniprot annotations for the filtered gene sets. 

## Utility
- Identification of potential protein targets 
- Iatrogenesis
- Synergistic drug combination 

# RUN 

## Start the docker daemon 

On linux: 
	sudo systemctl docker start 
	sudo service docker start 

On mac: 
	open -a Docker 

On windows:
	restart-service *docker*


## Build the image 

docker build -t dgib -f Dockerfile .

## Run the app 
docker run -p 8501:8501 dgib


