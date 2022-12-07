![alt text](https://github.com/mghezaiel/Drug-gene-interaction-browser/blob/master/5AF7A9B2-4183-4BBE-8454-D6A9A81EA79F_4_5005_c.jpeg)
# Drug gene interaction browser 

A web application to browse drug gene interactions, link genes to known protein/protein interactions and query protein annotations. 

The tool uses the Drug gene interaction database (DGIDB) to identify genes interacting with a target drug. 

Interacting genes can be filtered by interaction type, interaction score, and linked to protein/protein interactions from databases such as MINT. 

## What it can be used for
- Identification of potential protein targets 
- Identification of iatrogenesis
- Identification of Synergistic drug combination 

## Databases
- DGIDB 
- MINT 

## Tool 

### Start the docker daemon 

On linux: 
	```
	sudo systemctl docker start
	```
	or
	```
	sudo service docker start
	```
	
On mac: ```
	open -a Docker 
	```
	
On windows:
	```
	restart-service *docker*
	```

### Build the image 

```docker build -t dgib -f Dockerfile .```

### Run the app 
```docker run -p 8501:8501 dgib```


