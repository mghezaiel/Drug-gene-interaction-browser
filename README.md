![alt text](https://github.com/mghezaiel/Drug-gene-interaction-browser/blob/master/Banner.png)
# Drug gene interaction browser 

A web application to browse drug gene interactions, link genes to known protein/protein interactions and query uniprot annotations. 

The tool uses the Drug gene interaction database (DGIDB) to identify genes interacting with a target drug. 

Interacting genes can be filtered by interaction type, interaction score, and linked to protein/protein interactions from databases such as MINT. 

## Utility
- Identification of potential target proteins
- Evaluation of iatrogenic effects
- Identification of synergistic drug combinations 
- Drug repurposing

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

### Run without docker 
In /classes : ```streamlit run main.py```
