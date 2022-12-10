from abstract import LoaderABC, DataFrameUtilsABC, ParserABC, DrugInfoRetrieverABC, CacheABC
from data_visualization import ImageProcessing
import pandas as pd 
import os 
import re 
import numpy as np 
from tqdm import tqdm 
from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st 
import pickle
import networkx as nx 
from chembl_webresource_client.new_client import new_client 
from chembl_webresource_client.settings import Settings

#Settings.Instance().RESPECT_RATE_LIMIT = True
#Settings.Instance().TOTAL_RETRIES = 0
#Settings.Instance().MAX_LIMIT = 20

class Loader(LoaderABC):
    
    # Class to download the databases 
    
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)
        assert hasattr(self,"data_dir") 
        
        self.urls = {"categories.tsv":f"""curl -o {os.path.join(self.data_dir,"categories.tsv")} https://www.dgidb.org/data/monthly_tsvs/2022-Feb/categories.tsv""",
        "drugs.tsv":f"""curl -o {os.path.join(self.data_dir,"drugs.tsv")} https://www.dgidb.org/data/monthly_tsvs/2022-Feb/drugs.tsv -O""",
        "genes.tsv":f"""curl -o {os.path.join(self.data_dir,"genes.tsv")}  https://www.dgidb.org/data/monthly_tsvs/2022-Feb/genes.tsv""",
        "interactions.tsv":f"""curl -o {os.path.join(self.data_dir,"interactions.tsv")} https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv""",
        "mint.txt":f"""curl -o {os.path.join(self.data_dir,"mint.txt")} http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/species:human"""}
        
    def get_data(self):
        
        # DOWNLOAD THE DATABASES FROM URLS 
        for filename in self.urls: 
            if filename not in os.listdir(self.data_dir):
                os.system(self.urls[filename])

    def load_data(self):
        
        # LOAD DATABASES IN PD.DATAFRAMES AND SAVE IT AS CLASS ATTRIBUTES
        for filename in self.urls:
            if filename != "mint.txt": 
                df = pd.read_csv(os.path.join(self.data_dir,filename),sep = "\t")
            else: 
                df = pd.read_csv(os.path.join(self.data_dir,filename),sep = "\t", header = None)
            setattr(self,filename.split(".")[0],df)
               
class DataFrameUtils(DataFrameUtilsABC): 
                
    @staticmethod 
    def get_dgidb_df_from_uniprot_accessions(data,uniprot_accessions_path):
        
        # MERGE R SCRIPTS OUTPUT (UNIPROT ACCESSIONS) WITH DGIDB INTERACTION DF
        uniprot_accessions = pd.read_csv(uniprot_accessions_path)
        return pd.merge(data,uniprot_accessions,left_on = "gene_name", right_on = "hgnc_symbol")
    
    
    @staticmethod 
    def save_uniprot_accession_list(uniprot_accessions,output_dir):
        
        # SAVE UNIPROT ACCESSIONS 
        df = pd.DataFrame(uniprot_accessions) 
        df.columns = ["uniprot_accession","taxid"]
        df.to_csv(output_dir)
        
        
class Parser(ParserABC):
    
    # PARSER
    
    def __init__(self,**kwargs): 
        self.__dict__.update(kwargs) 
        assert hasattr(self,"data_dir")
        
    def parse_mint_interactions(self):
        
        # PARSE MINT INTERACTIONS 
        
        df = pd.read_csv(os.path.join(self.data_dir,"mint.txt"), sep = "\t", header = None) 
        df["uniprot_A"] = df[0].apply(lambda x: re.findall("(?<=uniprotkb:)[A-Z0-9]+",x))
        df["uniprot_B"] = df[1].apply(lambda x: re.findall("(?<=uniprotkb:)[A-Z0-9]+",x)) 
        df["mi_score"] = df[14].apply(lambda x: re.findall("[0-9.]+",x))
        df["uniprot_A"] = df["uniprot_A"].apply(lambda x: x[0] if len(x)>0 else np.nan)
        df["uniprot_B"] = df["uniprot_B"].apply(lambda x: x[0] if len(x)>0 else np.nan) 
        df["mi_score"] = df["mi_score"] .apply(lambda x: x[0] if len(x)>0 else np.nan)
        df["interaction_type"] = df[11].apply(lambda x: re.findall("(?<=\()[a-z ]+",x)[0])
        df["interaction_id"] = df[13].apply(lambda x: re.findall("EBI-[0-9]+",x)[0])
        df = df.dropna()
        df["A_mask"] = df["uniprot_A"].apply(lambda x: True if type(x[0])==str else False)
        df = df[df["A_mask"]]
        df["B_mask"] = df["uniprot_B"].apply(lambda x: True if type(x[0])==str else False)
        df = df[df["B_mask"]]
        df["taxid_A"] = df[9].apply(lambda x: re.findall("(?<=taxid:)[0-9]+",x))
        df["taxid_A"] = df["taxid_A"].apply(lambda x: x[0] if len(x)>0 else np.nan)
        df["taxid_B"] = df[10].apply(lambda x: re.findall("(?<=taxid:)[0-9]+",x))
        df["taxid_B"] = df["taxid_B"].apply(lambda x: x[0] if len(x)>0 else np.nan)
        df["PMID"] = df[8].apply(lambda x: re.findall("(?<=pubmed:)[0-9]+",x))
        df["PMID"] = df["PMID"].apply(lambda x: x[0] if len(x)>0 else "Unknown")
        
        return df.dropna()
    
    @staticmethod
    def parse_drug_identifiers(data):
        
        # GET CHEMBL-ID FROM COLUMNS 
        
        data["chembl_id"] = data["drug_concept_id"].apply(lambda x: re.findall("(?<=CHEMBL)[0-9]+",x) if type(x)==str else "Unknown")
        data["chembl_id"] = data["chembl_id"].apply(lambda x: x[0] if len(x)>0 else "Unknown")
        return data

class DrugInfoRetriever(DrugInfoRetrieverABC): 
    
    # RETRIEVE DRUG INFORMATION USING CHEMBL WEBRESSOURCE CLIENT
    
    @staticmethod
    def retriever(data_dir,drug_name):
        
        # RETRIEVE DRUG INFORMATION AND RENDER CHEMICAL STRUCTURE FROM SMILES

        molecule = new_client.molecule
        info = molecule.search(drug_name)[0]
        
        try:
            indication_class = info["indication_class"]
        except: 
            indication_class = np.nan 
        try:
            smiles = info["molecule_structures"]["canonical_smiles"]
        except: 
            smiles = np.nan 
        try:
            molecular_formula = info["molecule_structures"]["standard_inchi"]
            molecular_formula = molecular_formula.split("/")[1]
        except: 
            molecular_formula = np.nan 
        
        
        infos = {"name":drug_name,"smiles":smiles,"molecular_formula":molecular_formula,"indication_class":indication_class}
        image_path = os.path.join(data_dir,f'{drug_name}.svg')   
        
        if isinstance(smiles,str):  
            ImageProcessing.render_molecule_image(infos["smiles"],image_path)

        return infos
    
    @staticmethod
    def get_drug_infos(data_dir,drug_name):
        
        
        drug_info = DrugInfoRetriever.retriever(data_dir,drug_name) 
        return drug_info
        
        
class Cache(CacheABC): 
    
    # CACHE OBJECTS IN MEMORY 
    
    @staticmethod
    def save_object(path,object): 
        
        # SAVE OBJECT USING PICKLE
        
        if not os.path.exists(path): 
            pickle.dump(object,open(path,"wb"))
                 
    @staticmethod
    def load_object(path): 
        
        # UNPICKLE OBJECT
        
        return pickle.load(open(path,"rb"))
    