from abstract import GraphProcessingABC 
import networkx as nx 
import numpy as np 
import pandas as pd 
import streamlit as st 

class GraphProcessing(GraphProcessingABC):
    
    # LOAD DATABASE IN GRAPH FORMATS AND PROCESS GRAPHS 
    
    def __init__(self,**kwargs): 
        self.__dict__.update(kwargs)
        assert hasattr(self,"database")
        assert hasattr(self,"data")
    
    def get_edgelist(self,edges):
        
        # GET EDGELIST FROM DATABASE ACCORDING TO EDGES 
        
        self.edgelist =  self.data.filter(edges,axis = 1).dropna()
        self.data = self.data.loc[self.edgelist.index]
        
    def load_network(self): 
        
        # LOAD GRAPH FROM EDGELIST
        
        self.g = nx.from_edgelist(self.edgelist.values)
        
    def add_labels(self,data = False,graph = False, edges = None,labels = None):
        
        # ADD EDGE LABELS
        
        edge_attrs = {} 

        if isinstance(data,bool): 
            data = self.data
            
        for idr,row in data.iterrows(): 
            edge = (row[edges[0]],row[edges[1]])
        
            
            if edge not in edge_attrs: 
                edge_attrs[edge]= {label:list() for label in labels}
                
            for label in labels:
                edge_attrs[edge][label].append(row[label])
        
        if isinstance(graph,bool):
            nx.set_edge_attributes(self.g,edge_attrs)
        else: 
            nx.set_edge_attributes(graph,edge_attrs)
            return graph
        

    def get_k_hop_graph(self,nodes,k_hop = 1): 
        
        # GET K HOP GRAPH BY MERGING EGO GRAPHS
        
        k_hop_graph = nx.Graph()
        for node in nodes: 
            try:
                curr_ego_graph = nx.ego_graph(self.g,node,radius = k_hop).copy()
                edge_attrs = {edge:{} for edge in curr_ego_graph.edges()}
                
                for edge in curr_ego_graph.edges():
                    edge_attrs[edge]["dgidb_uniprot_accession"]= [node]
                    
                nx.set_edge_attributes(curr_ego_graph,edge_attrs)
                k_hop_graph.add_edges_from(curr_ego_graph.edges(data = True))
                
            except:
                continue 
                
        return k_hop_graph
                
    def get_ego_graph(self,node = "PREDNISONE"):
        
        # GET EGO GRAPH 
        
        return nx.ego_graph(self.g,node)
    

    def get_min_max_interaction_score(self,nodes = False):
        
        # GET MIN AND MAX INTERACTION SCORES ACCORDING TO NODES
        
        if self.database == "dgidb":
            interaction_score = "interaction_group_score"
            
        if self.database == "mint": 
            interaction_score = "mi_score"
            
        
        if nodes: 
            g = self.g.subgraph(nodes).copy()
        else: 
            g = self.g

        interaction_scores = [float(edge[2][interaction_score][0]) for edge in g.edges(data = True)]
        return min(interaction_scores),max(interaction_scores)
    

    def get_interaction_types_list(self,nodes = False): 
        
        # GET LIST OF INTERACTION TYPES ACCORDING TO NODES
        
        if self.database == "dgidb":
            
            interaction_type = "interaction_types"
            
        if self.database == "mint": 
            interaction_type = "interaction_type"
        
        if nodes: 
            g = self.g.subgraph(nodes).copy()
        else: 
            g = self.g
        
        interaction_types = [edge[2][interaction_type][0] for edge in g.edges(data = True)]
        return np.unique(interaction_types).tolist()
    

    def filter_on_interaction_score(self,graph,value = False):
        
        # FILTER EDGES ON INTERACTION SCORE
        
        if self.database == "dgidb":
            interaction_score = "interaction_group_score"
            
        if self.database == "mint": 
            interaction_score = "mi_score"
            
        if not graph: 
            graph = self.g
        
        if not value: 
            return graph 
        
        filtered_edges = [tuple(edge[:2]) for edge in graph.edges(data = True) if float(edge[2][interaction_score][0])>=value]

        return graph.edge_subgraph(filtered_edges).copy()
    

    def filter_on_interaction_type(self,graph,value = False):
    
        # FILTER EDGES ON INTERACTION TYPE
        
        if self.database == "dgidb":
            
            interaction_type = "interaction_types"
            
        if self.database == "mint": 
            interaction_type = "interaction_type"
        
        if not graph: 
            graph = self.g
        
        if not value: 
            return graph 
        
        filtered_edges = [tuple(edge[:2]) for edge in graph.edges(data = True) if edge[2][interaction_type][0]== value]
        return graph.edge_subgraph(filtered_edges).copy()