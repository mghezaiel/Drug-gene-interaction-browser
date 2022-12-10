from abc import ABC, abstractmethod 


class GraphProcessingABC(ABC): 
    
    @abstractmethod
    def load_network(): 
        pass 
    
    @abstractmethod 
    def get_edgelist(): 
        pass 
    
    @abstractmethod 
    def get_ego_graph(): 
        pass 
    
    @abstractmethod 
    def add_labels(): 
        pass
    
    @abstractmethod 
    def get_k_hop_graph(): 
        pass 
    
    @abstractmethod 
    def filter_on_interaction_score(): 
        pass 
    
    @abstractmethod 
    def filter_on_interaction_type(): 
        pass 
    
    
class LoaderABC(ABC): 
    
    @abstractmethod 
    def get_data(): 
        pass 

    
    @abstractmethod
    def load_data(): 
        pass 

class DataFrameUtilsABC(ABC): 
    
    @abstractmethod 
    def get_drug_chembl_id(): 
        pass 
    
    @abstractmethod 
    def get_dgidb_edgelist(): 
        pass 
    
    @abstractmethod 
    def get_mint_edgelist(): 
        pass 
    
    @abstractmethod 
    def get_gene_names_from_drug_names(): 
        pass 
    
    @abstractmethod
    def save_hgnc_list(): 
        pass 
    
    @abstractmethod 
    def get_dgidb_df_from_uniprot_ids(): 
        pass 
    
    @abstractmethod 
    def filter_on_dgidb_interaction_score(): 
        pass
    
    @abstractmethod
    def get_dgidb_min_max_interaction_score():
        pass 
    
    @abstractmethod
    def group_mint_k_hop_nodes_df(): 
        pass 
    
    @abstractmethod 
    def save_uniprot_accession_list(): 
        pass 
    
    @abstractmethod
    def get_min_max_interaction_score(): 
        pass 
    
    @abstractmethod 
    def get_interaction_types_list(): 
        pass 
    
    
class ParserABC(): 
    
    @abstractmethod 
    def parse_mint_interactions(): 
        pass 
    
    @abstractmethod 
    def parse_drug_identifiers(): 
        pass 
    
    
class DrugInfoRetrieverABC(): 
    
    @abstractmethod
    def retriever(): 
        pass 
    
    @abstractmethod 
    def get_drug_infos(): 
        pass 
    
class AppABC(): 
    
    @abstractmethod
    def load_dgidb_graph(): 
        pass
    
    @abstractmethod
    def load_mint_graph(): 
        pass
    
    @abstractmethod
    def get_dgidb_k_hop_graph(): 
        pass
    
    @abstractmethod
    def get_target_drug_info(): 
        pass
    
    @abstractmethod
    def dgidb_hgnc_symbols_to_uniprot_accession(): 
        pass
    
    @abstractmethod
    def filter_dgidb_graph(): 
        pass
    
    @abstractmethod
    def get_mint_k_hop_graph():
        pass
    
    @abstractmethod
    def filter_mint_graph(): 
        pass
    
class CacheABC(): 
    
    @abstractmethod 
    def save_object(): 
        pass 
    
    @abstractmethod 
    def load_object(): 
        pass 
    
class DrawNetworkABC(): 
    
    @abstractmethod 
    def get_edge_trace(): 
        pass 
    
    @abstractmethod 
    def get_node_trace(): 
        pass 
    
    @abstractmethod 
    def get_figure(): 
        pass 
    
class ImageProcessingABC(): 
    
    @abstractmethod 
    def render_molecule_image(): 
        pass
    