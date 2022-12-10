from utils import Loader, DataFrameUtils, Parser, DrugInfoRetriever,Cache
from data_visualization import DrawNetwork, ImageProcessing
from graph_processing import GraphProcessing
import os 
import numpy as np 
import itertools 
import pandas as pd 
import matplotlib.pyplot as plt 
import networkx as nx 
import streamlit as st
from PIL import Image
import time 
import datetime

def App(data_dir,cache_dir): 
        
    # EDGE NAMES AND EDGE LABELS FOR THE DATABASES
    dgidb_edges = ["gene_name","drug_name"]
    dgidb_labels = ["interaction_group_score","interaction_types","entrez_id","chembl_id"]
    mint_edges = ["uniprot_A","uniprot_B"]
    mint_labels = ["mi_score","interaction_type","interaction_id","taxid_A","taxid_B","PMID"]
        
    # PATHS FOR CACHED OBJECTS
    loader_path = os.path.join(cache_dir,"loader")
    parser_path = os.path.join(cache_dir,"parser")
    dgidb_interactions_df_path = os.path.join(cache_dir,"dgidb_interactions")
    dgidb_graph_path = os.path.join(cache_dir,"dgidb_graph")
    mint_interactions_df_path = dgidb_target_drug_metadata_path = os.path.join(cache_dir,"mint_interactions")
    mint_graph_path = os.path.join(cache_dir,"mint_graph")
    
    ## DATA LOADER (CACHED)
    if not os.path.exists(loader_path):
        with st.spinner(text="Downloading databases ..."):
            data = Loader(data_dir = data_dir)
            data.get_data()
            data.load_data()
            Cache.save_object(loader_path,data)
    else: 
        data = Cache.load_object(loader_path)
        
    ## PARSER (CACHED)
    if not os.path.exists(parser_path):
        parser = Parser(data_dir = data_dir)
        Cache.save_object(parser_path,parser)
    else: 
        parser = Cache.load_object(parser_path)
        
    ## DRUG/GENE INTERACTION DF (CACHED)
    if not os.path.exists(dgidb_interactions_df_path):
        dgidb_interactions = Parser.parse_drug_identifiers(data.interactions)
        Cache.save_object(dgidb_interactions_df_path,dgidb_interactions)
    else: 
        dgidb_interactions = Cache.load_object(dgidb_interactions_df_path)
    
    # SIDEBAR TO CHOOSE PARAMETERS
    sidebar = st.sidebar
    
    # TABS TO DISPLAY INTERACTION GRAPHS
    tab1, tab2 = st.tabs(["DGIDB","MINT"])
    
    ## DGIDB INTERACTIONS    
    with tab1:
        
        # DRUG NAMES
        drugs = np.unique([drug for drug in dgidb_interactions["drug_name"].values])
        default = drugs.tolist().index("CAFFEINE")
        
        # DGIDB GRAPH
        if not os.path.exists(dgidb_graph_path):
            dgidb_graph = GraphProcessing(data = dgidb_interactions, database = "dgidb")
            dgidb_graph.get_edgelist(edges = dgidb_edges)
            dgidb_graph.load_network()
            dgidb_graph.add_labels(edges = dgidb_edges,labels = dgidb_labels)
            dgidb_graph.g.remove_edges_from(nx.selfloop_edges(dgidb_graph.g))
            Cache.save_object(dgidb_graph_path,dgidb_graph)
        else: 
            dgidb_graph = Cache.load_object(dgidb_graph_path)

        # DRUG SELECTION 
        target_drug = sidebar.selectbox("Drug",drugs,index = default)
        
        # DGIDB TARGET DRUG GRAPH PATH
        dgidb_target_drug_graph_path = os.path.join(cache_dir,f"dgidb_{target_drug}_graph")
        
        # GET TARGET DRUG EGO GRAPH
        if not os.path.exists(dgidb_target_drug_graph_path):
            dgidb_target_drug_graph = dgidb_graph.get_ego_graph(node = target_drug)
            dgidb_target_drug_graph.remove_edges_from(nx.selfloop_edges(dgidb_target_drug_graph))
            Cache.save_object(dgidb_target_drug_graph_path,dgidb_target_drug_graph)
            
        else: 
            dgidb_target_drug_graph = Cache.load_object(dgidb_target_drug_graph_path)
        
        # DGIDB TARGET DRUG HGNC SYMBOLS LIST
        dgidb_target_drug_hgnc = [hgnc for hgnc in list(dgidb_target_drug_graph.nodes()) if hgnc!=target_drug] 

        # GET DRUG INFORMATION FROM CHEMBL WEBRESSOURCE CLIENT
        target_drug_chembl_id = int([edge[2]["chembl_id"][0] for edge in dgidb_target_drug_graph.edges(data = True)][0])
        drug_infos_path = os.path.join(cache_dir,f"{target_drug}_drug_info")
        
        if not os.path.exists(drug_infos_path):
            drug_infos = DrugInfoRetriever.get_drug_infos(cache_dir,target_drug)
            Cache.save_object(drug_infos_path,drug_infos)
            
        else: 
            drug_infos = Cache.load_object(drug_infos_path)

        # INPUT OUTPUT PATH FOR HGNC SYMBOL TO UNIPROT ACCESSION CONVERSION
        input_path = os.path.join(data_dir,
                              f"dgidb_{target_drug_chembl_id}_hgnc_symbols.csv")
        output_path = os.path.join(data_dir,
                               f"dgidb_{target_drug_chembl_id}_uniprot_accessions.csv") 

        # SAVE HGNC LIST 
        df = pd.DataFrame(dgidb_target_drug_hgnc)
        df.columns = ["hgnc_symbol"]
        df.to_csv(input_path)

        # CONVERT HGNC SYMBOLS TO UNIPROT ACCESSIONS
        if not os.path.exists(output_path):
            os.system(f"Rscript get_uniprot_accessions_from_hgnc_symbols.R -i {input_path} -o {output_path}")
        
        # ADD UNIPROT ACCESSION TO DGIDB TARGET DRUG INTERACTION DF
        dgidb_target_drug_metadata_path = os.path.join(cache_dir,f"dgidb_{target_drug}_metadata")
        
        if not os.path.exists(dgidb_target_drug_metadata_path):
            dgidb_target_drug_metadata= DataFrameUtils.get_dgidb_df_from_uniprot_accessions(data = dgidb_graph.data,
                                                                      uniprot_accessions_path = output_path)
            Cache.save_object(dgidb_target_drug_metadata_path,dgidb_target_drug_metadata)
        else: 
            dgidb_target_drug_metadata = Cache.load_object(dgidb_target_drug_metadata_path)

        # ADD UNIPROT ACCESSIONS TO DGIDB TARGET DRUG GRAPH LABEL LIST
        dgidb_target_drug_graph_with_uniprot_accessions_path = os.path.join(cache_dir,f"dgidb_{target_drug}_with_uniprot_accessions")
        
        if not os.path.exists(dgidb_target_drug_graph_with_uniprot_accessions_path):
            
            dgidb_target_drug_graph = dgidb_graph.add_labels(data = dgidb_target_drug_metadata, 
                                                       graph = dgidb_target_drug_graph, 
                                                        edges = dgidb_edges,labels = ["uniprotswissprot"]) 
            
            Cache.save_object(dgidb_target_drug_graph_with_uniprot_accessions_path,dgidb_target_drug_graph)
        else: 
            dgidb_target_drug_graph = Cache.load_object(dgidb_target_drug_graph_with_uniprot_accessions_path)
            
        # GET DGIDB TARGET DRUG INTERACTION TYPE LIST 
        dgidb_interaction_types_list = dgidb_graph.get_interaction_types_list(list(
            dgidb_target_drug_graph.nodes()))
        dgidb_interaction_types_list.append("All types")
        dgidb_interaction_types_list = np.unique([interaction_type for interaction_type in dgidb_interaction_types_list])
        dgidb_interaction_types_list = [interaction_type for interaction_type in dgidb_interaction_types_list
                                      if interaction_type!="nan"]
        
        # DISPLAY A SELECT BOX TO SELECT AN INTERACTION TYPE
        selected_interaction_type = sidebar.selectbox("Interaction type", options = dgidb_interaction_types_list)

        if selected_interaction_type=="All types": 
            selected_interaction_type = False 

        # FILTER DGIDB TARGET DRUG GRAPH TO THE SELECTED INTERACTION TYPE
        dgidb_target_drug_graph_filtered_on_interaction_type = dgidb_graph.filter_on_interaction_type(dgidb_target_drug_graph,
                                                                            value = selected_interaction_type)

        # GET MIN MAX INTERACTION SCORE FROM THE FILTERED DGIDB TARGET DRUG GRAPH
        dgidb_target_drug_min_max_interaction_scores_path = os.path.join(cache_dir,
                                                    f"dgidb_{target_drug}_min_max_interaction_score")
        
        if not os.path.exists(dgidb_target_drug_min_max_interaction_scores_path):
            
            # GET MIN AND MAX DGIDB INTERACTION SCORE 
            dgidb_min_interaction_score, dgidb_max_interaction_score = dgidb_graph.get_min_max_interaction_score(
            nodes = list(dgidb_target_drug_graph_filtered_on_interaction_type.nodes()))
            Cache.save_object(dgidb_target_drug_min_max_interaction_scores_path,
                              (dgidb_min_interaction_score, dgidb_max_interaction_score))
            
        else: 
            dgidb_min_interaction_score, dgidb_max_interaction_score = Cache.load_object(dgidb_target_drug_min_max_interaction_scores_path)

        # DISPLAY A SLIDER TO SELECT AN INTERACTION SCORE THRESHOLD
        step = np.round(0.1*np.abs(dgidb_max_interaction_score-dgidb_min_interaction_score),3)

        # FILTER THE DGIDB TARGET DRUG GRAPH TO INTERACTION SCORES > threshold
        if dgidb_min_interaction_score!=dgidb_max_interaction_score:
            sidebar.header("Filter Drug/gene interaction graph")
            interaction_score_threshold = sidebar.slider('Interaction score',min_value = dgidb_min_interaction_score,
                                                     max_value = dgidb_max_interaction_score,
                                                     )
            # DEFAULT THRESHOLD
            if interaction_score_threshold == dgidb_min_interaction_score:
                interaction_score_threshold = False

            # FILTER DGIDB TARGET DRUG GRAPH TO INTERACTION SCORE > THRESHOLD
            dgidb_target_drug_graph_filtered_on_interaction_score = dgidb_graph.filter_on_interaction_score(
                dgidb_target_drug_graph_filtered_on_interaction_type, 
                value = interaction_score_threshold)           
        else: 
            
            # DISPLAY DEFAULT DGIDB TARGET DRUG GRAPH INTERACTION SCORE (MIN)
            st.caption(f"Interaction score: {dgidb_max_interaction_score}")
            dgidb_target_drug_graph_filtered_on_interaction_score = dgidb_target_drug_graph_filtered_on_interaction_type
        
        # FILTERED DGIDB TARGET DRUG NODE LIST
        filtered_dgidb_target_drug_nodes = np.unique([edge[2]["uniprotswissprot"][0] 
                                            for edge in dgidb_target_drug_graph_filtered_on_interaction_type.edges(data = True) 
                                            if "uniprotswissprot" in edge[2].keys()]) 

        # DISPLAY DRUG INFORMATION
        with st.container():
            st.header("Target drug")
            c1,c2 = st.columns([1,1])
            drug_infos = pd.DataFrame.from_dict(drug_infos,orient = "index").fillna("Not listed")
            drug_infos[0] = drug_infos[0].apply(lambda x: x.upper())
            drug_infos.rename({0:" "},axis = 1,inplace = True)
            drug_infos.index = [idx.replace("_"," ") for idx in drug_infos.index]
            target_drug_chemical_structure_image = os.path.join(cache_dir,target_drug+".svg")
            blank_image = os.path.join(data_dir,"blank.svg")
            
            
            if not os.path.exists(target_drug_chemical_structure_image): 
                c1.caption("Chemical structure not available")
                c1.image(blank_image)
            else:
                c1.caption("Chemical structure")
                c1.image(target_drug_chemical_structure_image)
            c2.caption("Drug information")
            c2.dataframe(drug_infos)

        # DISPLAY TARGET DRUG/PROTEIN INTERACTION GRAPH
        with st.container():
            st.header("Interactions")
            # SAVE DGIDB K_HOP GRAPH AS IMAGE
            dgidb_target_drug_graph_image = os.path.join(cache_dir,target_drug+f"dgidb_graph.png")

            dgidb_layout = nx.kamada_kawai_layout(dgidb_target_drug_graph_filtered_on_interaction_type)
            node_colors = ["red" if node==target_drug else "blue" 
                           for node in dgidb_target_drug_graph_filtered_on_interaction_type.nodes()]

            dgidb_graph_figure = DrawNetwork.get_figure(dgidb_target_drug_graph_filtered_on_interaction_type,
                                                        dgidb_layout,node_colors, height = 500, width = 500)
            st.caption("Drug/gene interaction graph")
            f = st.plotly_chart(dgidb_graph_figure,use_container_width = True)
            plt.close()

    with tab2:
        
        sidebar.header("Filter MINT interaction graph")
        
        # SELECT THE NUMBER OF HOP 
        mint_k_hop = sidebar.slider("K-hop protein neighbors", min_value = 1, max_value = 5, step = 1)
        
        # GET MINT INTERACTION DF
        mint_interactions_df_path = dgidb_target_drug_metadata_path = os.path.join(cache_dir,"mint_interactions")
        if not os.path.exists(mint_interactions_df_path):
            mint_interactions = parser.parse_mint_interactions()
            Cache.save_object(mint_interactions_df_path,mint_interactions)
        else: 
            mint_interactions = Cache.load_object(mint_interactions_df_path)
        
        # GET MINT GRAPH 
        if not os.path.exists(mint_graph_path):
            mint_graph = GraphProcessing(data = mint_interactions, database = "mint")
            mint_graph.get_edgelist(mint_edges)
            mint_graph.load_network()
            mint_graph.add_labels(edges = mint_edges, labels = mint_labels)
            mint_graph.g.remove_edges_from(nx.selfloop_edges(mint_graph.g))
            Cache.save_object(mint_graph_path,mint_graph)
            
        else: 
            mint_graph = Cache.load_object(mint_graph_path)

        # GET MINT K HOP GRAPH
        mint_k_hop_graph_path = os.path.join(cache_dir,f"dgidb_{target_drug}_mint_graph_{mint_k_hop}_hop_graph")

        if not os.path.exists(mint_k_hop_graph_path):
            mint_k_hop_graph  = mint_graph.get_k_hop_graph(nodes = filtered_dgidb_target_drug_nodes, 
                                                                k_hop = mint_k_hop)
            Cache.save_object(mint_k_hop_graph_path,mint_k_hop_graph)
        else: 
            mint_k_hop_graph = Cache.load_object(mint_k_hop_graph_path)

        # GET DGIDB TARGET DRUG UNIPROT ID FROM MINT GRAPH EDGE LABELS
        dgidb_target_drug_node_list_path = os.path.join(cache_dir,f"dgidb_{target_drug}_mint_graph_dgidb_node_list")

        if not os.path.exists(dgidb_target_drug_node_list_path):
            dgidb_target_drug_node_list = np.unique([edge[2]["dgidb_uniprot_accession"][0] 
                                               for edge in mint_k_hop_graph.edges(data = True)])

            dgidb_target_drug_node_list = [(uniprot_accession,9606) for uniprot_accession in dgidb_target_drug_node_list]
            Cache.save_object(dgidb_target_drug_node_list_path,dgidb_target_drug_node_list)

        else: 
            dgidb_k_hop_node_list = Cache.load_object(dgidb_target_drug_node_list_path)

        # GET MINT K HOP GRAPH NODE LIST
        mint_k_hop_node_list_path = os.path.join(cache_dir,f"dgidb_{target_drug}_mint_graph_{mint_k_hop}_hop_mint_node_list")

        if not os.path.exists(mint_k_hop_node_list_path):

            mint_k_hop_node_list = np.array([(edge[0],edge[2]["taxid_A"][0],edge[1],edge[2]["taxid_B"][0]) 
                                for edge in mint_k_hop_graph.edges(data = True)]) 

            mint_k_hop_node_list = mint_k_hop_node_list.reshape(2*mint_k_hop_node_list.shape[0],2)
            Cache.save_object(mint_k_hop_node_list_path,mint_k_hop_node_list)

        else: 
            mint_k_hop_node_list = Cache.load_object(mint_k_hop_node_list_path)

        # GET MINT K_HOP INTERACTION TYPES
        mint_k_hop_interaction_types_list = mint_graph.get_interaction_types_list(
            list(mint_k_hop_graph.nodes()))
        mint_k_hop_interaction_types_list.insert(0,"All types")

        selected_interaction_type = sidebar.selectbox("Interaction type", mint_k_hop_interaction_types_list)

        if selected_interaction_type=="All types": 
            selected_interaction_type = False 

        # FILTER MINT K_HOP GRAPH ON INTERACTION TYPE 
        mint_k_hop_graph_filtered_on_interaction_type = mint_graph.filter_on_interaction_type(mint_k_hop_graph,
                                                                                value = selected_interaction_type)
        
        mint_k_hop_min_max_interaction_scores_path = os.path.join(cache_dir,
                                                            f"dgidb_{target_drug}_mint_graph_{mint_k_hop}_hop_mint_min_max_interaction_scores")

        # GET MINT K HOP MIN MAX INTERACTION SCORES (SLIDER)
        mint_min_interaction_score, mint_max_interaction_score = mint_graph.get_min_max_interaction_score(
            list(mint_k_hop_graph_filtered_on_interaction_type.nodes()))

        if mint_min_interaction_score!=mint_max_interaction_score:
            step = np.round(0.1*np.abs(mint_max_interaction_score-mint_min_interaction_score),3)

            interaction_score_threshold = sidebar.slider('Interaction score',min_value = mint_min_interaction_score,
                                                     max_value = mint_max_interaction_score,
                                                    )   
            if interaction_score_threshold == mint_min_interaction_score:
                interaction_score_threshold = False

            # FILTER MINT K_HOP GRAPH ON INTERACTION SCORE
            mint_k_hop_graph_filtered_on_interaction_score = mint_graph.filter_on_interaction_score(mint_k_hop_graph_filtered_on_interaction_type,
                                                                                     value = interaction_score_threshold)
        else: 
            sidebar.caption(f"Interaction score: {mint_min_interaction_score}")
            mint_k_hop_graph_filtered_on_interaction_score = mint_k_hop_graph_filtered_on_interaction_type

        mint_k_hop_graph_filtered_on_interaction_score.remove_edges_from(nx.selfloop_edges(
                         mint_k_hop_graph_filtered_on_interaction_score))
        
        # DISPLAY MINT K HOP GRAPH AS IMAGE
        mint_k_hop_graph_image = os.path.join(cache_dir,target_drug+f"mint_{mint_k_hop}_hop_graph.png")

        mint_layout = nx.kamada_kawai_layout(mint_k_hop_graph_filtered_on_interaction_score)
        mint_node_colors = ["blue" if node==target_drug else "blue" 
                       for node in mint_k_hop_graph_filtered_on_interaction_score.nodes() ]

        mint_graph_figure = DrawNetwork.get_figure(mint_k_hop_graph_filtered_on_interaction_score,
                                                    mint_layout,mint_node_colors,width = 700,height = 700)
        
        with st.container():
            st.header("Interactions")
            st.caption("Protein/protein interaction graph")
            f = st.plotly_chart(mint_graph_figure,use_container_width = True)
            plt.close()
                        
        # DOWNLOAD TARGET PROTEINS ANNOTATIONS FROM UNIPROT
        with st.container():
            if st.button("Get protein annotations"):
                with st.spinner("Downloading protein information from UniprotKB/Swiss-prot..."):

                    # GET DGIDB TARGET DRUG UNIPROT ID FROM MINT GRAPH EDGE LABELS                

                    dgidb_target_drug_node_list = np.unique([edge[2]["dgidb_uniprot_accession"][0] 
                                                       for edge in mint_k_hop_graph_filtered_on_interaction_score.edges(data = True)])

                    dgidb_target_drug_node_list = [(uniprot_accession,9606) for uniprot_accession in dgidb_target_drug_node_list]
                    
                    # GET MINT K HOP NODE LIST
                    mint_k_hop_node_list = np.array([(edge[0],edge[2]["taxid_A"][0],edge[1],edge[2]["taxid_B"][0]) 
                                            for edge in mint_k_hop_graph_filtered_on_interaction_score.edges(data = True)]) 

                    mint_k_hop_node_list = mint_k_hop_node_list.reshape(2*mint_k_hop_node_list.shape[0],2)

                    # PATH FOR MINT UNIPROT ACCESSION AND ANNOTATIONS
                    mint_k_hop_nodes_uniprot_path = os.path.join(data_dir,
                                    f"dgidb_{target_drug_chembl_id}_mint_{mint_k_hop}_hop_uniprot_accessions.csv")

                    mint_k_hop_nodes_annotations_path = os.path.join(data_dir,
                                    f"dgidb_{target_drug_chembl_id}_mint_{mint_k_hop}_hop_annotations.csv")

                    # SAVE MINT K HOP NODE LIST
                    DataFrameUtils.save_uniprot_accession_list(mint_k_hop_node_list,mint_k_hop_nodes_uniprot_path)

                    # RUN R SCRIPT TO QUERY UNIPROT
                    os.system(f"Rscript get_protein_infos.R -i {mint_k_hop_nodes_uniprot_path} -o {mint_k_hop_nodes_annotations_path}")
                    
                                
                # DISPLAY RESULTS 
                st.header("UniprotKB/Swiss-prot annotations")
                mint_k_hop_protein_annotations = pd.read_csv(mint_k_hop_nodes_annotations_path).fillna("Unknown")
                mint_k_hop_protein_annotations.set_index("From",inplace = True, drop = True)
                columns = ["Gene.Names","Organism","Protein.existence","Protein.families"]
                mint_k_hop_protein_annotations = mint_k_hop_protein_annotations.filter(columns, axis = 1)
                mint_k_hop_protein_annotations.set_axis(["Gene names", "Organism", "Protein existence", "Protein family"],
                                                     axis =1, 
                                                     inplace = True)                
                current_dateTime = datetime.datetime.now().strftime("%Y_%m_%d %HH%M")
                st.dataframe(mint_k_hop_protein_annotations)
                output = mint_k_hop_protein_annotations.to_csv()
                st.download_button(
                label="Download results as CSV",
                data=output,
                file_name=f"results_{current_dateTime}",
                mime='text/csv',
            )
                
                    
                
        