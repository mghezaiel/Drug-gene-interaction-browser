from abstract import DrawNetworkABC, ImageProcessingABC
import plotly 
import plotly.graph_objects as go
import numpy as np 
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import DrawingOptions

class DrawNetwork(DrawNetworkABC):
    
    # DRAW GRAPHS AS PLOTLY IMAGES
    
    def get_edge_trace(G,layout):
        
        # PLOT EDGES
        
        edge_x = [] 
        edge_y = [] 

        for edge in G.edges(): 
            x0,y0= layout[edge[0]]
            x1,y1 = layout[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

        return edge_trace
 
    
    def get_node_trace(G,layout,node_colors): 
        
        # PLOT NODES 
        
        node_x = [] 
        node_y = [] 

        for node in G.nodes(): 
            x,y= layout[node]

            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
        x=node_x, y=node_y,text = [node for node in G.nodes()],
        mode='markers',
        hoverinfo="text",
        marker=dict(
            showscale=False,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=node_colors,
            size=10,
            line_width=2))
        return node_trace 

    def get_figure(G,layout,node_colors,height,width,output_path = False): 
        
        # PLOT EDGES AND NODES
        edge_trace= DrawNetwork.get_edge_trace(G,layout)
        node_trace= DrawNetwork.get_node_trace(G,layout,node_colors)

        # WHITE BACKGROUND
        figure_layout = go.Layout(
        plot_bgcolor='rgb(255,255,255)',

        )
        
        # MERGE EDGES AND NODES IN A FIGURE
        fig = go.FigureWidget(data = [edge_trace, node_trace],layout = figure_layout)
        
        #fig.update_layout(hovermode="x")
        fig.update_layout(showlegend=False)

        #x axis
        fig.update_xaxes(visible=False)

        #y axis    
        fig.update_yaxes(visible=False)
        fig.update_layout(
        autosize=False,
        width=width,
        height=height)
        
        if output_path: 
            fig.write_image(output_path)
            
        else: 
            return fig

class ImageProcessing(ImageProcessingABC):
    
    # DRAW DRUG CHEMICAL STRUCTURES
    
    def render_molecule_image(smiles,path):
        
        # DRAW CHEMICAL STRUCTURE USING RDKIT 
        
        mol = Chem.MolFromSmiles(smiles)
        rdDepictor.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)

        drawer.drawOptions().useCDKAtomPalette()
        drawer.SetLineWidth(2)
        drawer.SetFontSize(1.0)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '')
        SVG(svg)
        
        # SAVE IMAGE
        with open(path, "w") as f:
            f.write(svg)
    