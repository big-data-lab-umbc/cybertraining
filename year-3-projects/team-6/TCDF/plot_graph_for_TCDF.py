import igraph
import numpy as np
        




def visualize_static(adja_mat,feature_names,plotname="graph_static.pdf"):

    if isinstance(adja_mat,np.ndarray):
        adja_mat=adja_mat.tolist()
        
    is_empty_graph=np.sum(np.abs(adja_mat))==0
        
    gr=igraph.Graph.Weighted_Adjacency(adja_mat, mode="directed", attr="weight", loops=False)
    gr.vs["name"]=feature_names
    layout = gr.layout_kamada_kawai()
    
    visual_style = {}
    visual_style["vertex_size"] = 40
    visual_style["vertex_label_size"]=8
    visual_style["vertex_color"] = [220,220,220]
    visual_style["vertex_shape"] = "rectangle"
    visual_style["vertex_label"] = gr.vs["name"]
    if not is_empty_graph:
        visual_style["edge_label"] =gr.es["weight"]
    visual_style["edge_width"] = 2
    visual_style["layout"] = layout
    visual_style["bbox"] = (500, 500)
    visual_style["margin"] = 30
    igraph.plot(gr,plotname,**visual_style)

def visualize_static_unweighted(adja_mat,feature_names,plotname="graph_static.pdf"):

    if isinstance(adja_mat,np.ndarray):
        adja_mat=adja_mat.tolist()

    is_empty_graph=np.sum(np.abs(adja_mat))==0

    gr=igraph.Graph.Weighted_Adjacency(adja_mat, mode="directed", attr="weight", loops=False)
    gr.vs["name"]=feature_names
    layout = gr.layout_kamada_kawai()

    visual_style = {}
    visual_style["vertex_size"] = 40
    visual_style["vertex_label_size"]=8
    visual_style["vertex_color"] = [220,220,220]
    visual_style["vertex_shape"] = "rectangle"
    visual_style["vertex_label"] = gr.vs["name"]
    visual_style["edge_width"] = 2
    visual_style["layout"] = layout
    visual_style["bbox"] = (500, 500)
    visual_style["margin"] = 30
    igraph.plot(gr,plotname,**visual_style)






# for testing
if __name__ == "__main__":

    aa=np.zeros((12,12))
    aa[1,10]=1
    aa[3,4]=1
    aa[3,11]=1
    aa[4,7]=1
    aa[7,0]=1
    aa[7,6]=1
    aa[10,2]=1
    aa[11,3]=1
    
    graph_h0k4=np.array([[0,0,0,0,0,0,0,0,0,0,0,0], # test 2: hidden_layer=0, kernel_size=4
                        [0,0,0,0,0,0,0,0,0,0,1,0],
                  	[0,0,0,0,0,0,0,0,0,0,0,0],
                  	[0,0,0,0,1,0,0,0,0,0,0,1],
                  	[1,0,0,0,0,0,0,0,0,0,0,0],
                  	[0,0,0,0,0,0,0,0,0,0,0,0],
                  	[0,0,0,0,0,0,0,1,0,0,0,0],
                  	[1,0,0,0,0,0,1,0,0,0,0,0],
                  	[0,0,0,0,0,0,0,0,0,0,0,0],
                  	[0,0,0,0,0,0,0,0,0,0,0,0],
                  	[0,0,1,0,0,0,0,0,0,0,0,0],
                  	[0,0,0,1,0,0,0,0,0,0,0,0]])

    feature_names=['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m','sea_ice','CC','CW','GH']

    visualize_static(aa,feature_names,"static.png")
    visualize_static_unweighted(aa,feature_names,"static_unweighted.png")
