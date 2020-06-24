import numpy as np
import pandas
import os.path
from compute_SHD import count_accuracy
from scipy.spatial.distance import hamming

distance_function='l1'
#distance_function='HD'

lambda_range=[0,0.1]
w_thresh_range=[0.2,0.3]

file_names=['combined_decomposed_drop_temp_all_norm_1980_2018']

file_names_for_table={'combined_decomposed_drop_temp_all_norm_1980_2018':'combined\\_decomposed\\_drop\\_temp\\_all\\_norm\\_1980\\_2018'}

path_to_temporal='/umbc/xfs1/cybertrn/cybertraining2020/team6/research/code_Matt/Final_June8/Temporal/'

feature_names_dict={'combined_decomposed_drop_temp_all_norm_1980_2018': ['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m','sea_ice','CC','CW','GH']}

full_feature_names_dict={'combined_decomposed_drop_temp_all_norm_1980_2018': ['Residual_heat_flux','Residual_shortwave','Residual_longwave','Residual_SLP',
                                                                                  'Residual_tot_precip','Residual_RH','Residual_u10m','Residual_v10m','Residual_sea_ice',
                                                                                   'Residual_cloud_cover','Residual_cloud_water','Residual_GH_mean']}





if distance_function=='l1':
    
    def dist_between_matrices(A,B):
        return np.sum(np.abs(A-B))


if distance_function=='SHD':

    def dist_between_matrices(A,B):
        
        _, _, _, shd, _ =count_accuracy(A,B)

        return shd


if distance_function=='HD':

    def dist_between_matrices(A,B):

        A_unweighted=np.zeros(A.shape)
        B_unweighted=np.zeros(B.shape)

        A_unweighted[A!=0]=1
        B_unweighted[B!=0]=1
        
        return hamming(A_unweighted.flatten(),B_unweighted.flatten())




def reduced_graph(adja_mat,file):

    df = pandas.read_csv('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/data/lagged_data_May_12/'+'lagged_'+file+'.csv', delimiter=',',header=0)
    all_variable_names=list(df.columns)
    feature_names=feature_names_dict[file] 
    full_feature_names=full_feature_names_dict[file]                                                     
    lag_range=range(1,12+1)                                                                                                                                           

    d=len(feature_names)
    adja_mat_reduced=np.zeros((d,d))

    n=len(all_variable_names)
    for ell in range(n):
        for mmm in range(n):
            if adja_mat[ell,mmm]!=0:
                
                temp1=all_variable_names[ell]
                temp2=all_variable_names[mmm]
                
                if not temp1 in full_feature_names:
                    if temp1[-3]=='-':
                        temp1=temp1[:-3]
                    else:
                        temp1=temp1[:-2]

                if not temp2 in full_feature_names:
                    if temp2[-3]=='-':
                        temp2=temp2[:-3]
                    else:
                        temp2=temp2[:-2]
                
                adja_mat_reduced[full_feature_names.index(temp1),full_feature_names.index(temp2)]=1

    return adja_mat_reduced



       
#Order of variables
#['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m','sea_ice','CC','CW','GH']
graph_Yiyi=np.array([[0,0,0,1,1,0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0,0,1,0,0,0],
                     [0,0,0,0,0,0,0,0,1,0,0,0],
                     [1,0,0,0,0,1,1,1,1,0,0,1],
                     [1,0,1,0,0,1,0,0,1,1,1,0],
                     [0,0,1,0,1,0,0,0,0,1,1,0],
                     [1,0,0,0,0,0,0,0,1,0,0,0],
                     [1,0,0,0,0,0,0,0,1,0,0,0],
                     [1,1,1,1,0,0,1,1,0,0,0,0],
                     [1,1,1,0,1,1,0,0,0,0,0,0],
                     [1,1,1,0,1,1,0,0,0,0,0,0],
                     [0,0,1,1,0,1,0,0,0,0,0,0]])



for file in file_names:

    dmat=np.zeros((len(lambda_range)*len(w_thresh_range),len(lambda_range)*len(w_thresh_range)))
    distance_to_ground_truth=np.zeros(len(lambda_range)*len(w_thresh_range))
    
    dmat_TEMPORAL=np.zeros((len(lambda_range)*len(w_thresh_range),len(lambda_range)*len(w_thresh_range)))
    distance_to_ground_truth_TEMPORAL=np.zeros(len(lambda_range)*len(w_thresh_range))

    dmat_STATIC_vs_TEMPORAL=np.zeros((len(lambda_range)*len(w_thresh_range),len(lambda_range)*len(w_thresh_range)))

    no_dag_indi=np.zeros(len(lambda_range)*len(w_thresh_range))
    no_dag_indi_TEMPORAL=np.zeros(len(lambda_range)*len(w_thresh_range))
    para_comb=[]
    para_comb_TEMPORAL=[]
   
    plotname="StaticModel__"+file
    plotname_TEMPORAL="TemporalModel__"+'lagged_'+file
    tablename="StaticAndTemporalModel --- "+file_names_for_table[file]
    

    #Static
    iter_A=0
    for lambda1_A in lambda_range:
        for w_threshold_A in w_thresh_range:
            
            iter_A+=1
            para_comb.append((lambda1_A,w_threshold_A))

            if os.path.isfile('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'!!!NODAG!!!.csv'):
                no_dag_A=True
                no_dag_indi[iter_A-1]=1
                A=np.loadtxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'!!!NODAG!!!.csv',delimiter=',')
            else:
                no_dag_A=False
                A=np.loadtxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'.csv',delimiter=',')
            
            
            distance_to_ground_truth[iter_A-1]=dist_between_matrices(graph_Yiyi,A)


            iter_B=0
            for lambda1_B in lambda_range:
                for w_threshold_B in w_thresh_range:

                    iter_B+=1

                    if os.path.isfile('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'!!!NODAG!!!.csv'):
                        no_dag_B=True
                        B=np.loadtxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'!!!NODAG!!!.csv',delimiter=',')
                    else:
                        no_dag_B=False
                        B=np.loadtxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'.csv',delimiter=',')

        
                    dmat[iter_A-1,iter_B-1]=dist_between_matrices(A,B)
    
                    

    #Temporal
    iter_A=0
    for lambda1_A in lambda_range:
        for w_threshold_A in w_thresh_range:

            iter_A+=1
            para_comb_TEMPORAL.append((lambda1_A,w_threshold_A))

            if os.path.isfile(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'!!!NODAG!!!.csv'):
                no_dag_A=True
                no_dag_indi_TEMPORAL[iter_A-1]=1
                A=np.loadtxt(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'!!!NODAG!!!.csv',delimiter=',')
            else:
                no_dag_A=False
                A=np.loadtxt(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'.csv',delimiter=',')
                
            A=reduced_graph(A,file)
            distance_to_ground_truth_TEMPORAL[iter_A-1]=dist_between_matrices(graph_Yiyi,A)


            iter_B=0
            for lambda1_B in lambda_range:
                for w_threshold_B in w_thresh_range:

                    iter_B+=1

                    if os.path.isfile(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'!!!NODAG!!!.csv'):
                        no_dag_B=True
                        B=np.loadtxt(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'!!!NODAG!!!.csv',delimiter=',')
                    else:
                        no_dag_B=False
                        B=np.loadtxt(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'.csv',delimiter=',')

                    B=reduced_graph(B,file)
                    dmat_TEMPORAL[iter_A-1,iter_B-1]=dist_between_matrices(A,B)



    #Static vs Temporal                                                                                                                         
    iter_A=0
    for lambda1_A in lambda_range:
        for w_threshold_A in w_thresh_range:

            iter_A+=1

            if os.path.isfile('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'!!!NODAG!!!.csv'):
                no_dag_A=True
                A=np.loadtxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'!!!NODAG!!!.csv',delimiter=',')
            else:
                no_dag_A=False
                A=np.loadtxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1_A)+'_Wthresh='+str(w_threshold_A)+'.csv',delimiter=',')


            iter_B=0
            for lambda1_B in lambda_range:
                for w_threshold_B in w_thresh_range:

                    iter_B+=1

                    if os.path.isfile(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'!!!NODAG!!!.csv'):
                        no_dag_B=True
                        B=np.loadtxt(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'!!!NODAG!!!.csv',delimiter=',')
                    else:
                        no_dag_B=False
                        B=np.loadtxt(path_to_temporal+'AdjMatrix_'+plotname_TEMPORAL+'_lambda1='+str(lambda1_B)+'_Wthresh='+str(w_threshold_B)+'.csv',delimiter=',')

                    B=reduced_graph(B,file)
                    dmat_STATIC_vs_TEMPORAL[iter_A-1,iter_B-1]=dist_between_matrices(A,B)



    print(plotname)
    print(dmat)
    print(dmat_TEMPORAL)
    print(dmat_STATIC_vs_TEMPORAL)

    
    
    
    #Produce table for report
    if distance_function=='SHD':
        digits_to_round=0
    else:
        digits_to_round=2

    f = open("FinalDistMatrix_distance="+distance_function+".txt", "w")
    f.write("\\begin{table}[h]\n")
    f.write("\caption{"+tablename+" --- distance="+distance_function+"}\n")
    f.write("\\vspace{4mm}\n")
    f.write("\centering\n")
    f.write("\\begin{small}\n")
    f.write("\\renewcommand{\\arraystretch}{1.2}\n")
    f.write("\\begin{tabular}{cc|"+''.join(["c"]*len(para_comb))+''.join(["c"]*len(para_comb_TEMPORAL))+"}\n")
    
    f.write("& & \multicolumn{"+str(len(para_comb))+"}{c}{\emph{Static}} & \multicolumn{"+str(len(para_comb_TEMPORAL))+"}{c}{\emph{Temporal}} \\\ \\rule{0pt}{\\abstA} \n")
    f.write("&")
    for cc,ell in enumerate(para_comb):
        f.write(" & $\lambda="+str(ell[0])+"$")
    for cc,ell in enumerate(para_comb_TEMPORAL):
        f.write(" & $\lambda="+str(ell[0])+"$")
    f.write("\\\ \n")
    f.write("&")
    for cc,ell in enumerate(para_comb):
        f.write(" & $t="+str(ell[1]))
        if no_dag_indi[cc]==1:
            f.write("^\star$")
        else:
            f.write("$")
    for cc,ell in enumerate(para_comb_TEMPORAL):
        f.write(" & $t="+str(ell[1]))
        if no_dag_indi_TEMPORAL[cc]==1:
            f.write("^\star$")
        else:
            f.write("$")
    f.write("\\\ \\hline \n")
    
    f.write("\\parbox[t]{2mm}{\multirow{"+str(len(para_comb))+"}{*}{\\rotatebox[origin=c]{90}{\emph{Static}}}}")
    for cc,ell in enumerate(para_comb):
        f.write(" & $\lambda="+str(ell[0])+", t="+str(ell[1]))
        if no_dag_indi[cc]==1:
            f.write("^\star$")
        else:
            f.write("$")
        for mmm in range(len(para_comb)):
            f.write(" & "+str(np.around(dmat[cc,mmm], digits_to_round)))
        for mmm in range(len(para_comb_TEMPORAL)):
            f.write(" & "+str(np.around(dmat_STATIC_vs_TEMPORAL[cc,mmm], digits_to_round)))
        f.write("\\\ \\rule{0pt}{\\abstA} \n")

    f.write("\\parbox[t]{2mm}{\multirow{"+str(len(para_comb_TEMPORAL))+"}{*}{\\rotatebox[origin=c]{90}{\emph{Temporal}}}}")
    for cc,ell in enumerate(para_comb_TEMPORAL):
        f.write(" & $\lambda="+str(ell[0])+", t="+str(ell[1]))
        if no_dag_indi_TEMPORAL[cc]==1:
            f.write("^\star$")
        else:
            f.write("$")
        for mmm in range(len(para_comb)):
            f.write(" & "+str(np.around(dmat_STATIC_vs_TEMPORAL[mmm,cc], digits_to_round)))
        for mmm in range(len(para_comb_TEMPORAL)):
            f.write(" & "+str(np.around(dmat_TEMPORAL[cc,mmm], digits_to_round)))
        f.write("\\\ \\rule{0pt}{\\abstA} \n")


    f.write("& Ground truth")
    for mmm in range(len(para_comb)):
        f.write(" & "+str(np.around(distance_to_ground_truth[mmm], digits_to_round)))
    for mmm in range(len(para_comb_TEMPORAL)):
        f.write(" & "+str(np.around(distance_to_ground_truth_TEMPORAL[mmm], digits_to_round)))
    f.write("\\\ \n")
    f.write("\end{tabular}\n")
    f.write("\end{small}\n")
    f.write("\end{table}")
    f.close()

