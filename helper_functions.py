import pandas as pd
from gprofiler import GProfiler
import re
import numpy as np
import math
import matplotlib.pyplot as plt


# pathway enrichment for genelist 'genes' probing 'databases'
def pathway_enrich(genes, databases):
   
    gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

    cluster_enrichment = gp.profile(organism='mmusculus', sources=databases, user_threshold=0.05,
                               significance_threshold_method='fdr', 
                               domain_scope ='annotated',
                               #background= 10000, 
                               query= genes) #"contains the list of enriched genes"

    cluster_enrichment_results = cluster_enrichment.set_index('native').sort_values('p_value').iloc[:,[2,5,7,10,1]]
    pd.set_option("display.max_colwidth", 800)
    return cluster_enrichment_results.iloc[:10,:]

def pathway_enrich_plot(genes, databases, title, background_genes, name_output, save: bool = False):
    """A function to plot the signature enrichment as a bargraph.  
    # Inputs:
    #    genes              - list of genes to be probed
    #    databases          - which databases to query, more information can be found here: https://biit.cs.ut.ee/gprofiler/page/apis
    #    title              - title for figure
    #    background_genes   - all the 
    #    save            - if I want to save the the figure
    # 
    """
    #Interpretation of differentially expressed genes in cluster 0 cells - g:profiler
    
    gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

    cluster_enrichment = gp.profile(organism='mmusculus', sources=databases, user_threshold=0.05,
                                   significance_threshold_method='fdr', 
                                   background= background_genes, 
                                   query=genes) #"contains the list of enriched genes"

    cluster_enrichment_results = cluster_enrichment.set_index('native').sort_values('p_value').iloc[:,[2,5,7,10,1]]

    # made new column with negative log p-value
    cluster_enrichment_results['-log10_p_value'] = cluster_enrichment_results['p_value'].map(lambda x: -math.log(x,10))
    
    if 'REAC:0000000' in cluster_enrichment_results.index.tolist():
        cluster_enrichment_results = cluster_enrichment_results.drop(labels='REAC:0000000', axis=0)

    plt.rcdefaults()
    fig, ax = plt.subplots()

    cluster_name = cluster_enrichment_results['name'].head(10)
    y_pos = np.arange(len(cluster_name))
    enrichment_value = cluster_enrichment_results['-log10_p_value'].head(10)

    ax.barh(y_pos, enrichment_value, align='center', color='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(cluster_name)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('-log10 p value')
    ax.set_title(title)

    if save:
        plt.savefig(name_output, format='pdf', bbox_inches = "tight")
    
    return plt.show()


# pathway enrichment for for genelist 'genes' probing 'databases', with genes that intersect with gene set
def pathway_enrich_genes(genes, databases):
    
    gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

    cluster_enrichment = gp.profile(organism='mmusculus', sources=databases, user_threshold=0.05,
                                   significance_threshold_method='fdr', 
                                   query=genes, #"contains the list of enriched genes"
                                   no_evidences=False)
    cluster_enrichment_results = cluster_enrichment.set_index('native').sort_values('p_value').iloc[:,[1,2,5,10,13]]

    pd.set_option("display.max_colwidth", 800)
    return cluster_enrichment_results.iloc[:10,:] 


# annoated list of genes ('genes') with genenames  - print the whole list
def gene_name_annotation_long(genes):
    gp = GProfiler(return_dataframe=True)
    gene_annot = gp.convert(organism='mmusculus',
             query= genes,
             target_namespace='ENTREZGENE_ACC')

    gene_annot['short_description'] = gene_annot['description'].map(lambda x: re.sub('\[.+\]', '', x)) # delete extra text between []
    gene_annot = gene_annot.drop(['description','name', 'converted','n_incoming','n_converted', 'namespaces', 'query'], axis=1)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # print all lines
        return display(gene_annot)
        
# annoated list of genes ('genes') with genenames  - print short list
def gene_name_annotation_short(genes):
    gp = GProfiler(return_dataframe=True)
    gene_annot = gp.convert(organism='mmusculus',
             query= genes,
             target_namespace='ENTREZGENE_ACC')

    gene_annot['short_description'] = gene_annot['description'].map(lambda x: re.sub('\[.+\]', '', x)) # delete extra text between []
    gene_annot = gene_annot.drop(['description','name', 'converted','n_incoming','n_converted', 'namespaces', 'query'], axis=1)
    return gene_annot


# calculate overlap between certain group - this drops all the false colums, so genes that are not present in 
def calc_overlap_old(samples, input_df, drop_false: bool = False):
    '''
    input is samples i.e. which HSC types
    input is raw df (csv)
    returns the DF with the multiindex - can be used with Upsetplot - 'output_df_all'
    returns the DF with the coefficients for only selected conditions - 'output_df_filt'
    '''
    output_df_filt = input_df.filter(regex='coef') # filter to only one column that has the LT subsets
    input_df_bin = output_df_filt.fillna(0) #replaces all NaNs with 0
    input_df_bin = (input_df_bin != 0) * 1 # sets all non zero values to one
    
    new_columns = ['HSC_0',
                 'HSC_1',
                 'HSC_2',
                 'HSC_3',
                 'HSC_4',
                 'HSC_5',
                 'HSC_all',
                 'MPP_surf',
                 'MPP1_surf',
                 'MPP2_surf',
                 'MPP3_surf',
                 'MPPs_0',
                 'MPPs_1',
                 'MPPs_2',
                 'MPPs_3',
                 'MPPs_4',
                 'MPPs_5',
                 'MPPs_6',
                 'MPPs_7',
                 'MPPs_all']
    #rename columns
    input_df_bin.columns = new_columns
    output_df_filt.columns = new_columns

    #select columns for only HSCs and MPPs
    input_df_bin = input_df_bin[samples]
    output_df_filt = output_df_filt[samples]

    #replace 0 and 1 with 'True' and 'False'
    input_df_bin.replace(0, False, inplace=True)
    input_df_bin.replace(1, True, inplace=True)

    #combine both df
    output_df_all = pd.concat([input_df_bin, input_df], axis=1) # combine thi primary
    output_df_filt = pd.concat([output_df_filt, input_df['primerid']], axis=1) # combine the filtered df

    #make multiindex
    output_df_all = output_df_all.set_index(list(input_df_bin.columns))
    
    # drop the genes that are absent in all of the datasets
    if drop_false:
        len_bool = (False, ) * len(samples) # generate a tuple with values that are False for the lenght of the sample list
        output_df_all = output_df_all.drop(output_df_all.loc[len_bool,:].index)
    
    return output_df_filt, output_df_all


# calculate overlap between certain group - this drops all the false colums, so genes that are not present in 
def calc_overlap(samples, input_df, drop_false: bool = False):
    '''
    input is samples i.e. which HSC types
    input is raw df (csv)
    returns the DF with the multiindex - can be used with Upsetplot - 'output_df_all'
    returns the DF with the coefficients for only selected conditions - 'output_df_filt'
    '''
    output_df_filt = input_df.filter(regex='coef') # filter to only one column that has the LT subsets
    input_df_bin = output_df_filt.fillna(0) #replaces all NaNs with 0
    input_df_bin = (input_df_bin != 0) * 1 # sets all non zero values to one
    
    new_columns = ['HSC_Metabo',
                 'HSC_Quiesc',
                 'HSC_Activated',
                 'HSC_Interferon',
                 'HSC_Acute-Activation',
                 'HSC_Cell-cycle',
                 'HSC_all',
                 'MPP',
                 'MPP1',
                 'MPP2',
                 'MPP3/4',
                 'LSK_Primed',
                 'LSK_Metabo',
                 'LSK_Prog',
                 'LSK_Cell-cycle',
                 'LSK_Acute-Activation',
                 'LSK_Interferon',
                 'LSK_Interferon cell-cycle',
                 'LSK_Myeloid',
                 'LSK_all']
    #rename columns
    input_df_bin.columns = new_columns
    output_df_filt.columns = new_columns

    #select columns for only HSCs and MPPs
    input_df_bin = input_df_bin[samples]
    output_df_filt = output_df_filt[samples]

    #replace 0 and 1 with 'True' and 'False'
    input_df_bin.replace(0, False, inplace=True)
    input_df_bin.replace(1, True, inplace=True)

    #combine both df
    output_df_all = pd.concat([input_df_bin, input_df], axis=1) # combine thi primary
    output_df_filt = pd.concat([output_df_filt, input_df['primerid']], axis=1) # combine the filtered df

    #make multiindex
    output_df_all = output_df_all.set_index(list(input_df_bin.columns))
    
    # drop the genes that are absent in all of the datasets
    if drop_false:
        len_bool = (False, ) * len(samples) # generate a tuple with values that are False for the lenght of the sample list
        output_df_all = output_df_all.drop(output_df_all.loc[len_bool,:].index)
    
    return output_df_filt, output_df_all


def average_expression(anndata, marker_dict, gene_symbol_key=None, partition_key='leiden'):
    """A function go get mean expressions of genes per cluster or group
    this method was adapted from here https://github.com/theislab/scanpy/issues/181
    removed the calculation of z-score since my data is scaled already anyways, also changed to extract raw data
    # 
    # Inputs:
    #    anndata         - An AnnData object containing the data set and a partition
    #    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
    #                      an anndata.var field with the key given by the gene_symbol_key input
    #    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
    #                      genes
    #    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
    #                      'leiden' """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        

    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_exp['specific'] = pd.Series({}, dtype='str')
    marker_names = []
    
    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        for gene in marker_dict[group]:
            ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
            if np.sum(ens_idx) == 0:
                continue
            else:
                anndata.obs[ens_idx[0]] = anndata.raw.X[:,ens_idx].mean(1) #works for both single and multiple mapping
                ens_idx = ens_idx[0]

            clust_marker_exp = anndata.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
            clust_marker_exp.append(group)
            marker_exp.loc[i] = clust_marker_exp
            marker_names.append(gene)
            i+=1

    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

# get genes that are expressed at least in one subtype of LTs
def make_categories_HSC(row):
    if (row['coef'] > 0) | (row['coefLT_1'] > 0) | (row['coefLT_2'] > 0) | (row['coefLT_3'] > 0) | (row['coefLT_5'] > 0) | (row['coefLT_all'] > 0):
        return 1
    if (row['coef'] < 0) | (row['coefLT_1'] <  0) | (row['coefLT_2'] <  0) | (row['coefLT_3'] <  0) | (row['coefLT_5'] <  0) | (row['coefLT_all'] <  0):
        return -1 

# get genes that are expressed at least in one subtype of MPPs
def make_categories_MPP(row):
    if (row['coefMPPs_0'] > 0) | (row['coefMPPs_1'] > 0) | (row['coefMPPs_2'] > 0) | (row['coefMPPs_3'] > 0) | (row['coefMPPs_all'] > 0):
        return 1
    if (row['coefMPPs_0'] < 0) | (row['coefMPPs_1'] <  0) | (row['coefMPPs_2'] <  0) | (row['coefMPPs_3'] <  0) | (row['coefMPPs_all'] <  0):
        return -1 
    
    
# get genes that are expressed at least in one subtype of LTs in male
def make_categories_HSC_m(row):
    if (row['coef'] > 0) | (row['coefLT_1male'] > 0) | (row['coefLT_2male'] > 0) | (row['coefLT_3male'] > 0) | (row['coefLT_5male'] > 0) | (row['coefLT_allmale'] > 0):
        return 1
    if (row['coef'] < 0) | (row['coefLT_1male'] <  0) | (row['coefLT_2male'] <  0) | (row['coefLT_3male'] <  0) | (row['coefLT_5male'] <  0) | (row['coefLT_allmale'] <  0):
        return -1 

# get genes that are expressed at least in one subtype of LTs in female
def make_categories_HSC_f(row):
    if (row['coefLT_0female'] > 0) | (row['coefLT_1female'] > 0) | (row['coefLT_2female'] > 0) | (row['coefLT_3female'] > 0) | (row['coefLT_5female'] > 0) | (row['coefLT_allfemale'] > 0):
        return 1
    if (row['coefLT_0female'] < 0) | (row['coefLT_1female'] <  0) | (row['coefLT_2female'] <  0) | (row['coefLT_3female'] <  0) | (row['coefLT_5female'] <  0) | (row['coefLT_allfemale'] <  0):
        return -1 

    
# get genes that are expressed at least in one subtype of MPPs in male
def make_categories_MPP_m(row):
    if (row['coefMPP_0male'] > 0) | (row['coefMPP_1male'] > 0) | (row['coefMPP_2male'] > 0) | (row['coefMPP_3male'] > 0) | (row['coefMPP_allmale'] > 0):
        return 1
    if (row['coefMPP_0male'] < 0) | (row['coefMPP_1male'] <  0) | (row['coefMPP_2male'] <  0) | (row['coefMPP_3male'] <  0) | (row['coefMPP_allmale'] <  0):
        return -1 

# get genes that are expressed at least in one subtype of MPPs in female
def make_categories_MPP_f(row):
    if (row['coefMPP_0female'] > 0) | (row['coefMPP_1female'] > 0) | (row['coefMPP_2female'] > 0) | (row['coefMPP_3female'] > 0) | (row['coefMPP_allfemale'] > 0):
        return 1
    if (row['coefMPP_0female'] < 0) | (row['coefMPP_1female'] <  0) | (row['coefMPP_2female'] <  0) | (row['coefMPP_3female'] <  0) | (row['coefMPP_allfemale'] <  0):
        return -1 

    
    
# define the overlap column
def make_overlap_column(row):
    if (row['HSC_any'] != 0) & (row['MPP_any'] != 0) & (row['overlap'] != 0):
        return 'overlap'
    if (row['HSC_any'] != 0) & (row['MPP_any'] == 0) & (row['overlap'] == 0):
        return 'HSC_only'
    if (row['HSC_any'] == 0) & (row['MPP_any'] != 0) & (row['overlap'] == 0):
        return 'MPP_only'
    
    
# define the overlap column by sex
def make_overlap_column_HSC(row):
    if (row['HSC_any_m'] != 0) & (row['HSC_any_f'] != 0) & (row['overlap_HSC'] != 0):
        return 'overlap_HSC'
    if (row['HSC_any_m'] != 0) & (row['HSC_any_f'] == 0) & (row['overlap_HSC'] == 0):
        return 'male_only_HSC'
    if (row['HSC_any_m'] == 0) & (row['HSC_any_f'] != 0) & (row['overlap_HSC'] == 0):
        return 'female_only_HSC'
    
    
def make_overlap_column_MPP(row):
    if (row['MPP_any_m'] != 0) & (row['MPP_any_f'] != 0) & (row['overlap_MPP'] != 0):
        return 'overlap_MPP'
    if (row['MPP_any_m'] != 0) & (row['MPP_any_f'] == 0) & (row['overlap_MPP'] == 0):
        return 'male_only_MPP'
    if (row['MPP_any_m'] == 0) & (row['MPP_any_f'] != 0) & (row['overlap_MPP'] == 0):
        return 'female_only_MPP'