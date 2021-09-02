import pandas as pd
from gprofiler import GProfiler
import re
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors

# for making custom colormaps
def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


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


# pathway enrichment for for genelist 'genes' probing 'databases', with genes that intersect with gene set
def pathway_enrich_genes_new(genes, databases):
    
    gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

    cluster_enrichment = gp.profile(organism='mmusculus', sources=databases, user_threshold=0.05,
                                   significance_threshold_method='fdr', 
                                   query=genes, #"contains the list of enriched genes"
                                   no_evidences=False)
    cluster_enrichment_results = cluster_enrichment.set_index('native').sort_values('p_value').iloc[:,[1,2,5,7,10,13]]

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
    #                      'leiden' 
    """

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
    
    
def multiple_testing_correction(ps, alpha=0.05, method='benjamini-hochberg', **kwargs):
    """ correct pvalues for multiple testing and add corrected `q` value
    :param ps: list of pvalues
    :param alpha: significance level default : 0.05
    :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
    :returns (q, rej): two lists of q-values and rejected nodes
    """
    # Implement copy from GOATools.
    _p = np.array(ps)
    q = _p.copy()
    rej = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    if method == 'bonferroni':
        q[mask] = p * len(p)
        rej[mask] = q[mask] < alpha
    elif method == 'benjamini-hochberg':
        _rej, _q = fdrcorrection(p, alpha)
        rej[mask] = _rej
        q[mask] = _q
    else:
        raise ValueError(method)
    return q, rej


def fdrcorrection(pvals, alpha=0.05):
    """ benjamini hocheberg fdr correction. inspired by statsmodels 
    """
    # Implement copy from GOATools.
    pvals = np.asarray(pvals)
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)

    ecdffactor = _ecdf(pvals_sorted)
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True
    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    del pvals_corrected
    reject_ = np.empty_like(reject)
    reject_[pvals_sortind] = reject
    return reject_, pvals_corrected_

def _ecdf(x):
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)


def enrichr(cluster_df, cluster_list, ref_gene_set, cutoff_adj_pval = 0.01, cutoff_logfch = 0.58, number_of_genes = 100, bckgr_genes = 14000):
    """A function to test enrichment of several genesets.  
    # Inputs:
    #    cluster_df            - cluster enrichment dataframe
    #    cluster_list          - list with names of the clusters
    #    ref_gene_set          - list or dictionary of reference gene sets
    #    cutoff_adj_pval       - cutoff for the adjusted p-value
    #    cutoff_logfch         - cutoff for the log fold change (has to be larger), default is 1.5 fold
    #    number_of_genes       - number of how many genes to test enrichment with
    # 
    """
    
    gsea_res = dict()
    pred = dict()

    for cl in cluster_list:
        
        #filter by only positive values + FDR < 0.01
        column_pvalue = cl + '_p'
        column_fch = cl + '_l'
        column_name = cl + '_n'
        sort_column = cl + '_s'
        
        glist = cluster_df[(cluster_df[column_pvalue] < cutoff_adj_pval) & (cluster_df[column_fch] > cutoff_logfch)][column_name].tolist()
       
        enr_res = gp.enrichr(gene_list=glist[:number_of_genes],
                         organism='Mouse',
                         gene_sets= ref_gene_set,
                         background = bckgr_genes,    
                         description='pathway',
                         cutoff = 0.05)
        if enr_res.results.shape[0] == 0:
            pred[cl] = "Unass"
        
        else:
            enr_res.results.sort_values(by="P-value",axis=0, ascending=True, inplace=True)
            gsea_res[cl]=enr_res
            pred[cl] = enr_res.results["Term"][0]
    return gsea_res

def enrichr_hclust(cluster_dic, ref_gene_set):
    """A function to test enrichment of several genesets.  
    # Inputs:
    #    cluster_dic            - cluster enrichment dataframe
    #    cluster_list           - list with names of the clusters
    #    ref_gene_set           - list or dictionary of reference gene sets
    #    cutoff_adj_pval        - cutoff for the adjusted p-value
    #    cutoff_logfch          - cutoff for the log fold change (has to be larger), default is 1.5 fold
    #    number_of_genes        - number of how many genes to test enrichment with
    # 
    """
    
    gsea_res = dict()
    pred = dict()

    for cl in list(clust_dic.keys()):

        enr_res = gp.enrichr(gene_list=cluster_dic[cl],
                         organism='Mouse',
                         gene_sets= ref_gene_set,
                         #background = bckgr_genes,    
                         description='pathway',
                         cutoff = 0.05)
        if enr_res.results.shape[0] == 0:
            pred[cl] = "Unass"
        
        else:
            enr_res.results.sort_values(by="P-value",axis=0, ascending=True, inplace=True)
            gsea_res[cl]=enr_res
            pred[cl] = enr_res.results["Term"][0]
    return gsea_res

def ribofilter(genelist):
    filtered_genelist = [item for item in genelist if 'Rpl' not in item]
    filtered_genelist = [item for item in filtered_genelist if 'Rps' not in item]
    return filtered_genelist

def col_select(temp_df):
    temp_df['Gene Overlap'] = temp_df['intersection_size'].astype(str) + '/' + temp_df['term_size'].astype(str)
    temp_df['Adjusted P-value'] = temp_df['p_value']
    temp_df['Pathway'] = temp_df['name']
    temp_df['Origin'] = temp_df.index
    temp_df['Gene names'] = temp_df['intersections']
    temp_df = temp_df[['Origin', 'Pathway', 'Gene Overlap', 'Adjusted P-value', 'Gene names']]
    return temp_df


def render_mpl_table(data, row_height=10, font_size=12,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        #size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=(7,8))
        ax.axis('off')
    mpl_table = ax.table(cellText=data.values, colLabels=data.columns, loc='center', cellLoc='center', **kwargs)
    mpl_table.auto_set_column_width(col=list(range(len(data.columns)))) # Provide integer list of columns to adjust
    mpl_table.auto_set_font_size(True)
    mpl_table.set_fontsize(font_size)
    cellDict=mpl_table.get_celld()
    cellDict[(0,0)].set_width(0.5)

    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return ax.get_figure(), ax