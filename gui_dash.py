#Author: Tess Marvin (tmarvin@nd.edu)
#Usage: pythonw gui_dash.py
#This function can take from 1 to 7 csv files and create one scatter plot of transcription over a time course
#It would be useful if your file naming convention is malariaisolatelog2anythinghere.csv
#e.g. NHP4026log2imputedtimecourse.csv
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from argparse import ArgumentParser
import numpy as np
from gooey import Gooey
from gooey import GooeyParser
def gene_fun(gene, csv_files):
    #first check to make sure the gene is in every CSV file provided
    #if the gene name has been misspelled, then the user will be notified
    for file in csv_files:
        data_c = pd.read_csv(file)
        if gene not in data_c:
            print('Please enter correct gene name')
            return(None)
#This is the plotting function -- must be abstract to accept up to 7 csv files at once
    color_list = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    num = 0

    been_thru = False
    for file in csv_files:
        date = pd.read_csv(file)
        #this is how to splice the malaria isolate name from the csv file name
        top, bot = os.path.split(file)
        isolate, b =bot.split('log2')
        if not been_thru:
            pl = date.plot(kind='scatter', x='Time', y=gene, color=color_list[num], figsize = [7,7], title = gene, label = isolate)
            num += 1
            been_thru = True
        else:
            name = date.plot(kind='scatter', x='Time', y=gene, color=color_list[num], label = isolate, ax = pl)
            name.set(xlabel = 'Time (hrs)', ylabel= 'Log2 (Fold Change)')
            num += 1
    plt.show()

def sample_and_geneqc(expdata):
    '''
    Here, we'll remove poor samples and poor genes, and QC low read #s
    Low reads are any gene/sample reading < 5 counts; those are zeroed out
    Poor samples are defined as those with < 3000 genes with reads
    Poor genes are those that appear in < 20% of samples
    '''

    # Curate Samples
    #If there are less than 5 reads -- zero it out
    expdata = expdata.mask(expdata < 5, 0)
    #so np.count_nonzero with axis = 1 counts the number of nonzero counts for a gene across all columns (samples)
    #basically, this counts how many samples this gene is in
    genecounts = pd.Series(data = np.count_nonzero(expdata, axis = 1),
                             index = expdata.index)
    #so np.count_nonzero with axis = 0 counts the number of nonzero counts for a sample across all rows (genes)
    #basically this counts how many genes with reads each sample has
    samplecounts = pd.Series(np.count_nonzero(expdata,axis = 0),
                          index = expdata.columns)
    #this ensures that only the genes that are in more than 20% of samples are included
    goodgenes = genecounts[genecounts/samplecounts.size > 0.2]

    goodsamples = samplecounts[samplecounts > 3000]
    #the .loc function creates a dataframe that only contains the goodgenes and goodsamples
    allcur = expdata.loc[goodgenes.index, goodsamples.index]

    return allcur

def tpm_norm(expdata, probefile, gctdata, fraglength = 250):
    '''
    TPM is (reads / effective gene length (in kb)) / (sample reads / 1e6)
    The elaborate nonsense below just calculates this, it's a bit of a
    mess but hey that's programming
    '''

    # Read in probe metadata, subtract frag length,
    # and convert to dict


    probeinfo = pd.read_csv(probefile, usecols = [1,6])
    probeinfo.index = probeinfo['Geneid']
    probeinfo.drop(columns = 'Geneid', inplace = True)
    probeinfo = probeinfo[probeinfo.index.isin(gctdata.index)]
    probeinfo['Length'] = probeinfo['Length'].apply(lambda x: x-fraglength)

    # Build output frame from input frame (cheating), calc length and
    # build lookup from samplename to total counts

    tpm = expdata.copy()
    curated_tpm = expdata.copy()

    # Iterate over rows- if non sample rows, just copy data, if sample
    # rows, calculate TPM

    # Note: it really doesn't like the way I did this and will throw
    # warnings- it does work, it's fine, don't worry about it, if you
    # are smarter than I am go ahead and rewrite it

    for (col, data) in expdata.iteritems():

        numreads = data.sum()

        tempframe = pd.concat([probeinfo, data], axis = 1)

        tempframe['rpkb'] = tempframe.iloc[:,1].divide(tempframe['Length'].divide(1000))
        tempframe['tpm'] = (tempframe['rpkb'] / data.sum()) * 1e6
        tempframe['tpm'] = tempframe['tpm'].clip(lower = 0)
        curated_tpm.loc[:,col] = tempframe['tpm']
    # Removing genes whose expression is ~0 in > 80% of samples
    goodgenes=[]
    for index, row in curated_tpm.iterrows():

        if (np.count_nonzero(row > 0.1) / len(row) > 0.2):
            goodgenes.append(index)

    curated_tpm = curated_tpm[curated_tpm.index.isin(goodgenes)]
    return curated_tpm

def fix_names(countdata):
    countdata = countdata.T
    '''
    This fixes the name specific issues in the NF54GFPxNHP4026 cross.
    The names have changed several times and been recorded in different
    formats, so I'll fix them with this.
    '''
    countdata.index = countdata.index.str.replace('\/','', regex = True)
    countdata.index = countdata.index.str.replace('ND5A5', 'AC075', regex = True)
    countdata.index = countdata.index.str.replace('ND6G8', 'AC125', regex = True)
    countdata.index = countdata.index.str.replace('N1', '', regex = True)
    countdata.index = countdata.index.str.replace('\\.', '', regex = True)
    countdata.index = countdata.index.str.replace('_4026', '_NHP4026', regex = True)
    countdata.index = countdata.index.str.replace('^4026', 'NHP4026', regex = True)
    countdata.index = countdata.index.str.replace('2H9', 'AC030', regex = True)
    countdata.index = countdata.index.str.replace('6E5', 'AC033', regex = True)
     # This is a hot mess of regular expression that just converts the entire
    # sample name to strain_##, where ## is the sampling timepoint. The formatting
    # is very inconsistent throughout, so a mess of replacements need to be carefully
    # made. This mess of substitutions makes them.
    countdata.index = countdata.index.str.replace('^GF_PL[\d]+[a,b,c]{0,1}_', '',regex = True)
    countdata.index = countdata.index.str.replace('[A,B,C]_', '', regex = True)
    countdata.index = countdata.index.str.replace('_[d]+$', '', regex = True)
    countdata.index = countdata.index.str.replace('_S.*', '', regex = True)
    countdata.index = countdata.index.str.replace('_[0-9]{3,4}$', '', regex = True)
    countdata.index = countdata.index.str.replace('hpi', '', regex = True)
    countdata.index = countdata.index.str.replace('_T', '_', regex = True)
    return countdata.T
def parentalboxes(gene, df_correctnames):
    '''
    So this function makes box plots for the parental parasites at the three time points (4, 30, 44 hrs)
    the parents are: NF54gfp x NHP4026
    The data can be accessed like this normalized_fixed_names.loc['gene', 'rep']
    e.g. normalized_fixed_names.loc['PF3D7_0100100', 'NF54gfp_44']
    '''
    timepoints = ['4', '30', '44']
    NF54data= []
    NHP4026data = []
    for time in timepoints:
        sample_name_NF = str('NF54gfp_' + time)
        l = list(df_correctnames.loc[gene, sample_name_NF])
        NF54data.append(l)
        sample_name_NH = str('NHP4026_' + time)
        m = list(df_correctnames.loc[gene, sample_name_NH])
        NHP4026data.append(m)
    data = {}
    data['4 hpi'] = {}
    data['30 hpi'] = {}
    data['44 hpi'] = {}
    num = 0
    for k,v in data.items():
        v['NF54gfp']= NF54data[num]
        v['NHP4026']= NHP4026data[num]
        num = num + 1
    fig, axes = plt.subplots(ncols=3, sharey=True)
    fig.subplots_adjust(wspace=0)
    for ax, name in zip(axes, ['4 hpi', '30 hpi', '44 hpi']):
        ax.boxplot([data[name][item] for item in ['NF54gfp', 'NHP4026']])
        ax.set(xticklabels=['NF54gfp', 'NHP4026'], xlabel=name)
        if name == '4 hpi':
            ax.set(ylabel='Gene Expression (TPM)')
        if name == '30 hpi':
            ax.set(title=gene)
        ax.margins(0.05)
    plt.show()

def progenybar(gene, df_correctnames):
    '''
    So this function makes a scatter plot depicting the gene expression for a gene of interst across all progeny (each replicate is a point)
    This is specifically for data at time points 4, 30, and 44
    '''
    how_many_reps = {}
    names4 = {}
    names30 = {}
    names44 = {}
    for name in df_correctnames.loc[gene].index:
        if name not in how_many_reps.keys():
            how_many_reps[name] = 1
        else:
            how_many_reps[name] += 1
    for name in df_correctnames.loc[gene].index:
        name_str = str(name)
        if 'NF54gfp' in name_str:
            continue
        if 'NHP4026' in name_str:
            continue
        if '_44' in name_str:
            if name_str in names44.keys():
                continue
            else:
                if how_many_reps[name] == 1:
                    m = []
                    m.append(df_correctnames.loc[gene, name_str])
                    names44[name] = m
                else:
                    m = list(df_correctnames.loc[gene, name_str])
                    names44[name] = m
        if '_30' in name:
            if name_str in names30.keys():
                continue
            else:
                if how_many_reps[name] == 1:
                    m = []
                    m.append(df_correctnames.loc[gene, name_str])
                    names30[name] = m
                else:
                    m = list(df_correctnames.loc[gene, name_str])
                    names30[name] = m
        if '_4' in name:
            if '_44' in name:
                continue
            if name_str in names4.keys():
                continue
            else:
                if how_many_reps[name] == 1:
                    m = []
                    m.append(df_correctnames.loc[gene, name_str])
                    names4[name] = m
                else:
                    m = list(df_correctnames.loc[gene, name_str])
                    names4[name] = m
    data = {}
    data['4 hpi'] = names4
    data['30 hpi'] = names30
    data['44 hpi'] = names44
    fig, axs = plt.subplots(3,1)
    fig.subplots_adjust(hspace=0.5)
    fig.suptitle(gene)
    for ax, name in zip(axs, ['4 hpi', '30 hpi', '44 hpi']):
        for k in data[name]:
            ax.scatter([k]*len(data[name][k]), data[name][k])
            samp = 'Samples at ' + str(name)
            ax.set(xlabel= samp, ylabel='Gene Expression (TPM)')
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize = 6)
    plt.show()

@Gooey(
    program_name='Plasmodium falciparum Genomic Analysis',
    menu=[{
    'name':'File',
    'items': [{
            'type': 'AboutDialog',
            'menuTitle': 'About',
            'name': 'Genomic Analysis of Plasmodium falciparum',
            'description': 'A tool to probe drug resistance in malaria',
            'version': '1.0',
            'copyright': '2020',
            'website': 'https://github.com/TessMarvin',
            'developer': 'Tess Marvin (tmarvin@nd.edu)',
            'license': 'University of Notre Dame'
    }]
    }]
)
def main():
    #So first we will handle the arguments that are "required" -- the files and the gene of interest
    #Here we give our GUI a title
    parser = GooeyParser(description="Dashboard for Gene Expression Analysis of Plasmodium falciparum")
    #Here we allow for the entry of the gene of interest's name
    parser.add_argument('genename', help='Enter the gene you wish to analyze.')
    #Next, we have a section where the user can select how they would like to analyze the data
    data_processing = parser.add_argument_group("Data Processing Options", "Customize Your Analysis")
    #First way to analyze the data: a scatterplot of the gene expression over time
    data_processing.add_argument('-scat', help='Enable Scatter Plot', action='store_true', widget='BlockCheckbox')
    #Here we allow for the selection of the data files to analyze
    data_processing.add_argument("-file_chooser", nargs='*', help = 'Choose Files to include in scatterplot analysis.', widget='MultiFileChooser')
    #Second way to analyze the data: parental box subplots
    data_processing.add_argument('-box', help='Enable Box Plots for Parental Replicates', action='store_true', widget='BlockCheckbox')
    #Third way to analyze data: progeny scatter plot of gene expression
    data_processing.add_argument('-progeny', help='Enable Scatter Plots for Progeny Replicates', action='store_true', widget='BlockCheckbox')
    #Now we need a file chooser for the probe meta data
    data_processing.add_argument("-probe_chooser", help = 'Choose Probe MetadataFile to include in progeny/parental analysis.', widget='FileChooser')
    #We also need a file chooser to pick out the raw count data
    data_processing.add_argument("-count_data", help = 'Choose Raw Count Data to include in progeny/parental analysis.', widget='FileChooser')

    #Now we parse all these arguments
    args = parser.parse_args()
    #Save the entry of the gene name
    gene = args.genename
    #Save whether or not the user would like to produce a scatter plot
    scatter= args.scat
    #Save the file pathways that the user has selected
    csvfiles= args.file_chooser
    #Save whether or not the user would like to produce a parental box plot
    boxer = args.box
    #Save whether or not the user would like to produce a progeny scatter plot
    prog = args.progeny
    #Save the probe MetadataFile
    probefile = args.probe_chooser
    #Save the raw count data
    raw_counts_data = args.count_data
    #make raw_couts_data a dataframe
    gctdata = pd.read_csv(raw_counts_data, sep = '\t', index_col = 0)

    #If the user would like a scatter plot, produce one
    if(scatter):
        #if no CSV files are found and we want to scatter, fail gracefully and request the CSV files be provided
        if(len(csvfiles) == 0):
            print("Please ensure that you select files to analyze")
            return(None)
        else:
            gene_fun(gene, csvfiles)
    if(boxer):
        #Quality control the reads
        qc_counts_data= sample_and_geneqc(gctdata)
        #Normalize the QC data
        normalized_counts= tpm_norm(qc_counts_data, probefile, gctdata)
        #Fix the naming conventions
        normalized_fixed_names = fix_names(normalized_counts)
        #plot the box
        parentalboxes(gene, normalized_fixed_names)
    if(prog):
        #Quality control the reads
        qc_counts_data= sample_and_geneqc(gctdata)
        #Normalize the QC data
        normalized_counts= tpm_norm(qc_counts_data, probefile, gctdata)
        #Fix the naming conventions
        normalized_fixed_names = fix_names(normalized_counts)
        #plot the box
        progenybar(gene, normalized_fixed_names)
if __name__ == '__main__':
   main()
