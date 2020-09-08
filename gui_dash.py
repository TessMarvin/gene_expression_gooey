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
#the first argument is the gene name
#next, the program searches recurisively thru the current directory for .csv files
#or the files can be passed in by stdin as part of a pipeline
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
    #Here we allow for the selection of the data files to analyze
    parser.add_argument("file_chooser", nargs='*', help = 'Choose Files to include in analysis.', widget='MultiFileChooser')
    #Next, we have a section where the user can select how they would like to analyze the data
    data_processing = parser.add_argument_group("Data Processing Options", "Customize Your Analysis")
    #At the moment, there is only one way to analyze the data: a scatterplot of the gene expression over time
    data_processing.add_argument('-scat', help='Enable Scatter Plot', action='store_true', widget='BlockCheckbox')
    #Now we parse all these arguments
    args = parser.parse_args()
    #Save the entry of the gene name
    gene = args.genename
    #Save whether or not the user would like to produce a scatter plot
    scatter= args.scat
    #Save the file pathways that the user has selected
    csvfiles= args.file_chooser
    #if no CSV files are found, fail gracefully and request the CSV files be in working directory
    if(len(csvfiles) == 0):
        print("Please ensure that you select files to analyze")
        return(None)
    else:
        #If the user would like a scatter plot, produce one
        if(scatter):
            gene_fun(gene, csvfiles)

if __name__ == '__main__':
   main()
