import pandas as pd
import os,sys,subprocess,re, argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
#import data
parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input file')
parser.add_argument('-a', '--annot', help='annotation file')
parser.add_argument('-c', '--categories', help='NOG categories you would like to include')

args=parser.parse_args()
data=args.input
annot=args.annot
cat= args.categories

df= pd.read_csv(data)
df.head()
annotation= pd.read_csv(annot, sep='\t')

#map annotation to df
#split contigs
df['Contig']= df['Contig'].str.split(' ').str[0]    
df['Nog_Category']=df['Contig'].map(annotation.set_index('query')['nog_categories'])
df['Function_GVOG']=df['Contig'].map(annotation.set_index('query')['function'])
df['Nog_desc']=df['Contig'].map(annotation.set_index('query')['nog_desc'])
df['Nog_desc']= df['Nog_desc'].str.split('\|').str[0]
#split nog categories
df['Nog_Category']= df['Nog_Category'].str.split(' ').str[0]
#adding in nog desc
for i in range(len(df)):
    if df['Function_GVOG'][i]== 'no_annotation':
        df['Function_GVOG'][i]= df['Nog_desc'][i]
#split at |
df['Function_GVOG']= df['Function_GVOG'].str.split('\|').str[0]
        #droping some cols
df.to_csv(re.sub('.csv','_annotated.csv',data),index=False)
df['Nog_Category'].value_counts()
#write the nog categories to a file
df['Nog_Category'].value_counts().to_csv(re.sub('.csv','_nog_categories.csv',data),index=True)
#plotting
df2= df.drop('Contig',axis=1)
df2=df2.drop('Nog_Category',axis=1)
df2=df2.drop('Function_GVOG',axis=1)
df2=df2.drop('Nog_desc',axis=1)
#log transform
#make df2 numeric
df2= df2.apply(pd.to_numeric)
df2= np.log2(df2+1)
sns.clustermap(df2, col_cluster=False, cmap='PuBu', figsize=(10, 10))
#save figure
plt.savefig(re.sub('.csv','_clustermap_total.png',data))
# print the top 20 most abundant gene means only taking means if column is numeric
df['mean_rpkm'] = df.iloc[:,1:41].mean(axis=1)
top=df['mean_rpkm'].sort_values(ascending=False).head(50)
#make data frame of top 50 including mean rpkm and Function_GVOG
top=pd.DataFrame(top)
top=top.merge(df[['Function_GVOG']], left_index=True, right_index=True)
top.to_csv(re.sub('.csv','_top50.csv',data))

#transpose the data frame
df2=df2.transpose()
#sum the rows
df2['sum_rpkm']=df2.sum(axis=1)
#make first column not the index
df2=df2.reset_index()
#make a new data frame with the sum_rpkm and date
sums=pd.DataFrame(df2[['sum_rpkm','index']])
#sort by date
sums=sums.sort_values(by='index')
#change index to datetime
sums['index']=pd.to_datetime(sums['index'])
#plot the sums and make x axis marks every 10 days
#make this plot a variable so that it can be saved

sums.plot(x='index', y='sum_rpkm', figsize=(10, 10))
plt.xticks(rotation=45)
plt.xlabel('Date')
plt.ylabel('Sum of RPKM')
plt.title('Sum of RPKM over time')
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=5))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
#save the plot
plt.savefig(re.sub('.csv','_sums.png',data))


#Graphing by category
df3=df[df['Nog_Category']=='L']
#if df3 is empty then skip
if 'L' in cat:
        g= sns.clustermap(df3.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df3['Function_GVOG'])
        g.ax_heatmap.set_title('Replication and Repair')
        plt.savefig(re.sub('.csv','_replication_repair.png',data))
        
#add a title

df_K=df[df['Nog_Category']=='K']
if 'K' in cat:
        df_K.head()
        g= sns.clustermap(df_K.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_K['Function_GVOG'])
        g.ax_heatmap.set_title('Transcription')
        plt.savefig(re.sub('.csv','_transcription.png',data))

df_J=df[df['Nog_Category']=='J']
if 'J' in cat:
        g= sns.clustermap(df_J.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_J['Function_GVOG'])
        g.ax_heatmap.set_title('Translation')
        plt.savefig(re.sub('.csv','_translation.png',data))

df_M=df[df['Nog_Category']=='M']
if 'M' in cat:
        g= sns.clustermap(df_M.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_M['Function_GVOG'])
        g.ax_heatmap.set_title('Cell membrane biogenesis')
        plt.savefig(re.sub('.csv','_cell_membrane.png',data))

df_H=df[df['Nog_Category']=='H']
if 'H' in cat:
        g= sns.clustermap(df_H.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_H['Function_GVOG'])
        g.ax_heatmap.set_title('Coenzyme metabolism')
        plt.savefig(re.sub('.csv','_coenzyme.png',data))

df_Q=df[df['Nog_Category']=='Q']
if 'Q' in cat:
        g= sns.clustermap(df_Q.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_Q['Function_GVOG'])
        g.ax_heatmap.set_title('Secondary Metabolites')
        plt.savefig(re.sub('.csv','_secondary_metabolites.png',data))

df_A=df[df['Nog_Category']=='A']
if 'A' in cat:
        df_A['Function_GVOG']=df_A['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_A.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_A['Function_GVOG'])
        g.ax_heatmap.set_title('Rna Processing')
        plt.savefig(re.sub('.csv','_rna_processing.png',data))

df_O=df[df['Nog_Category']=='O']
if 'O' in cat:
        g= sns.clustermap(df_O.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_O['Function_GVOG'])
        g.ax_heatmap.set_title('Post Translational Modification')
        plt.savefig(re.sub('.csv','_post_trans_modification.png',data))

df_E=df[df['Nog_Category']=='E']
if 'E' in cat:
        df_E['Function_GVOG']=df_E['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_E.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_E['Function_GVOG'])
        g.ax_heatmap.set_title('Amino Acid Metabolism')
        plt.savefig(re.sub('.csv','_aa_metabolism.png',data))

df_C=df[df['Nog_Category']=='C']
if 'C' in cat:
        df_C['Function_GVOG']=df_C['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_C.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_C['Function_GVOG'])
        g.ax_heatmap.set_title('Energy Production')
        plt.savefig(re.sub('.csv','_energy_production.png',data))

df_B=df[df['Nog_Category']=='B']
if 'B' in cat:
        df_B['Function_GVOG']=df_B['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_B.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_B['Function_GVOG'])
        g.ax_heatmap.set_title('Chromatin Structure')
        plt.savefig(re.sub('.csv','_chromatin.png',data))

df_D=df[df['Nog_Category']=='D']
if 'D' in cat:
        df_D['Function_GVOG']=df_D['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_D.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_D['Function_GVOG'])
        g.ax_heatmap.set_title('Cell Cycle')
        plt.savefig(re.sub('.csv','_cell_cycle.png',data))

df_N=df[df['Nog_Category']=='N']
if 'N' in cat:
        df_N['Function_GVOG']=df_N['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_N.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_N['Function_GVOG'])
        g.ax_heatmap.set_title('Cell Motility')
        plt.savefig(re.sub('.csv','_cell_motility.png',data))

df_T=df[df['Nog_Category']=='T']
if 'T' in cat:
        df_T['Function_GVOG']=df_T['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_T.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_T['Function_GVOG'])
        g.ax_heatmap.set_title('Signal Transduction')
        plt.savefig(re.sub('.csv','_signal_transduction.png',data))

df_U=df[df['Nog_Category']=='U']
if 'U' in cat:
        df_U['Function_GVOG']=df_U['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_U.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_U['Function_GVOG'])
        g.ax_heatmap.set_title('Intracellular Trafficking')
        plt.savefig(re.sub('.csv','_intracellular_trafficking.png',data))

df_V=df[df['Nog_Category']=='V']
if 'V' in cat:
        df_V['Function_GVOG']=df_V['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_V.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_V['Function_GVOG'])
        g.ax_heatmap.set_title('Defense Response')
        plt.savefig(re.sub('.csv','_defense_response.png',data))

df_W=df[df['Nog_Category']=='W']
if 'W' in cat:
        df_W['Function_GVOG']=df_W['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_W.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_W['Function_GVOG'])
        g.ax_heatmap.set_title('Extracellular Structure')
        plt.savefig(re.sub('.csv','_extracellular_structure.png',data))

df_Y=df[df['Nog_Category']=='Y']
if 'Y' in cat:
        df_Y['Function_GVOG']=df_Y['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_Y.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_Y['Function_GVOG'])
        g.ax_heatmap.set_title('Nuclear Structure')
        plt.savefig(re.sub('.csv','_nuclear_structure.png',data))

df_Z=df[df['Nog_Category']=='Z']
if 'Z' in cat:
        df_Z['Function_GVOG']=df_Z['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_Z.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_Z['Function_GVOG'])
        g.ax_heatmap.set_title('Cytoskeleton')
        plt.savefig(re.sub('.csv','_cytoskeleton.png',data))

df_F=df[df['Nog_Category']=='F']
if 'F' in cat:
        df_F['Function_GVOG']=df_F['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_F.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_F['Function_GVOG'])
        g.ax_heatmap.set_title('Nucleotide Transport')
        plt.savefig(re.sub('.csv','_nucleotide_transport.png',data))

df_G=df[df['Nog_Category']=='G']
if 'G' in cat:
        df_G['Function_GVOG']=df_G['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_G.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_G['Function_GVOG'])
        g.ax_heatmap.set_title('Carbohydrate Metabolism')
        plt.savefig(re.sub('.csv','_carbohydrate_metabolism.png',data))

df_I=df[df['Nog_Category']=='I']
if 'I' in cat:
        df_I['Function_GVOG']=df_I['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_I.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_I['Function_GVOG'])
        g.ax_heatmap.set_title('Lipid Metabolism')
        plt.savefig(re.sub('.csv','_lipid_metabolism.png',data))

df_P=df[df['Nog_Category']=='P']
if 'P' in cat:
        df_P['Function_GVOG']=df_P['Function_GVOG'].str.split('|').str[0]
        g= sns.clustermap(df_P.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_P['Function_GVOG'])
        g.ax_heatmap.set_title('Inorganic Ion Transport')
        plt.savefig(re.sub('.csv','_inorganic_ion_transport.png',data))
