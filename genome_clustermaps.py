import pandas as pd
import os,sys,subprocess,re, argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
#import data
data=sys.argv[1]
df= pd.read_csv(data)
df.head()
annotation= pd.read_csv(sys.argv[2], sep='\t')
annotation.head()

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


#Graphing by category
df3=df[df['Nog_Category']=='L']
df3.head()
g= sns.clustermap(df3.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df3['Function_GVOG'])
g.ax_heatmap.set_title('Replication and Repair')
plt.savefig(re.sub('.csv','replication_repair.png',data))
#add a title

df_K=df[df['Nog_Category']=='K']
df_K.head()
g= sns.clustermap(df_K.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_K['Function_GVOG'])
g.ax_heatmap.set_title('Transcription')
plt.savefig(re.sub('.csv','transcription.png',data))

df_J=df[df['Nog_Category']=='J']
g= sns.clustermap(df_J.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_J['Function_GVOG'])
g.ax_heatmap.set_title('Translation')
plt.savefig(re.sub('.csv','translation.png',data))

df_M=df[df['Nog_Category']=='M']
g= sns.clustermap(df_M.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_M['Function_GVOG'])
g.ax_heatmap.set_title('Cell membrane biogenesis')
plt.savefig(re.sub('.csv','cell_membrane.png',data))

df_H=df[df['Nog_Category']=='H']
g= sns.clustermap(df_H.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_H['Function_GVOG'])
g.ax_heatmap.set_title('Coenzyme metabolism')
plt.savefig(re.sub('.csv','coenzyme.png',data))

df_Q=df[df['Nog_Category']=='Q']
g= sns.clustermap(df_Q.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_Q['Function_GVOG'])
g.ax_heatmap.set_title('Secondary Metabolites')
plt.savefig(re.sub('.csv','secpmdary_metabolites.png',data))

df_A=df[df['Nog_Category']=='A']
df_A['Function_GVOG']=df_A['Function_GVOG'].str.split('|').str[0]
g= sns.clustermap(df_A.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_A['Function_GVOG'])
g.ax_heatmap.set_title('Rna Processing')
plt.savefig(re.sub('.csv','rna_processing.png',data))

df_O=df[df['Nog_Category']=='O']
g= sns.clustermap(df_O.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_O['Function_GVOG'])
g.ax_heatmap.set_title('Post Translational Modification')
plt.savefig(re.sub('.csv','post_trans_modification.png',data))

df_E=df[df['Nog_Category']=='E']
df_E['Function_GVOG']=df_E['Function_GVOG'].str.split('|').str[0]
g= sns.clustermap(df_E.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_E['Function_GVOG'])
g.ax_heatmap.set_title('Amino Acid Metabolism')
plt.savefig(re.sub('.csv','aa_metabolism.png',data))

df_C=df[df['Nog_Category']=='C']
df_C['Function_GVOG']=df_C['Function_GVOG'].str.split('|').str[0]
g= sns.clustermap(df_C.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_C['Function_GVOG'])
g.ax_heatmap.set_title('Energy Production')
plt.savefig(re.sub('.csv','energy_production.png',data))

df_B=df[df['Nog_Category']=='B']
df_B['Function_GVOG']=df_B['Function_GVOG'].str.split('|').str[0]
g= sns.clustermap(df_B.iloc[:,1:41], col_cluster=False, cmap='PuBu', figsize=(10, 10), yticklabels=df_B['Function_GVOG'])
g.ax_heatmap.set_title('Chromatin Structure')
plt.savefig(re.sub('.csv','chromatin.png',data))

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

print(sums)