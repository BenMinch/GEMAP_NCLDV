#genome annotations
import pandas as pd
import sys,re,os,subprocess
annotation= pd.read_csv(sys.argv[2], sep='\t')
listdir=sys.argv[1]
for i in os.listdir(listdir):
    path=os.path.join(listdir,i)
    data=pd.read_csv(path,sep=' ')
#only keep first 2 columns
    data= data.iloc[:,0:2]
    #remove rows with NA in the second column
    data=data.dropna(subset=['query'])
    data['function']= ''
    data['nog_desc']= ''
    data['nog_categories']= ''
    data['pfam_descs']= ''
    data['function']= data['Hit'].map(annotation.set_index('GVOG')['NCVOG_descs'])
    data['nog_desc']= data['Hit'].map(annotation.set_index('GVOG')['nog_descs'])
    data['nog_categories']= data['Hit'].map(annotation.set_index('GVOG')['nog_categories'])
    data['pfam_descs']= data['Hit'].map(annotation.set_index('GVOG')['pfam_descs'])
#writte the data to a file
    data.to_csv(i+'.tsv', sep='\t', index=False)
