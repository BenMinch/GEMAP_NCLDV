#Renaming files in a directory
import pandas as pd
import os,sys,subprocess,re, argparse
import numpy as np
import numbers


parser=argparse.ArgumentParser()
parser.add_argument('-c', '--coverm', help='path to coverm directory')
parser.add_argument('-r', '--renamed', help='path to renamed directory')
parser.add_argument('-ref', '--reference', help='path to reference file')
parser.add_argument('-com','--combined', help='path to combined directory')
parser.add_argument('-o','--output', help='path to output directory')
args=parser.parse_args()
coverm=args.coverm
renamed=args.renamed
reference=args.reference
combined=args.combined
output=args.output
reference=pd.read_csv(reference)

#rename files
for i in os.listdir(coverm):
    name=i.split('_')[0]
    for x in range(len(reference)):
        if name==reference['SRR_run'][x]:
            newname=reference['Date'][x]+ '_'+reference['Sample_Name'][x]+'.coverm'
            path=os.path.join(coverm,i)
            newpath=os.path.join(renamed,newname)
            move= 'cp '+path+' '+newpath
            subprocess.call(move, shell=True)

#concatinate all second column of all files in a folder
files = os.listdir(renamed)
for i in files:
    name=i.split('.')[0]
    path=os.path.join(renamed,i)
    df=pd.read_csv(path,sep='\t')
    #change the 2nd column name to the file name
    df.columns.values[1]=name
    df=df.iloc[:,1]
    df.to_csv(os.path.join(combined,i),sep='\t',index=False)

#paste all the files together
paste='paste '+combined+'/* > '+combined+'/all.combined'
subprocess.call(paste, shell=True)

#averaging the days
inputs=combined+'/all.combined'
data= pd.read_csv(inputs, sep='\t')
#split column names at _ and take the first element
data.columns= data.columns.str.split('_').str[0]
#if two columns have the same name, average them

data2= data.groupby(by=data.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])


outname= output
data2.to_csv(outname,index=False)
