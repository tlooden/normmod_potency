# -*- coding: utf-8 -*-
"""
Analyzing the norm_mod output
"""

#import statsmodels.stats.multitest
import scipy.stats
import nispat
import numpy as np
import seaborn as sns
import statsmodels.stats.multitest
import pandas as pd
import matplotlib.pyplot as plt

tasks=['flanker','hariri','reward_m','reward_s','tom']
#%%

"""
Function that takes as input Z.txt files (features x observations) from a 
normative modeling analysis.
Returns two binarized matrices. One containing the significant positive 
deviations, and one containing the significant negative deviations - in 
numpy format.


Parameters -----------------

inputfile:  Enter the location of your Z.txt file

threshold:  you can use a predetermined numerical Z-threshold, otherwise  
            it will default to fdr-correction.
            
persub:     Enter 'True' for this parameter for the function to return one 
            value per participant wrt their # deviant edges.

"""



def threshmat(inputfile, threshold='fdr', variant='comb'):
    
    inputfile=np.array(nispat.fileio.load(inputfile))
    
    if threshold=='fdr':
        #preload matrices for new entries
        posthreshmat=np.zeros(np.shape(inputfile))
        negthreshmat=np.zeros(np.shape(inputfile))
        
        #FDR correction is done on a per feature basis
        for index,row in enumerate(inputfile):
            
            #Transform to p-values for positive and negative deviations
            pospvalrow= scipy.stats.norm.sf(row)
            negpvalrow= scipy.stats.norm.sf(-row)
            
            #FDR-correcting the p-values and returning binarized matrix
            posthreshmat[index] = statsmodels.stats.multitest.fdrcorrection(pospvalrow)[0]
            negthreshmat[index] = statsmodels.stats.multitest.fdrcorrection(negpvalrow)[0]
            
    else:
        #Check simply if Z-value exceeds predetermined threshold
        posthreshmat=np.where(inputfile>threshold, 1, 0)
        negthreshmat=np.where(inputfile<(threshold*-1), 1, 0)
    
    if variant=='pos':
        return posthreshmat
    elif variant=='neg':
        return negthreshmat
    elif variant=='comb':
        # Combine both matrices and sum over edges to get one value per subject
        combthreshmat=posthreshmat+negthreshmat
        return combthreshmat
        

#%%

# For TD vs ASD

threshold=2
variant='comb'


# Generate long-form arrays for DataFrame

for i,task in enumerate(tasks):
    inputfile_control='/project/3015045.07/norm_mod/leap_potency_abs_harm/'+task+'/crossval/Z.txt'
    inputfile_case='/project/3015045.07/norm_mod/leap_potency_abs_harm/'+task+'/prediction/Z.txt'
    
    control = threshmat(inputfile_control, 2, variant=variant)
    case = threshmat(inputfile_case, 2, variant=variant)
    
    sumcontrol=np.mean(control,0)*100
    sumcase=np.mean(case,0)*100
    
    
    a = np.append(np.zeros(len(sumcontrol)), np.ones(len(sumcase)))
    b = np.array([i for x in range(len(a))])
    c = np.append(sumcontrol,sumcase)
    d = np.stack([a,b,c])
    
    
    if i==0:
        e=d
    else:
        e=np.hstack([e,d])


#%%

# For TD vs ASDadhd- vs ASDadhd+

threshold=2
variant='pos'


# Generate long-form arrays for DataFrame

for i,task in enumerate(tasks):
    inputfile_control='/project/3015045.07/norm_mod/leap_potency_abs_harm/'+task+'/crossval/Z.txt'
    inputfile_case='/project/3015045.07/norm_mod/leap_potency_abs_harm/'+task+'/prediction/Z.txt'
    
    control = threshmat(inputfile_control, threshold, variant=variant)
    case = threshmat(inputfile_case, threshold, variant=variant)
    
    # Sum over edges and multiply by 100 to get percentage
    sumcontrol=np.mean(control,0)*100
    sumcase=np.mean(case,0)*100
    
    # Generate dummy for diagnostic groups
    dumcontrol  = np.zeros(len(sumcontrol))
    dumasd      = np.ones(len(sumcase))
    dumasd[list1[i]]=2
    
    
    a = np.hstack([dumcontrol, dumasd])
    b = np.array([i for x in range(len(a))])
    c = np.append(sumcontrol,sumcase)
    d = np.stack([a,b,c])
    
    
    if i==0:
        e=d
    else:
        e=np.hstack([e,d])
        
        
#%%

# Build into datframe
df_neg=pd.DataFrame(e.T, columns=[ 'diagnosis', 'task','% deviant edges'])

# Rename variables inside the dataframe to meaningful terms.
df_neg['task']=df_neg['task'].replace({0: 'flanker', 1: 'hariri', 2: 'reward_m', 3: 'reward_s', 4 : 'tom'})
df_neg['diagnosis']=df_neg['diagnosis'].replace({0: 'TD', 1: 'ASD', 2:'ASDadhd'})

#%%

result=scipy.stats.ttest_ind(d,d1)

mean1=np.mean(d)
mean2=np.mean(d1)
std1=np.std(d)
std2=np.std(d1)

#%%

# Plotting
sns.set(style="whitegrid")
sns.violinplot(x='task',y='% deviant edges', hue='diagnosis',split=False, inner="quartile", data=df_pos)
plt.legend(loc='lower right')

plt.savefig("/home/mrstats/triloo/Workspace/Images/norm_mod//normmod_abs_pos.pdf")


#%%


matsub_all=np.load('/home/mrstats/triloo/Workspace/Analysis/Aggregate_potency_matrices/matsub_all.npy', allow_pickle=True)

matsub_ADHD=np.loadtxt('/home/mrstats/triloo/Workspace/LEAP_info/ASDadhd+.txt')



#%%

#print(matsub_onlyADHD[0][4])

print(matsub_all[0][0])

#print(matsub_all[0][4][0][0])

#%%

# Generate a list of indexes for the subjects with high ADHD scores corresponding
# to each of the tasks.

list1=[[]for i in range(5)]


for t in range(5):
    for n,i in enumerate(matsub_ADHD):
        i=str(int(i))
        
        try:
            list1[t].append(matsub_all[0][t].index(i))
        except:
            pass
    

#%%


print(len(list1[1]))

np.ones(len(list1[i]))


a=np.append([np.ones(3), np.zeros(3)], np.ones(3))

np.full(len(list1[0]),2)








