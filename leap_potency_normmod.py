"""
    This function is a motherfunction that executes all parallel normative
    modelling routines. Different specifications are possible using the sub-
    functions.
    ** Input:
        * processing_dir     -> Full path to the processing directory. 
                                This will be the directory which you have   
                                cloned from Rindkind/test_normative_modeling.
        * python_path        -> Full path to the python distribution
        * normative_path     -> Full path to the normative.py
        * job_name           -> Name for the bash script that is the output of
                                this function
        * covfile_path       -> Full path to a .txt file that contains all
                                covariats (subjects x covariates) for the
                                responsefile
        * respfile_path      -> Full path to a .txt that contains all features
                                (subjects x features)
        * batch_size         -> Number of features in each batch
        * memory             -> Memory requirements written as string
                                for example 4gb or 500mb
        * duation            -> The approximate duration of the job, a string
                                with HH:MM:SS for example 01:01:01
        * cv_folds           -> Number of cross validations
        * testcovfile_path   -> Full path to a .txt file that contains all
                                covariats (subjects x covariates) for the
                                testresponse file
        * testrespfile_path  -> Full path to a .txt file that contains all
                                test features
        * log_path           -> Pathfor saving log files
        * binary             -> If True uses binary format for response file
                                otherwise it is text
    written by (primarily) T Wolfers, (adapted) SM Kia
"""


# run normative parallel - set up the scripts to call the LEAP potency files-
# setup




#%%

# Basic setup

import nispat
import os

base_dir=           '/project/3015045.07/norm_mod/leap_potency_abs_harm/'
python_path=        '/home/mrstats/triloo/software/anaconda3/envs/norm_mod/bin/python'
normative_path=     '/home/mrstats/triloo/software/anaconda3/envs/norm_mod/lib/python3.7/site-packages/nispat-0.12-py3.7.egg/nispat/normative.py'
cv_folds=           10
viscovfile_path=    '/project/3015045.07/norm_mod/leap_potency_harm/visualcovar.txt'
batch_size=         1403 
memory=             '4gb'
duration=           '01:50:00'

#%%

# run normative modeling script for prediction

for task in ['hariri','flanker', 'reward_s','reward_m', 'tom']:

    home_dir=           base_dir+task+'/'
    
    if not os.path.exists(home_dir+'prediction/'):
        os.mkdir(home_dir+'prediction/')
    
    processing_dir=     home_dir+'prediction/'
    job_name=           task+'_pr'
    covfile_path=       home_dir+'covar_TD_'+task+'.txt'
    respfile_path=      home_dir+'feats_TD_'+task+'.txt'
    testcovfile_path=   home_dir+'covar_ASD_'+task+'.txt'
    testrespfile_path=  home_dir+'feats_ASD_'+task+'.txt'
    
    print(testrespfile_path)

    nispat.normative_parallel.execute_nm(processing_dir, 
                                         python_path, 
                                         normative_path, 
                                         job_name, 
                                         covfile_path, 
                                         respfile_path, 
                                         batch_size, 
                                         memory, 
                                         duration,
                                         testcovfile_path = testcovfile_path,
                                         testrespfile_path = testrespfile_path
                                         )


#%%
# run normative modeling script for cross validation


for task in ['hariri','flanker', 'reward_s','reward_m', 'tom']:

    home_dir=           base_dir+task+'/'
    
    if not os.path.exists(home_dir+'crossval/'):
        os.mkdir(home_dir+'crossval/')
        
        
    processing_dir=     home_dir+'crossval/'
    job_name=           task+'_cv'
    covfile_path=       home_dir+'covar_TD_'+task+'.txt'
    respfile_path=      home_dir+'feats_TD_'+task+'.txt'    
    testcovfile_path= None
    testrespfile_path= None

    nispat.normative_parallel.execute_nm(processing_dir, 
                                         python_path, 
                                         normative_path, 
                                         job_name, 
                                         covfile_path, 
                                         respfile_path, 
                                         batch_size, 
                                         memory, 
                                         duration,
                                         testcovfile_path = testcovfile_path,
                                         testrespfile_path = testrespfile_path,
                                         cv_folds=cv_folds
                                         )


#%%
# run normative modeling script for forward model

cv_folds = None

for task in ['hariri','flanker', 'reward_s','reward_m', 'tom']:

    home_dir=           base_dir+task+'/'
    
    if not os.path.exists(home_dir+'forward/'):
        os.mkdir(home_dir+'forward/')
        
    processing_dir=     home_dir+'forward/'
    job_name=           task+'_fo'
    covfile_path=       home_dir+'covar_TD_'+task+'.txt'
    respfile_path=      home_dir+'feats_TD_'+task+'.txt'
    testcovfile_path= viscovfile_path
    testrespfile_path = None
    
    nispat.normative_parallel.execute_nm(processing_dir, 
                                         python_path, 
                                         normative_path, 
                                         job_name, 
                                         covfile_path, 
                                         respfile_path, 
                                         batch_size, 
                                         memory, 
                                         duration,
                                         cv_folds = cv_folds,
                                         testcovfile_path = testcovfile_path,
                                         testrespfile_path = testrespfile_path
                                         )

#%%
# combine batches for all tasks and kinds

for task in ['hariri','flanker', 'reward_s','reward_m', 'tom']:
    for kind in ['prediction','crossval','forward']:
        
        print(task, kind)

        processing_dir='/project/3015045.07/norm_mod/leap_potency_abs_harm/'+task+'/'+kind+'/'
        nispat.normative_parallel.collect_nm(processing_dir, collect=True)


#%%
        
        
        
#testing
test=nispat.fileio.load('/project/3015045.07/norm_mod/leap_potency_harm/flanker/covars_ASD_hariri.txt')        
        
        
        
