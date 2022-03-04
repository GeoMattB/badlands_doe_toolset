import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plot
from scipy.ndimage import gaussian_filter1d
import scipy.interpolate as interp

import similaritymeasures
from scipy.spatial.distance import directed_hausdorff
import h5py

import badlands_doe_toolset.postproc_utils as ppu
import os
from multiprocessing import Pool
import time
from functools import partial




def exp_thickness_doeval(model_extract,wells_file,expattrib='cumdiff',outdir='results_final'):
    #load the extracted data files
    model_extract=h5py.File(model_extract,'r')
    explist=list(model_extract.keys())
    #load the file with the well truth information 
    wells=ppu.Welldata()
    wells.loadWell(wells_file)
        
    expname_in=[]
    wellname_in=[]
    trth_in=[]
    expres_in=[]   
   
    for i in explist:
        for j in range(0,len(wells.well_name)):
            expname_in.append(i)
            wellname_in.append(wells.well_name[j])
            trth_in.append(wells.well_tru[j])
            expres_in.append(model_extract[i+'/'+str(wells.well_name[j])+'/'+expattrib][0][0]) #his is ugly
    
  # make a dataframe for the outputs with columns that DoEval expects    

    results=pd.DataFrame()
    results['Experiment']=expname_in
    results['Well']=wellname_in
    results['Y Truth']=trth_in
    results['Y Exp']=expres_in
                
    # We need an Nexp column and a Y Label column for DoEval to work, strip the non-numerics from the experiment name, add a Y label column with Thickness as string.
    Nexp=[]
    for i in results['Experiment']:
        ns = re.sub("[^0-9]", "", i)
        Nexp.append(ns)
    results['Nexp']=Nexp
    results['Y Label']=['Thickness']*(len(results))
    
    # save this to a csv for evaluation.
    results.to_csv(outdir+'/doeval_well_thickness.csv',index=False) 
    return(results)   

# here we generate outputs for evaluation of the well data in DoEval. A well log/curve similarity to compare stratigraphy and a well thickness ratio 
# are used to classify how well the experiment matches a measured result.
# wellfile is the file with the well locations as specified for postproc_utils (note this file also contains links to the well logs to be compared)

def exp_attr_sim_DoEgen(model_extract,obslog,expattrib,wells_file,minthkRatio,exp_srm=5,log_smth=10,outdir='results_final'):
    model_file_loc=model_extract # hdf5 arrays aren't picklable as the hdf5 attributes are read on demand, send the file loc as a string to the child functions. 
    if type(model_extract)==str: # this is read like this to keep consistency 
        model_extract=h5py.File(model_extract,'r')
    
    explist=list(model_extract.keys())
    proc = os.cpu_count() #just use all the cpus to calculate the similarities.
    print ('Starting well sim comparisons in parallel: '+str(proc)+str(' processes'))
    print(str(len(explist))+' experiments')
    with Pool(processes = proc) as p:
        start = time.time()
        similarity_all = p.map_async(partial(exp_attr_sim,model_file_loc,wells_file,obslog,expattrib,exp_srm,log_smth),explist) #the multiproc start bit
        p.close()
        p.join()
        print('Completed; total processing time: ' + str(round((time.time() - start)/60, 2)) + ' minutes')
    print('writing outputs')
    
 #### Sort out the dataframe for the outputs in a format DoEval expects
    wells=ppu.Welldata()
    wells.loadWell(wells_file)  
    for i in range(0,(len(explist))):
        for j in range(0,len(wells.well_name)):
            # print(i,j)
            dd=pd.DataFrame()
            c=similarity_all.get()[i][j]
            dd=dd.append(c).transpose()
            dd.columns=['Exp_name','Wellname','sumthck','df','dh']
            if i==0 and j==0: #initialise the first row of the dataframe
                sim_results=dd
            else:
                sim_results=pd.concat([dd,sim_results])
    # the index doesn't get created correctly above with pd.concat, reset
    sim_results.reset_index(inplace = True, drop = True)
    # add the colums with the comparison values from wells to the dataframe
    measThick=[0]*(len(sim_results))
    
    #Loop through results array and add the measured well thickness column
    for j in range(0,len(wells.well_name)):
        search=wells.well_name[j]
        trthk=wells.well_tru[j]
        temp=sim_results.loc[sim_results.isin([search]).any(axis=1)].index.tolist()
        for i in temp:
            measThick[i]=trthk
    # We need a numeric Nexp (that matches the DoEgen csv) column for DoEval to work, strip the non-numerics from the experiment name.
    Nexp=[]
    for i in sim_results['Exp_name']:
        ns = re.sub("[^0-9]", "", i)
        Nexp.append(ns)
    sim_results['Nexp']=Nexp
   
    sim_results['measThick']=measThick
    sim_results['thickRatio']=sim_results['sumthck']/sim_results['measThick']
    sim_results['thickDiff']=sim_results['measThick']-sim_results['sumthck']
    sim_results['tru_val']=[1]*(len(sim_results)) # perfect df similarity is 0, this is close enough for our comparison om DoEval.
    #sim_results.fillna(999, inplace=True)
    sim_results.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # now create a similarity column (sim_results_clean) with the low thickness ratio wells tagged as outliers / bad matches in this case using 999.
    sim_results_clean=[]
    tempindex=0 # 
    for i in sim_results['thickRatio']:
        if i < minthkRatio or sim_results['df'][tempindex]==np.nan:
            a=999
        else:a=sim_results['df'][tempindex]
        tempindex +=1
        sim_results_clean.append(a)
        
    sim_results.replace([np.inf, -np.inf], np.nan, inplace=True)
    sim_results['tru_val_similarity']=sim_results_clean
 
    # configure the "Y Label" and "Y Truth" columns for the thicknesses
    sim_results['Y Label']=['ThicknessRatio']*(len(sim_results))
    sim_results['Y Truth']=sim_results['tru_val']
    sim_results['Y Exp']=sim_results['thickRatio']
    # save results to file for DoEval for the thickness evaluations
    sim_results.to_csv(outdir+'/doeval_well_thickRatio.csv',na_rep='NULL',index=False)
 
    # configure the "Y Label" and "Y Truth" columns for the well similarity
    sim_results['Y Label']=sim_results['Wellname']
    sim_results['Y Truth']=sim_results['tru_val']
    sim_results['Y Exp']=sim_results['tru_val_similarity']
    # save results to file for DoEval for the experiment vs log similarity DoEval evaluations
    sim_results.to_csv(outdir+'/doeval_well_similarity.csv',na_rep='NULL',index=False)
   
    return (sim_results)
    

# Exp is the name of the experiment group in the model extract hdf5 file as as string, ie 'experiment_11'. Here it is a default setting for futire use but also to make integration with multiprocessing.pool simpler, to do with the order of args when combining partial & pool.
    
def exp_attr_sim(model_extract,wells_file,obslog,expattrib,exp_srm=5,log_smth=10,Exp='Experiment'):
    #if model_extract is the hdf5 data file then load that, otherwise it's assumed to be already loaded and we can continue
    print(Exp)
    if type(model_extract)==str:
        model_extract=h5py.File(model_extract,'r')
    
    wells=ppu.Welldata()
    wells.loadWell(wells_file) 
        
    Nwells=range(0,len(wells.well_name))
    WellSim=[]
    for i in Nwells:
        well_log_file=str('data/'+wells.well_log[i])
        well_name=str(wells.well_name[i])
        WellSim.append(well_attr_sim(model_extract,Exp,well_name,obslog,well_log_file,expattrib,exp_srm=5,log_smth=10))
    return (WellSim)

def well_attr_sim(model_extract,Exp,Well,obslog,well_log_file,expattrib,exp_srm=5,log_smth=10):
#if model_extract is the hdf5 data file then load that, otherwise it's assumed to be already loaded and we can continue
    if type(model_extract)==str:
        model_extract=h5py.File(model_extract,'r')

    #Read in log data
    log_df=pd.read_csv(well_log_file)
    # check if obslog is in the well log file
    
    loglist=list(log_df.columns)
    if obslog not in loglist:
        print('the attribute you are comparing is not in the observed / measured log file. Select from: '+str(loglist))
        return('no attribute')
    
    obs_log=log_df[obslog]
    obs_depth=log_df['Depth']
    obs_log_smooth = gaussian_filter1d(obs_log, log_smth)

    Exp_name = Exp 
    Wellname = Well 

    Depth=model_extract['/'+Exp_name+'/'+Wellname+'/layDepth'][:] 
    Depth=Depth.ravel()
    
    Thick=model_extract['/'+Exp_name+'/'+Wellname+'/layThick'][:] 
    Thick=Thick.ravel()

    exp_log=model_extract['/'+Exp_name+'/'+Wellname+'/'+expattrib][:] # get the attribute from the model
    exp_log=exp_log.ravel()

    np.seterr(divide='ignore', invalid='ignore') # stops on errors when theres a NaN, just keep going
    
    #cutoff for deepest / shallowest model_extracts
    MaxDepth=(round(np.min(Depth),0))
    MinDepth=int(round(np.max(Depth),0))
    TotThick=MinDepth-MaxDepth
    #MinDepth=0
    #print(MaxDepth, MinDepth)
         
    # Model & measured array values need have the same shape, interpolate both to same using depths.

    f = interp.interp1d(Depth,exp_log, kind='linear', fill_value='extrapolate')
    expDepth=np.arange(MaxDepth,MinDepth,exp_srm)
    expLog=f(expDepth)
                
    g = interp.interp1d(obs_depth,obs_log_smooth,  kind='linear', fill_value='extrapolate')
    obsDepth_interp=np.arange(MaxDepth,MinDepth,exp_srm) #these is the depth the logs get resampled/interpolated to
    obsLog_interpol=g(obsDepth_interp)
        
    x1=expDepth.tolist()
    x2=obsDepth_interp.tolist()
    y1=expLog.tolist()
    y2=obsLog_interpol.tolist()
    P = np.array([x1, y1]).T
    Q = np.array([x2, y2]).T


    #print('calculate sims')
    sumthck=np.sum(Thick)
    try:
        df = similaritymeasures.frechet_dist(P, Q)
        dh, ind1, ind2 = directed_hausdorff(P, Q)
    except (IndexError, RuntimeError, TypeError):
        df='err'
        dh='err'
    #dtw, d = similaritymeasures.dtw(P, Q)
    #pcm = similaritymeasures.pcm(P, Q)
    #area = similaritymeasures.area_between_two_curves(P, Q)
    #cl = similaritymeasures.curve_length_measure(P, Q)
    
    sim = [Exp_name,Wellname,sumthck,df,dh]

    return(sim)
    
def well_attr_plot(model_extract,Exp,Well,obslog,well_log_file,expattrib,exp_srm=5,log_smth=10):    
    #if model_extract is a file then load that, otherwise it's assumed to be already loaded and we can continue
    # check for the output directory
    path='results_final/plots'
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    if type(model_extract)==str:
        model_extract=h5py.File(model_extract,'r')    
    #Read in log data
    log_df=pd.read_csv(well_log_file)
    # check if obslog is in the well log file
    
    loglist=list(log_df.columns)
    if obslog in loglist:
        obs_log=log_df[obslog]
        obs_depth=log_df['Depth']
        obs_log_smooth = gaussian_filter1d(obs_log, log_smth)
    else:
        print('no matching attribute found measured log file. Select from: '+str(loglist))
        #return('no attribute')
        # just create some dummy data in case you ust want to display the experimental plots.
        obs_log=[0,0]
        obs_depth=[0,0]
        obs_log_smooth = [0,0]


    Depth=model_extract['/'+Exp+'/'+Well+'/layDepth'][:]   #/experiment_10/Cast01/Depth')[:]
    Depth=Depth.ravel()
      
    DThick=model_extract['/'+Exp+'/'+Well+'/layinstdiff'][:] 
    DThick=DThick.ravel()
    
    
    if expattrib in model_extract['/'+Exp+'/'+Well+'/'].keys():
        exp_log=model_extract['/'+Exp+'/'+Well+'/'+expattrib][:]
        exp_log=exp_log.ravel()
    else:
        print('no matching attribute in experiment file select from: ',(list(model_extract['/'+Exp+'/'+Well+'/'].keys())))
        exp_log=Depth.ravel()*0
    
    tStart=int(model_extract[Exp].attrs['tStart'])
    tEnd=int(model_extract[Exp].attrs['tEnd'])
    tDisplay=int(model_extract[Exp].attrs['tDisplay'])
    TimeSteps=list(range(tStart,tEnd+tDisplay,tDisplay)) 
    #print(TimeSteps)
    
    np.seterr(divide='ignore', invalid='ignore') # stops on errors when theres a NaN, just keep going
    
    #cutoff for deepest / shallowest model_extracts
    MaxDepth=(round(np.min(Depth),0))
    MinDepth=int(round(np.max(Depth),0))
         
    # Model & measured array values need to have the same shape, interpolate both to same using depths.
    
    f = interp.interp1d(Depth,exp_log, kind='linear', fill_value='extrapolate')
    expDepth=np.arange(MaxDepth,MinDepth,exp_srm)
    expLog_interpol=f(expDepth)
    
    dt = interp.interp1d(Depth,TimeSteps, kind='linear', fill_value='extrapolate')
    expTimesteps=dt(expDepth)
         
    g = interp.interp1d(obs_depth,obs_log_smooth,  kind='linear', fill_value='extrapolate')
    obsDepth_interp=np.arange(MaxDepth,MinDepth,exp_srm) #these is the depth the logs get resampled/interpolated to yes, its a duplicate of the previous stuff
    obsLog_interpol=g(obsDepth_interp)
        
    fig, (ax1, ax2,) = plot.subplots(1, 2, figsize=(3,6))
    
    title=(str(Exp) + ' : \n ' + str(Well))
    fig.suptitle(title, fontsize=16, x=0.05)
    plot.subplots_adjust(right=3)
    
    ax1.set_title('Depth at well location')
    ax1.set(ylabel="Well Depth")
    ax1.set(xlabel=obslog+' / '+expattrib)
    #ax1.set_xlim([-50, 250])
    ax1.plot(expLog_interpol,expDepth,'b-', label=str(expattrib))
    ax1.plot(obsLog_interpol,obsDepth_interp,'g--',label=str(obslog))
    ax1.legend(loc="center right")
    ax1a=ax1.twiny()
    #ax1a.set_xlim([-25,125])
    ax1a.set_xticks(range(0,125,25))
    ax1a.set_xticklabels(range(0,125,25))
    ax1a.set(xlabel='Rate of dep m/100kyr')
    ax1a.plot(DThick,Depth,'y--',alpha=0.7, label='Rate of Dep.')
    ax1a.legend(loc="upper right")
        
                    
    ax2.set_title('Time / model step')
    ax2.set(ylabel='Time years')
    ax2.set_ylim([tStart, tEnd])
    ax2.set(xlabel=obslog+' / '+expattrib)
    #ax2.set_xlim([-50, 250])
    
    ax2.plot(expLog_interpol,expTimesteps,'b-', label=str(expattrib))
    ax2.plot(obsLog_interpol,expTimesteps,'g--',label=str(obslog))
    ax2.legend(loc="center right")
        
    textstr = 'avg: '+str(round(np.mean(DThick),1))+' m/100kyr' 
    ax2.text(0.6, 0.13, textstr, transform=ax2.transAxes, fontsize=12, horizontalalignment='left')
             
    ax2a=ax2.twiny()
    #ax2a.set_xlim([-25, 125])
    ax2a.set_xticks(range(0,125,25))
    ax2a.set_xticklabels(range(0,125,25))
    ax2a.set(xlabel='Rate of dep m/100kyr')
    ax2a.plot(DThick,TimeSteps,'y--',alpha=0.7, label='Rate of Dep.') 
    ax2a.legend(loc="upper right")
    
#    if writefig==1:   
    plot.savefig('results_final/plots/'+str(Exp)+'_'+str(Well)+ '.png', dpi=125,bbox_inches='tight')
    plot.show()


    

def exp_attr_plot(model_extract,Exp,obslog,well_loc,explog,exp_srm=5,log_smth=10):
    if type(model_extract)==str: # this is read like this to keep consistency
        model_extract=h5py.File(model_extract,'r')
    wells=ppu.Welldata()
    wells.loadWell(well_loc)  
    for i in range(0,len(wells.well_name)):
        well_log_file='data/'+wells.well_log[i]
        Well = wells.well_name[i]
        well_attr_plot(model_extract,Exp,Well,obslog,well_log_file,explog,exp_srm,log_smth)
     
def well_attr_plot_DoEgen(model_extract,obslog,well_loc,explog,exp_srm=5,log_smth=10):
#    model_file_loc=model_extract # hdf5 arrays aren't picklable as the hdf5 attributes are read on demand, send the file loc as a string to the child functions. 
    #if type(model_extract)==str: # this is read like this to keep consistency 
    model_extract=h5py.File(model_extract,'r')
    explist=list(model_extract.keys())
    wells=ppu.Welldata()
    wells.loadWell(well_loc)  
    #wellList=wells.well_name
    for Exp in explist:
        for i in range(0,len(wells.well_name)):
            well_log_file='data/'+wells.well_log[i]
            Well = wells.well_name[i]
            well_attr_plot(model_extract,Exp,Well,obslog,well_log_file,explog,exp_srm,log_smth)
