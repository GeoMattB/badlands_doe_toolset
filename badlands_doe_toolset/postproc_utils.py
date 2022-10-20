import h5py
import numpy as np
from scipy import spatial
import pandas as pd
import networkx as nx
import scipy.interpolate as interp
from badlands.model import xmlParser
import badlands_companion.stratalMesh as mesh
from multiprocessing import Pool
import multiprocessing
import threading
import os
import sys
import time
import badlands_doe_toolset.badlands_multiproc_run as mpr
import ntpath
from functools import partial
import matplotlib.pyplot as plt
from pyevtk.hl import gridToVTK
from pathlib import Path

""" this Stratdata class might be removed and the Stratigraphy class from StratalMesh used instead, however I don't need to reshape the arrays here"""
class Stratadata:
    """
    Class for loading data from stratal files.
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.dx = None
        self.nx = None
        self.ny = None
        
    
    def loadStrat(self,stratfile):
        """
        Read the Badlands HDF5 stratal file (sed.timeXXX.hdf5).
        """
        df = h5py.File(stratfile, 'r+')
        coords = np.array((df['/coords']))
        x,y = np.hsplit(coords, 2)
        self.x=x
        self.y=y
        self.dx = x[1]-x[0]
        self.nx = int((x.max() - x.min())/self.dx+1)
        self.ny = int((y.max() - y.min())/self.dx+1)
        self.nz=(df['layDepth']).shape[1]
        #Load all of the attributes with 'lay' in their name, these are the layer properties.
        self.strat_attribs=[]
        for i in df.keys():
            if 'lay' in i:
                self.strat_attribs.append(i)
        self.strat_layers={}
        for j in self.strat_attribs:
            self.strat_layers[j] = df['/'+j]
           
    # in some cases only the coordinates are required and this is faster to load (modified from badlands companion)   
    def loadStratCoords(self,stratfile): 
        """
        Read the Badlands HDF5 stratal file (sed.timeXXX.hdf5) minimal information.
        """
        df = h5py.File(stratfile, 'r+')
        coords = np.array((df['/coords']))
        x,y = np.hsplit(coords, 2)
        self.x=x
        self.y=y
        self.dx = x[1]-x[0]
        self.nx = int((x.max() - x.min())/self.dx+1)
        self.ny = int((y.max() - y.min())/self.dx+1)
        self.nz=(df['layDepth']).shape[1]
       
        
class TINfile:
    """
    Class for loading data from tin.timeXX files.
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.area=None
        self.cells = None
        self.cumdiff = None
        self.cumfail = None
        self.discharge = None
        self.lake = None
        
    
    def loadTIN(self,tinfile):
        """
        Read the Badlands HDF5 tin.timeXX file.
        """
        tin_df = h5py.File(tinfile, 'r+')
        coords = np.array((tin_df['/coords']))
        self.x=coords[:,0]
        self.y=coords[:,1]
        self.z=coords[:,2]
        self.area=np.array((tin_df['/area']))
        self.cells = np.array((tin_df['/cells']))-1

        self.tin_keylist=[]
        for i in tin_df.keys():
            if i not in ['area','cells','coords']:
                self.tin_keylist.append(i)
        self.tin_layer={}
        for j in self.tin_keylist:
            self.tin_layer[j] = tin_df['/'+j]
        
class Flowfile:
    """
    Class for loading data from flow.timeXX files.
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.basin=None
        self.chi = None
        self.connect = None
        self.discharge = None
        self.flowdensity = None
        self.sedload = None
        self.fetch = None
    
    def loadFlow(self,flowfile):
        """
        Read the Badlands HDF5 flow.timeXX file.
        """
        flw_df = h5py.File(flowfile, 'r+')
        self.coords = np.array((flw_df['/coords']))
        self.x=self.coords[:,0]
        self.y=self.coords[:,1]
        self.z=self.coords[:,2]
        self.connect = np.array((flw_df['/connect']))-1
        self.basin=np.array((flw_df['/basin']))
       
        self.flw_keylist=[]
        for i in flw_df.keys():
            if i not in ['connect','coords']:
                self.flw_keylist.append(i)
        self.flw_layer={}
        for j in self.flw_keylist:
            self.flw_layer[j] = flw_df['/'+j]

"""
Functions to generate post-processing outputs go here.
interp_from_tin interpolates a single tin.timeXX.hdf5 onto a higher resolution mesh.
strat_tin_write reads a the final stratal file mesh coordinates and calls interp_from_tin for each timestep, then writes the new attribute to the startal file 
strat_tin_write_doegen reads a list of experiment xml files from a common directory and calls strat_tin_write for each experiment.
"""

### interpolates data from 1 tin file onto the XY coords of a stratal sed.time grid
def interp_from_tin(X,Y,nx,ny,tin_file,attribs):
    #use the tin class to load the tin array
    tin=TINfile()
    tin.loadTIN(tin_file)
    attribs_inter=[]
    #interpolate the data onto the strata array
    for i in attribs:
        #print(attribs[i])
        #attrib_temp = interp.griddata((tin.x, tin.y), getattr(tin,attribs[i]) , (X, Y), method='nearest') #
        attrib_temp = interp.griddata((tin.x, tin.y), tin.tin_layer[i][:] , (X, Y), method='nearest')
        attrib_temp = attrib_temp.flatten()
        attribs_inter.append(attrib_temp)
    return(attribs_inter)


#This function writes the instdiff parameter to the tin.time files. Useful for time-step visualisation in paraview, 
# could probably be done entirely in PV as a filter but I couldn't work out how to do it dynamically for animations in paraview.
# to do use multiproc to speed this up
def tin_write_instdiff(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)

    #loop through the tin.timeXX loading each as we go.
    for i in range(0,maxSteps):
        tin_N=TINfile()
        tin_N1=TINfile()
        tin_N.loadTIN(modelh5dir+'/tin.time'+str(i)+'.hdf5')
        #return(tin_N)
        tin_N1.loadTIN(modelh5dir+'/tin.time'+str(i+1)+'.hdf5')
        
        #this is a bit not-that-useful, it just adds an attribute to the first (zero) layer to keep outputs consistent / tidy.
        if i==0:
            instdiff=tin_N.tin_layer['cumdiff'][:]
            with h5py.File((modelh5dir+'/tin.time'+str(i)+'.hdf5'), 'a') as t:
                print('tin.time'+str(i))
                if ('instdiff') in t.keys():
                    del t['instdiff'] # if the value exists already remove it
                t.require_dataset('instdiff', data=instdiff, shape=instdiff.shape, dtype="f",compression="gzip" )
                t.close()
        
        # subtracts (layer N) from (layer N+1) and adds it to the tin file, this is the main bit.
        else:
            instdiff=(tin_N1.tin_layer['cumdiff'][:])-(tin_N.tin_layer['cumdiff'][:])
        with h5py.File((modelh5dir+'/tin.time'+str(i+1)+'.hdf5'), 'a') as t:
            print('tin.time'+str(i+1))
            if ('instdiff') in t.keys():
                del t['instdiff'] # if the value exists already remove it
            t.require_dataset('instdiff', data=instdiff, shape=instdiff.shape, dtype="f",compression="gzip" )
            t.close()


#This function writes the list of supplied parameters in a models tin.time files to the final sed.time file.
def strat_tin_write(attribs,modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    strat_file=modelh5dir+'/sed.time'+str(maxSteps)+'.hdf5'
    
    #load in the XY coords from the stratal file.
    mesh = Stratadata()
    mesh.loadStratCoords(strat_file) #get the sed.time file coordinates only to interpolate onto (faster).
    X=mesh.x
    Y=mesh.y
    ny=mesh.ny
    nx=mesh.nx
    # create the dataframes for the outputs from the list of required attributes
    dataout = {}
    for i in attribs:
        dataout[i] = pd.DataFrame()
    print('adding '+str(attribs)+' to '+str(strat_file))
    #loop through the tin.timeXX files to interpolate the attributes onto the strata.
    for i in range(0,maxSteps+1):
        attributes_interpolate_results=interp_from_tin(X,Y,nx,ny,modelh5dir+'/tin.time'+str(i)+'.hdf5',attribs)
        # loop over the results to populate the dict/dataframes
        for j in range(0,len(attributes_interpolate_results)):       
            dataout[attribs[j]][i]=attributes_interpolate_results[j]
    
    # dataout now needs to be written to the final sed.time file, open the strat file, append each attribute
    # note, the require_dataset function won't overwrite the attribute if it already exists, so the attribute is first removed if it exists using the if / in statement for each.     
    with h5py.File(strat_file, 'a') as g:
        print('writing to file')
        for i in attribs:
            if i =='cumdiff': 
                if ('layinstdiff') in g.keys():
                    del g['layinstdiff'] # if the value exists already remove it
                print('cumdiff specified, generating instdiff as well')
                dataout['instdiff']=np.copy(dataout['cumdiff'])
                dataout['instdiff'][:, 1:] -= dataout['instdiff'][:, :-1]  #cumdiff is more useful as a rate of depostion for that step, subtract previous column (timestep) to get deposition/erosion per step.
                dataout['instdiff'][:, 0] = 0
                g.require_dataset(str('layinstdiff'), data=dataout['instdiff'], shape=dataout['instdiff'].shape, dtype="f",compression="gzip" )

            if i =='cumhill': 
                if ('layinsthill') in g.keys():
                    del g['layinsthill'] # if the value exists already remove it
                print('cumhill specified, generating insthill as well')
                dataout['insthill']=np.copy(dataout['cumhill'])
                dataout['insthill'][:, 1:] -= dataout['insthill'][:, :-1]  #cumhill is more useful as a rate of depostion for that step, subtract previous column (timestep) to get deposition/erosion per step.
                dataout['insthill'][:, 0] = 0
                g.require_dataset(str('layinsthill'), data=dataout['insthill'], shape=dataout['insthill'].shape, dtype="f",compression="gzip" )

            else:
                if ('lay'+str(i)) in g.keys():
                    del g['lay'+str(i)] # if the value exists already remove it
                g.require_dataset('lay'+str(i), data=dataout[i], shape=dataout[i].shape, dtype="f",compression="gzip" )
        g.close()
        print('completed '+strat_file)


#This function writes the specified attributes from tin files to the final strata output in a set of DoEgen experiments
def strat_tin_write_doegen(xml_dir,attribs,proc=None):
#    import lib.badlands_multiproc_run as mpr
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
    # xml_file_loc=[]
    # for i in xm.xml_list:
    #    a = xml_dir+'/'+i
    #    xml_file_loc.append(a)
    #    print (a)     
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount)) 
    if proc == None:
        if cpuCount >= 2:
            proc = int(cpuCount/2) #Number of Processes to run simultaneously X-2 basically uses all the threads available this has been tested to ~20 threads. Bottlenecks   could be the size of the strat file and hard drive r/w speed.
        else:
            proc = 1 # simple error check.
    print ("number of processes assigned: "+str(proc))
    with Pool(processes = proc) as p:
        start = time.time()
        async_result = p.map_async(partial(strat_tin_write,attribs), xm.xml_list)
        p.close()
        p.join()
        print('All outputs written to final stratigraphy file in depth: ' + str(round((time.time() - start)/60, 2)) + ' minutes')
#################

# As above for tin interp, this interpolates a single file or timesteps flow attributes onto a stratal grid
def interp_from_flow(X,Y,nx,ny,flow_file,attribs):
    #use the Flow class to load the flow array
    flow=Flowfile()
    flow.loadFlow(flow_file)
    attribs_inter=[]
    #interpolate the data onto the strata array
    for i in attribs:
        #attrib_temp = interp.griddata((flow.x, flow.y), getattr(flow,attribs[i]) , (X, Y), method='nearest') #
        attrib_temp = interp.griddata((flow.x, flow.y), flow.flw_layer[i][:] , (X, Y), method='nearest')
        attrib_temp = attrib_temp.flatten()
        attribs_inter.append(attrib_temp)
    return(attribs_inter)

#This writes the interp onto the final stratal file for all of flow files in a model.
def strat_flow_write(attribs,modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    strat_file=modelh5dir+'/sed.time'+str(maxSteps)+'.hdf5'
    
    #load in the XY coords from the stratal file.
    mesh = Stratadata()
    mesh.loadStratCoords(strat_file) #get the sed.time file coordinates only to interpolate onto (faster).
    X=mesh.x
    Y=mesh.y
    ny=mesh.ny
    nx=mesh.nx
    # create the dataframes for the outputs from the list of required attributes
    dataout = {}
    for i in attribs:
        dataout[i] = pd.DataFrame()
    print('adding '+str(attribs)+' to '+str(strat_file))
    #loop through the tin.timeXX files to interpolate the attributes onto the strata.
    for i in range(0,maxSteps+1):
        attributes_interpolate_results=interp_from_flow(X,Y,nx,ny,modelh5dir+'/flow.time'+str(i)+'.hdf5',attribs)
        # loop over the results to populate the dict/dataframes
        for j in range(0,len(attributes_interpolate_results)):       
            dataout[attribs[j]][i]=attributes_interpolate_results[j]
    
    # dataout now needs to be written to the final sed.time file, open the strat file, append each attribute
    # note, this is an append function and if the attribute already exists it won't be overwritten.    
    with h5py.File(strat_file, 'r+') as g:
        print('writing to file')
        for i in attribs:
            if ('lay'+str(i)) in g.keys():
                del g['lay'+str(i)] # if the value exists already remove it
            g.require_dataset('lay'+str(i), data=dataout[i], shape=dataout[i].shape, dtype="f",compression="gzip" )
        g.close()
        print('completed '+strat_file)


#This function writes the flow files to the final strata output for a set of DoEgen experiments
def strat_flow_write_doegen(xml_dir,attribs):
#    import lib.badlands_multiproc_run as mpr
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
   
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount)) 
    if cpuCount >= 2:
        proc = int(cpuCount-1) #Number of Processes to run simultaneously X-2 basically uses all the threads available this has been tested to ~20 threads. Bottlenecks   could be the size of the strat file and hard drive r/w speed.
    else:
        proc = 1 # simple error check.
    print ("number of processes assigned: "+str(proc))
    with Pool(processes = proc) as p:
        start = time.time()
        async_result = p.map_async(partial(strat_flow_write,attribs), xm.xml_list)
        p.close()
        p.join()
        print('All outputs written to final stratigraphy file in depth: ' + str(round((time.time() - start)/60, 2)) + ' minutes')


# This uses Poro0 and PoroC to decompact the thickness in the final strala (sed.time) file output.
def strat_depothick_write(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    poro0=model.poro0
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    strat_file=modelh5dir+'/sed.time'+str(maxSteps)+'.hdf5'
    strat=Stratadata()
    strat.loadStrat(strat_file)
    p=strat.strat_layers['layPoro']
    th=strat.strat_layers['layThick']
    d=(1+(1-(p/poro0)))*th # 
    with h5py.File(strat_file, 'a') as g:
        g.require_dataset('laydepoThick', data=d, shape=d.shape, dtype="f", compression="gzip"  )
        g.close()

def strat_depothick_write_doegen(xml_dir):
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount)) 
    if cpuCount >= 2:
        proc = int(cpuCount-1) #Number of Processes to run simultaneously X-2 basically uses all the threads available this has been tested to ~20 threads. Bottlenecks   could be the size of the strat file and hard drive r/w speed.
    else:
        proc = 1 # simple error check.
    print ("number of processes assigned: "+str(proc))
    with Pool(processes = proc) as p:
        start = time.time()
        async_result = p.map_async(partial(strat_depothick_write), xm.xml_list)
        p.close()
        p.join()
        print('depothick written to final stratigraphy file in depth: ' + str(round((time.time() - start)/60, 2)) + ' minutes')

#cumdiff is total, cumulative additional deposition. We want a step wise accumulation rate something like cumdiff@200kyrs minus cumdiff@100kyrs. It's useful to have this in both the strat and 
#def instdiff_calc(xml_dir):
    


"""
Terrain_roughness interpolates the elevation onto the regular, stratalMesh grid             
then uses a sliding window to calculate Rileys (Riley 1999) Terrain Roughness Index.
Scale is important for this indicator and whle this uses a 3x3 sliding window it may be that larger say, 5x5 or 9x9 sliding windows are useful
this could be modified to vary the sliding window size later.

"""

def strat_TerrRough_write_doegen(xml_dir):
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount)) 
    if cpuCount >= 2:
        proc = int(cpuCount-1) #Number of Processes to run simultaneously X-1 basically uses all the threads available this has been tested to ~20 threads. Bottlenecks   could be the size of the strat file and hard drive r/w speed.
    else:
        proc = 1 # simple error check.
    print ("number of processes assigned: "+str(proc))
    with Pool(processes = proc) as p:
        start = time.time()
        async_result = p.map_async(partial(strat_TerrRough_write), xm.xml_list)
        p.close()
        p.join()
        print('Terrain roughness complete in: ' + str(round((time.time() - start)/60, 2)) + ' minutes')


def strat_TerrRough_write(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    #maxSteps = 4 # debugging line with a smaller file
    strat_file=modelh5dir+'/sed.time'+str(maxSteps)+'.hdf5'
    strat=Stratadata()
    strat.loadStratCoords(strat_file)
    
    #just use h5py to load the elev data from the stratal file using the strat class for this loads too many attributes
    df = h5py.File(strat_file, 'r+')
    elev = np.array((df['/layElev']))
    TRI = np.full(elev.shape, -1.0)
    for i in range(0,maxSteps+1):
        elev_temp=(elev[:,i].reshape(strat.ny,strat.nx))# reshape the elevation data as function requires nx / ny array. 
        TRI_temp=np.full(elev_temp.shape, -1.0)
        TRI_temp[1:-1, 1:-1] =(((elev_temp[1:-1, 1:-1])-(elev_temp[:-2, 1:-1]))**2 +  
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[2:, 1:-1]))**2 +  
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[1:-1, :-2]))**2 +
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[1:-1, 2:]))**2 + 
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[2:, 2:]))**2 + 
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[:-2, :-2]))**2 + 
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[2:, :-2]))**2 +
                 ((elev_temp[1:-1, 1:-1]) - (elev_temp[:-2, 2:]))**2)**0.5
        #print(TRI_temp)
        TRI[:,i]=TRI_temp.flatten()
    df.close()    
    #write the results to the strat file
    with h5py.File(strat_file, 'r+') as g:
        if 'layTerrRough' in g.keys():
            del g['layTerrRough'] # if the value exists already remove it
        g.require_dataset('layTerrRough', data=TRI, shape=TRI.shape, dtype="f", compression="gzip" )
        g.close()
        print(str(model.outDir) + '  final strat file written')

#Writes a QGIS compatible TIN mesh ascii file, 
# write the header and then the triangle cell / node index values
#assign an index to the cells as well, this is required for the format
def TIN_hdf_to_2dm(outpath,hdf5_input):
    g=TINfile()
    g.loadTIN(hdf5_input)
    outname=Path(hdf5_input).stem
    #print(outname)
    outfile2dm=outpath+'/'+(Path(hdf5_input).stem)+'.2dm'

    with open (outfile2dm,'w') as f:
            f.write('MESH2D')
            f.write('\n')
            f.write(f"MESHNAME, {outname}")
            f.write('\n')
            for i in range(0,len(g.cells)):
                f.write('E3T')
                f.write(' ')
                f.write(str(i))
                f.write(' ')
                f.write(str(g.cells[i][0]))
                f.write(' ')
                f.write(str(g.cells[i][1]))
                f.write(' ')
                f.write(str(g.cells[i][2]))
                f.write('\n')
            f.close()
    #append coordinate values to the text file (round XY coords to cm / 2 decimal places)
    with open(outfile2dm,'a') as f:
        for i in range(0,len(g.x)):
            f.write('ND')
            f.write(' ')
            f.write(str(i))
            f.write(' ')
            f.write(str(round(g.x[i],2)))
            f.write(' ')
            f.write(str(round(g.y[i],2)))
            f.write(' ')
            f.write(str(round(g.z[i],2)))
            f.write('\n')
        f.close()

def experiment_TIN_hdf_to_2dm(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    model2dmdir=model.outDir+'/2dm'
    isExist = os.path.exists(model2dmdir)
    if not isExist:
        os.makedirs(model2dmdir)
    #loop through the hdf5 tin files creating the output for each
    for i in range(0,maxSteps+1):
        hdf5_input=modelh5dir+'/tin.time'+str(i)+'.hdf5'
        TIN_hdf_to_2dm(model2dmdir,hdf5_input)
        
## multiproc version of the experiment level conversion
def MP_TIN_hdf_to_2dm(modelfile):
#    import lib.badlands_multiproc_run as mpr
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5/'
    maxSteps=int(model.tEnd/model.tDisplay)
    model2dmdir=model.outDir+'/2dm/'
    isExist = os.path.exists(model2dmdir)
    if not isExist:
        os.makedirs(model2dmdir)
    TINlist=[]
    for i in range(0,maxSteps+1):
        TINlist.append(modelh5dir+'tin.time'+str(i)+'.hdf5')
    print (str(len(TINlist))+' surfaces to convert from hdf5 to 2dm for QGIS')
    cpuCount = os.cpu_count() #lets use all the threads
    print ("number of threads available here: "+str(cpuCount)) 
    with Pool(processes = cpuCount) as p:
        async_result = p.map_async(partial(TIN_hdf_to_2dm,model2dmdir),TINlist)
        p.close()
        p.join()
        print('All outputs written')
#################

class Welldata:
    """
    Class for loading data from well csv file.
    """
    def __init__(self):
        self.well_name = None
        self.well_xy = None
    
    def loadWell(self,wellfile):
        """
        Read in the text file for the well location/s, this is an example of the format
        ---------------start file-----------don't use this line
        name,X,Y,thickness,logfile
        Well01,434030,5866930,3500,logfile1.csv
        Well02,514000,5836000,1314,logfile2.csv
        
        -------------end file---------------don't use this line
        """
        well_df = pd.read_csv(wellfile)
        well_arr_in = well_df.to_numpy()
        self.well_xy  = (well_arr_in[:,1:3]) 
        self.well_name= (well_arr_in[:,0])
        self.well_tru = (well_arr_in[:,3]) # truth / total thickness for well
        self.well_log = (well_arr_in[:,4])



class Wellmodel():
    """
    Class for extracting psuedo well data from model.

    """
    def __init__(self):
        self.exp_name=None
        self.well_name = None
        self.well_xy = None

    def extractWellsExperimentTIN(self,modelfile,wellfile,outfile_loc,srch_dist='nearest'): #was extractWellsModel
        """Get the well locations from the wellfile then extract 
        Then extract the properties from the TIN file at that location
        Setting srch_dist to something other than NULL will return all (or none) of the daa points within a radius.
        """
        h5file=h5py.File(str(outfile_loc),'a')
        
        model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
        modelh5dir=model.outDir+'/h5'
        maxSteps=int(model.tEnd/model.tDisplay)
        tin_file=modelh5dir+'/tin.time'+str(maxSteps)+'.hdf5'
        
        #get the experiment parameters that we think are needed
        self.exp_name = str(ntpath.basename(model.outDir))
        self.total_layers=maxSteps
        self.tDisplay=model.tDisplay
        self.tStart=model.tStart
        self.tEnd=model.tEnd
        self.poro0=model.poro0
        self.poroC=model.poroC
        
        wells=Welldata()
        wells.loadWell(wellfile)
        #self.well_name=[]
        
        tinData=TINfile()
        print('load data from: '+tin_file)
        tinData.loadTIN(tin_file)
        tingrid = pd.DataFrame({'X':tinData.x.flatten(),'Y':tinData.y.flatten()}).to_numpy() 
        
        # make the experiment group in the h5file
        h5file.create_group(self.exp_name)
        # add some useful experiment attributes to the hdf5
        en=str(self.exp_name)+'/'
        h5file[en].attrs['tDisplay']=self.tDisplay
        h5file[en].attrs['total_layers']=self.total_layers
        h5file[en].attrs['tStart']=self.tStart
        h5file[en].attrs['tEnd']=self.tEnd
        h5file[en].attrs['poro0']=self.poro0
        h5file[en].attrs['poroC']=self.poroC    
        
        #Loop through each psuedo well location in this model, find the nearest node in the stratgrid, extract and append data
        for i in range(0,(len(wells.well_xy))):
            tree = spatial.cKDTree([wells.well_xy[i]]) # build KD tree from well location
            dd, ii = tree.query(tingrid, k=1) # calculate distance from well location to all points in the sedfile
            if srch_dist == 'nearest':
                well_nodes=np.unique(np.where(dd==np.min(dd))) # select the nearest node only by smallest dd value
            if type(srch_dist) == float or type(srch_dist) ==int:
                well_nodes=np.unique(np.where(dd<=srch_dist)) # alternate selection using search radius, returns multiple values.
            if type(srch_dist) != float and type(srch_dist) !=int and srch_dist != 'nearest':
                print('srch_dist should only be blank/empty or a number')            
            #make the well sub-group down for the well things
            enwn=str(self.exp_name)+'/'+str(wells.well_name[i])+'/'
            h5file.create_group(enwn)
                       
            # just grab all of the attributes from the tin_attrs list.
            for j in tinData.tin_keylist:
               h5file.create_dataset(str(enwn)+'/'+j,data=tinData.tin_layer[j][well_nodes])

        print (str(len(wells.well_name)) + ' psuedo well locations added to results')
        h5file.close()
        return('data written to '+str(outfile_loc))
        
        
    
    def extractWellsExperiment(self,modelfile,wellfile,outfile_loc): #was extractWellsModel
        """
        Get the well locations from the wellfile then extract the 
        Then extract the properties at the well
        """
        h5file=h5py.File(str(outfile_loc),'a')
        
        model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
        modelh5dir=model.outDir+'/h5'
        maxSteps=int(model.tEnd/model.tDisplay)
        strat_file=modelh5dir+'/sed.time'+str(maxSteps)+'.hdf5'
        
        #get the experiment parameters that we think are needed
        self.exp_name = str(ntpath.basename(model.outDir))
        self.total_layers=maxSteps
        self.tDisplay=model.tDisplay
        self.tStart=model.tStart
        self.tEnd=model.tEnd
        self.poro0=model.poro0
        self.poroC=model.poroC
        
        wells=Welldata()
        wells.loadWell(wellfile)
        #self.well_name=[]

        # load stratigraphy file
        strat=Stratadata()
        print('load data from: '+strat_file)
        strat.loadStrat(strat_file)
        stratgrid = pd.DataFrame({'X':strat.x.flatten(),'Y':strat.y.flatten()}).to_numpy() # format sed.time input xy input for kdtree
        
        # make the experiment group in the h5file
        h5file.create_group(self.exp_name)
        # add some useful experiment attributes to the hdf5
        en=str(self.exp_name)+'/'
        h5file[en].attrs['tDisplay']=self.tDisplay
        h5file[en].attrs['total_layers']=self.total_layers
        h5file[en].attrs['tStart']=self.tStart
        h5file[en].attrs['tEnd']=self.tEnd
        h5file[en].attrs['poro0']=self.poro0
        h5file[en].attrs['poroC']=self.poroC       
        
        #Loop through each psuedo well location in this model, find the nearest node in the stratgrid, extract and append data
        for i in range(0,(len(wells.well_xy))):
            # print('extract data at: ' + wells.well_name[i]+ ' location') # debugging line remove later
            tree = spatial.cKDTree([wells.well_xy[i]]) # build KD tree from well location
            dd, ii = tree.query(stratgrid, k=1) # calculate distance from well location to all points in the sedfile
            well_nodes=np.unique(np.where(dd==np.min(dd))) # select the nearest node only by smallest dd value
            #well_nodes=np.unique(np.where(dd<=srch_dist)) # alternate selection using search radius, returns multiple values then reduce to 1 using average in append component.
            
            #make the well sub-group down for the well things
            enwn=str(self.exp_name)+'/'+str(wells.well_name[i])+'/'
            h5file.create_group(enwn)
                       
            # just grab all of the attributes from the strat_attrs list.
            for j in strat.strat_attribs:
               h5file.create_dataset(str(enwn)+'/'+j,data=strat.strat_layers[j][well_nodes])

        print (str(len(wells.well_name)) + ' psuedo well locations added to results')
        h5file.close()
        return('data written to '+str(outfile_loc))

def extract_WellLayer_all(xml_dir,wellfile,outfile='experiment_well_extracts',outdir='results_final/'):
    isExist = os.path.exists(outdir)
    if not isExist:
        os.makedirs(outdir)
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
    wells=Welldata()
    wells.loadWell(wellfile)
    modelList=xm.xml_list
        
    for j in xm.xml_list:
#        print ('extract data from '+j)
        exp_loop=Wellmodel()
        exp_loop.extractWellsExperiment(j,wellfile,str(outdir+outfile))

def extract_WellThick_ALL_TIN(xml_dir,wellfile,outfile='experiment_well_extracts',outdir='results_final/',srch_dist='nearest'):
    isExist = os.path.exists(outdir)
    if not isExist:
        os.makedirs(outdir)
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
    wells=Welldata()
    wells.loadWell(wellfile)
    modelList=xm.xml_list
        
    for j in xm.xml_list:
#        print ('extract data from '+j)
        exp_loop=Wellmodel()
        exp_loop.extractWellsExperimentTIN(j,wellfile,str(outdir+outfile),srch_dist)



"""
Below are a set of functions that determine the upstream fetch areas from a set of starting points, nodes near a lake shoreline are what I need,but the same thing should work for an ocean coastline.
I initially created this as a sub-basin ID (thus the sbid references) but fetch is a much better name
so this is how the outputs are written to the flow file.
In the case it was initially set up for, additional complexity is created by some flow paths running sub parallel to the basin boundary and multipe checks are made to ensure that all nodes are classified.
This might be simplified and run faster cases where not all nodes are required or a closed polygon around a depocenter is not used.
"""

class Fetchnodes:
    """
    Class for loading data from fetch area node initial points csv file.
    """
    def __init__(self):
        self.fnodes = None
    
    def loadFnodes(self,nodefile):
        """
        Read in the text file for the node locations, this is an example of the format, 
        NodeID should ne numeric, preferably an int
        
        X,Y NodeID
        ---------------start file-----------don't use this line
        526046,5856329,9
        528074,5855388,10
        530102,5854446,11
        -------------end file---------------don't use this line
        """
        self.fnodes=np.loadtxt(open(nodefile), delimiter=",")
        return (self)


# calculate the area for all steps in a set of multiple experiments    
def fetch_area_doegen(xml_dir,fetchnodefile,interp_dist=None):
    xm=mpr.XmlList()
    xm.loadXml(xml_dir)
    for f in xm.xml_list:
        fetchrun=fetch_area_experiment(f,fetchnodefile,interp_dist)    
    return(fetchrun)    

#calculate the fetch area for all steps i a single experiment 
def fetch_area_experiment(modelfile,fetchnodefile, interp_dist=None):
    fnds=Fetchnodes()
    fnds.loadFnodes(fetchnodefile)
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    if interp_dist==None:
        interp_dist=model.stratdx*2 # the distance to limit flow node search near the fetch node initial points. Set here to double the model resolution, anything over 1 should* be fine.
    
    # generate list of hdf5 flowfiles
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    fflst=[]
    proc = os.cpu_count()-1
    
    ##or just os.listdir for flow.time** files?
    for i in range(0,maxSteps+1):
        flow_file=modelh5dir+'/flow.time'+str(i)+'.hdf5'
        fflst.append(flow_file)

    print('Determining fetch areas for '+model.outDir )   
   #print ('Starting parallel calucations using multiprocessing:')
    with Pool(processes = proc) as p:
        async_result = p.map_async(partial(fetch_area_calculate,fnds,interp_dist),fflst)
        p.close()
        p.join()
    
    print('Fetch area calcs complete for '+model.outDir )
    return (async_result) 

#determine the fetch areas for a single flowfile (flow.timeXX.hdf5)
# To account for loops in the flow paths that are sub-parallel to the boundary of the depocenter the oints within the depocetner are removed for the initlal calculations. If your experiment/model does not contain those 
# you may find that a line describing a shoreline may work better, but this isn't implemented here.
def fetch_area_calculate(fnodes_arr,interp_dist,flow_file,depth=10):

    flw=Flowfile()
    flw.loadFlow(flow_file)
    fnodes=fnodes_arr.fnodes
#    fnodes=np.loadtxt(open(nodefile), delimiter=",")
    connect = flw.connect #-1 to adjust indices--see hydroGrid.py in companion
    basin=flw.basin
    fetch=(np.full_like(basin,-9999)) 
  
    x=flw.x
    y=flw.y
#networkX creation
    G = nx.DiGraph() #create initial newtorkx graph - data as directional to allow upstream
    GN = nx.Graph() # create a non-directed Graph of the same dataset
    G.clear()
    G.add_edges_from(connect) # build the edges /connections from the flow connections
    GN.add_edges_from(connect)
    H = G.copy()
    
# remove nodes from the graph that are within the depocenter using matplotlib 
    ptsdf=pd.DataFrame()
    ptsdf['x']=x
    ptsdf['y']=y
    pts=ptsdf[['x','y']]
    
    limitline=plt.plot(fnodes[:,0],fnodes[:,1])
    mask = limitline[0].get_path().contains_points(pts)
    inside=np.where(mask)[0].tolist()
    plt.close()
    
#   remove the interior/inside points from the G graph
    G.remove_nodes_from(inside)
    
# KDtree for distance to limits
    tree = spatial.cKDTree(fnodes[0:, 0:2])
    flowgrid = pd.DataFrame({'X':x.flatten(),'Y':y.flatten()}).to_numpy()
    dd, ii = tree.query(flowgrid, k=1) #dd should be the distance
    def_limit_nodes=np.unique(np.where(dd<interp_dist)) # initial nodes on the flow grid to feed the filters
    node_num=len(def_limit_nodes)
#interpolate the fnode values onto sub-basin array
    dlxy=(fnodes[0:, 0:2])
    dldist=(fnodes[0:, 2])
    sbid = interp.griddata(dlxy, \
                       dldist, \
                       flowgrid, \
                       method='nearest')

    
    #create the list of nodes to work with  
    
    filt_list= def_limit_nodes.tolist()
# tag the nodes within the depocenter
    for i in range(0,len(filt_list)):
         if filt_list[i] in inside:
            filt_list[i]=1
 
# make a new list with the tagged nodes removed
    shortlist=list.copy(filt_list)
    shortlist = [i for i in shortlist if i != 1]

#run the upstream node function and update the fetch ID with the sbid at downstream source on each loop
    for i in shortlist:
        nd=upstream_node(i,G)
        for f in nd:
            fetch[f]=sbid[i] #finds all upstream nodes from limits (erosional domain if you've selected a depocenter)

#create new shortlist from unique upstream fetch areas
#removed as there were some issues with this returning an index that isn't tracking well. Fix later for minor speed improvement.
#    uniq_upnodes=np.unique(fetch,return_index=True)[1].tolist()


#fill in the remaining nodes inside the depocenter
    for i in shortlist:
        ds=downstream_node(i,H)
        for j in ds:
            if j in inside:
                fetch[j]=-1*sbid[i]
    for i in inside:
       if fetch[i]==-9999:
            fetch[i]=0        

# write catchment ID back to hdf5
    with h5py.File(flow_file,'r+') as h:
        if 'fetch' in h.keys():
            del h['fetch'] # if the value exists already remove it
        h.require_dataset('fetch', data=fetch, shape=fetch.shape, dtype="int32", compression="gzip")
        #print (flow_file+" modified")
    h.close()
  
    return(shortlist)
    

#  networkX  finds upstream nodes. Starting point is the list of nodes near basin limits interpolated using the kdtree.
def upstream_node(dsinit,G): #initial downstream node, catchment number, network
    alist=[]
    node=list(nx.edge_dfs(G,dsinit, orientation='reverse'))
    if len(node)==0 or len(node)==1:
        return(1,1) #send something back if there's no upstream nodes just to stop errors in case filters were bad
    else:
        alist=np.array(list(node))
    arr=np.append((np.array(alist)[0:,0]),(np.array(alist)[0:,1])) # get 0-1 columns append as single column
    numbers=list(map(int, arr)) # array / also list of connections, remove duplicates
# int conversion as output is list and string
    upnodes=np.unique(numbers) # output is
    return(upnodes)
    
#  networkX  finds downsteam nodes. Starting point is the list of nodes near basin limits.    
def downstream_node(dsinit,H): #initial downstream node, networkx unmasked array
    alist=[]
    node=list(nx.edge_dfs(H,dsinit))
    if len(node)<=1:
        return(1,1) #send something back if there's no downstream nodes just to stop errors in case filters were bad
    else:
        alist=np.array(list(node))
    arr=np.append((np.array(alist)[0:,0]),(np.array(alist)[0:,1])) # get 0-1 columns append as single column
    numbers=list(map(int, arr)) # array / int conversion as output is list and string
    downnodes=np.unique(numbers) # output is also list of connections, remove duplicates
    return(downnodes) 

## modified version of badlands-companions vtk/vts creator.
## reshapes rather than iterates to create arrays with correct data and shapes.
## creates properties in the output of the properties in the sed.time file.
def gridtovts(stratfile,vtsoutfile):
    strat=Stratadata()
    strat.loadStrat(stratfile)
    
    #reshape the xy arrays to the grid dimensions, then use np.repeat to make an xy for each layer, reshape z
    # this is more efficient than layered iterating
    xrs=strat.x.reshape(strat.ny,strat.nx)
    yrs=strat.y.reshape(strat.ny,strat.nx)
    x=np.repeat(xrs[:, :, np.newaxis],strat.nz,axis=2)
    y=np.repeat(yrs[:, :, np.newaxis],strat.nz,axis=2)
    z = (strat.strat_layers['layDepth'][:]).reshape(strat.ny,strat.nx,strat.nz)
    
    # Create the dictionary for the point attributes data and reshape the inputs to fit.
    layers={}
    for i in strat.strat_attribs:
        if 'layDepth' not in i:
            layers[i] = (strat.strat_layers[i][:]).reshape(strat.ny,strat.nx,strat.nz)
    
    #pass pyevtk gridToVTK function the data it needs
    gridToVTK(vtsoutfile, x, y, z, pointData = layers)


