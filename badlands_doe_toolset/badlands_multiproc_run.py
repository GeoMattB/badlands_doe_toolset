from multiprocessing import Pool
import multiprocessing
import threading
import os
import sys
import time
from contextlib import redirect_stdout
import ntpath
import time
import re

from badlands.model import Model as badlandsModel
from badlands.model import xmlParser



# class to load all of the model xml files generated by doegen in the xml directory to a list that is sent to the multiprocessor
class XmlList:
        def __init__(self):
            self.xml_list = []
        def loadXml(self,xml_dir):
            for f in os.listdir(xml_dir):
                if f.endswith(".xml"):
                    self.xml_list.append(xml_dir+'/'+f)
            self.xml_list.sort(key=lambda xml_file : list(map(int, re.findall(r'\d+', xml_file)))[0])
            # Here we do a simple check that all the xml experiment configuration files can be read by the XmlParser without errors BEFORE we do the multiprocessing.                   
            print ('Checking XML files are valid. If there is an error after this then there is something wrong with your configuration xml.')
            for i in self.xml_list:
                test=xmlParser.xmlParser(i,'false')
                #print (i + ' out directory is '+ test.outDir)
            test=None    
            print ('All good, simple xml check passed, number of experiment configuration files loaded: '+ str(len(self.xml_list)))
           

def badlands_experiment(xml_file):
    starting.acquire() # lock to only start 1 experiment process every -see next line- seconds
    threading.Timer(0.2, starting.release).start() # release so the next of the processes can start.
    print("Processing : ", str(xml_file))
    g=xmlParser.xmlParser(xml_file, 'false') # read the xml file to tag names and locations 'false' exploits a bug, needs to be fixed
    model = badlandsModel()
    model.load_xml(xml_file)
    if not os.path.exists('log/'): 
        os.makedirs('log/')
    with open("log/"+str(ntpath.basename(g.outDir))+'.log', 'w') as f: #redirect the print progress outputs to a file
        with redirect_stdout(f):
            rstart = time.time()
            print(rstart)
            model.run_to_time(g.tEnd)
            print('Runtime:, ' + str(round((time.time() - rstart)/60, 2)) + ' ,minutes')
        f.close()
    print(str(xml_file)+' - complete,   Runtime: ' + str(round((time.time() - rstart)/60, 2)) + ' minutes')
    return (str(xml_file)+' Runtime:,' + str(round((time.time() - rstart)/60, 2)) + ',minutes')

def init(lock):    # set up the lock function so there's a delay between starting each badlands experiment.
    global starting
    starting = lock

def badlands_multi_experiment(xml_dir,proc=None):
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount)) 
    if proc==None:
        if cpuCount >= 2:
            proc = int(cpuCount/2) #Number of Processes to run simultaneously X/2 is conservative and based on an oldish PC, HPC and modern CPU's require more testing.
        else:
            proc = 1 #default for safety
    print ("number of processes assigned: "+str(proc))
    xm=XmlList()
    xm.loadXml(xml_dir)
    #print (xm.xml_list)  
    print ('Starting parallel experiment runs using multiprocessing:')
    with Pool(processes = proc,initializer=init, initargs=[multiprocessing.Lock()]) as p:
        start = time.time()
        #p.map_async(badlands_experiment, xm.xml_list)
        async_result = p.map_async(badlands_experiment, xm.xml_list) #the multiproc start bit
        p.close()
        p.join()
        print('Simulations completed; total processing time: ' + str(round((time.time() - start)/60, 2)) + ' minutes')
        return (async_result)
    

# slightly extended function to just run tests without outputs.
def testXml(xml_dir):
    xml_list=[]
    for f in os.listdir(xml_dir):
        if f.endswith(".xml"):
            xml_list.append(f)
    # Here we do a simple check that the xml configuration can be read by the XmlParser without errors.                  
    print ('Checking XML files are valid. If there is an error after this then there is something wrong with your configuration xml.')
    for i in xml_list:
        print ('checking : '+i)
        test=xmlParser.xmlParser(xml_dir+'/'+i,'false')
        print ('passed : ' + i + ' out directory is '+ test.outDir)

### This is essentially the debugging part but it might come in useful.
def badlands_multi_test(xml_dir,proc=None):
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount)) 
    if proc==None:
        if cpuCount >= 2:
            proc = int(cpuCount/2) #Number of Processes to run simultaneously X/2 is conservative and based on an oldish PC, HPC and modern CPU's require more testing.
        else:
            proc = 1 #default for safety
    print ("number of processes assigned: "+str(proc))
    xm=XmlList()
    xm.loadXml(xml_dir)
    #print (xm.xml_list)  
    print ('Starting parallel experiment runs using multiprocessing:')
    with Pool(processes = proc,initializer=init, initargs=[multiprocessing.Lock()]) as p:
        start = time.time()
        #p.map_async(badlands_experiment, xm.xml_list)
        async_result = p.map_async(testXmlFile, xm.xml_list) #the multiproc start bit
        p.close()
        p.join()
        print('Simulations completed; total processing time: ' + str(round((time.time() - start)/60, 2)) + ' minutes')
        return (async_result)

def testXmlFile(xml_file):
    starting.acquire() # lock to only start 1 experiment process every -see next line- seconds
    threading.Timer(0.2, starting.release).start() # release so the next of the processes can start.
    rstart = time.time()
    print("testing : ", str(xml_file))
    #time.sleep(0.2)
    return(xml_file+str(' , ')+str(rstart-time.time()))
    
    
#### Tests to estimate how much faster multiproc will run multiple models.
def badlands_optimise_multiproc(xml_dir,procStep=2):
    cpuCount = os.cpu_count()
    print ("number of threads available here: "+str(cpuCount))
    speedtest=[]
    for i in range(2,cpuCount+procStep,procStep):
        print('test proc no.: '+str(i))
        exptest=badlands_multi_experiment_optimise(xml_dir,i)
        speedtest.append(exptest)
    return(speedtest)

def badlands_multi_experiment_optimise(xml_dir,proc):
    xm=XmlList()
    xm.loadXml(xml_dir)
    #print (xm.xml_list)  
    print ('Starting parallel experiment test runs using multiprocessing: '+str(proc)+str(' processes'))
    with Pool(processes = proc,initializer=init, initargs=[multiprocessing.Lock()]) as p:
        start = time.time()
        #p.map_async(badlands_experiment, xm.xml_list)
        async_result = p.map_async(badlands_experiment_optimise, xm.xml_list) #the multiproc start bit
        p.close()
        p.join()
        print('Simulations completed; total processing time: ' + str(round((time.time() - start)/60, 2)) + ' minutes')
        return (proc,(round((time.time() - start)/60, 2)))
    
def badlands_experiment_optimise(xml_file):
    starting.acquire() # lock to only start 1 experiment process every -see next line- seconds
    threading.Timer(0.1, starting.release).start() # release time in seconds so the next of the processes can start.
    # print("Processing : ", str(xml_file))
    g=xmlParser.xmlParser(xml_file, 'false') # read the xml file to tag names and locations 'false' exploits a bug, needs to be fixed
    model = badlandsModel()
    model.load_xml(xml_file)
    if not os.path.exists('log_optim/'): 
        os.makedirs('log_optim/')
    with open("log_optim/"+str(ntpath.basename(g.outDir))+'.log', 'w') as f: #redirect the print progress outputs to a file
        with redirect_stdout(f):
            rstart = time.time()
            print(rstart)
            model.run_to_time(g.tDisplay*2)
            print('Runtime:, ' + str(round((time.time() - rstart)/60, 2)) + ' ,minutes')
        f.close()
    #print(str(xml_file)+' - complete,   Runtime: ' + str(round((time.time() - rstart)/60, 2)) + ' minutes')
    return (str(xml_file)+' Runtime:,' + str(round((time.time() - rstart)/60, 2)) + ',minutes')    