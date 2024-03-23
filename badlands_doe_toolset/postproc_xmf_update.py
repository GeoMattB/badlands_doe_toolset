##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file has been added to the Badlands companion by Matt Boyd.                 ##
##                                                                                   ##
##  For full license and copyright information contact the author.                  ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This recreates the xmf files in the flow.time and tin.time series that are required for display in paraview.
It recreates the standard badlands outputs to remain compatible with existing paraview / xmf outputs 
by reading the hdf5 files then adds any new attributes it finds.
It will overwite / create a new xmf file so if you've lost your xmf files this could possibly recreate them.
"""

import h5py
import pandas as pd
from badlands import xmlParser

def flowfile_xmf(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    modelXmfdir=model.outDir+'/xmf'
    for i in range(0,maxSteps+1):
        hdf_file=modelh5dir+'/flow.time'+str(i)+'.hdf5'
        print ('reading '+hdf_file)
        xmf_file=str(modelXmfdir)+'/flow.time'+str(i)+'.xmf'
        print ('output will be '+xmf_file)
        f=h5py.File(hdf_file, 'r')
        tmstep=float(model.tStart+(i*model.tDisplay))
        conshp=len(f['connect'])
        geom=len(f['coords'])
#        row=pd.DataFrame()
#        row['conshp']=conshp,
#        row['geom']=geom,
        texfile=open(xmf_file,'w')
        texfile.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
 <Domain>
    <Grid GridType="Collection" CollectionType="Spatial">
      <Time Type="Single" Value="{tmstep}"/>
      <Grid Name="Block.0">
         <Topology Type="Polyline" NodesPerElement="2" NumberOfElements="{int(conshp)}" BaseOffset="1">
          <DataItem Format="HDF" DataType="Int" Dimensions="{int(conshp)} 2">h5/flow.time{i}.hdf5:/connect</DataItem>
         </Topology>
         <Geometry Type="XYZ">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 3">h5/flow.time{i}.hdf5:/coords</DataItem>
         </Geometry>
         <Attribute Type="Scalar" Center="Node" Name="BasinID">
          <DataItem Format="HDF" NumberType="Integer" Precision="4" Dimensions="{int(geom)} 1">h5/flow.time{i}.hdf5:/basin</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="Discharge [m3/s]">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/flow.time{i}.hdf5:/discharge</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="Chi">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/flow.time{i}.hdf5:/chi</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="sedload [m3/s]">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/flow.time{i}.hdf5:/sedload</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="flowdensity adim">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/flow.time{i}.hdf5:/flowdensity</DataItem>
         </Attribute>""")
        for g in (f.keys()): ## if there are any new parameters add them to the xmf for paraview.
            if g !='coords' and g !='connect' and g !='basin' and g !='discharge' and g !='chi' and g !='sedload' and g !='flowdensity':
                texfile.write(f"""
         <Attribute Type="Scalar" Center="Node" Name="{str(g)}">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/flow.time{i}.hdf5:/{str(g)}</DataItem>
         </Attribute>""")
        texfile.write(f""" 
      </Grid>
    </Grid>
 </Domain>
</Xdmf>
""")
    texfile.close()
    
def tinfile_xmf(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    modelXmfdir=model.outDir+'/xmf'
    SeaLvl=model.seapos
    for i in range(0,maxSteps+1):
        hdf_file=modelh5dir+'/tin.time'+str(i)+'.hdf5'
        print ('reading '+hdf_file)
        xmf_file=str(modelXmfdir)+'/tin.time'+str(i)+'.xmf'
        print ('output will be '+xmf_file)
        f=h5py.File(hdf_file, 'r')
        tmstep=float(model.tStart+(i*model.tDisplay))
        cells=len(f['cells'])
        geom=len(f['coords'])
        row=pd.DataFrame()
        row['cells']=cells,
        row['geom']=geom,
        texfile=open(xmf_file,'w')
        texfile.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
 <Domain>
    <Grid GridType="Collection" CollectionType="Spatial">
      <Time Type="Single" Value="{tmstep}"/>
      <Grid Name="Block.0">
         <Topology Type="Triangle" NumberOfElements="{(int(row['cells']))}" BaseOffset="1">
          <DataItem Format="HDF" DataType="Int" Dimensions="{int(cells)} 3">h5/tin.time{i}.hdf5:/cells</DataItem>
         </Topology>
         <Geometry Type="XYZ">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 3">h5/tin.time{i}.hdf5:/coords</DataItem>
         </Geometry>
         <Attribute Type="Scalar" Center="Node" Name="lake">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/lake</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="Discharge">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/discharge</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="EroDep">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumdiff</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="EroDep hillslope">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumhill</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="EroDep failure">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumfail</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="Sealevel">
          <DataItem ItemType="Function" Function="$0 * 0.00000000001 + {SeaLvl}" Dimensions="{int(geom)} 1">
           <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumdiff</DataItem>
          </DataItem>
         </Attribute>
         """)
        for g in (f.keys()): ## if there are any new parameters add them to the xmf for paraview.
            if g !='coords' and g !='area' and g !='cells' and g !='lake' and g !='discharge' and g !='cumdiff' and g !='cumhill' and g !='cumfail' and g !='lake':
                texfile.write(f"""
         <Attribute Type="Scalar" Center="Node" Name="{str(g)}">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/{str(g)}</DataItem>
         </Attribute>""")
        texfile.write(f""" 
      </Grid>
    </Grid>
 </Domain>
</Xdmf>
""")
    texfile.close()
 
# XMF for sed time strat files
# single xmf build write the xmf to view in paraview
def sedfile_xmf(hdf_file,modelinc,tDisplay=100000):  
        print ('load '+hdf_file)
        strat=ppu.Stratadata()
        strat.loadStrat(hdf_file)
        f=h5py.File(hdf_file, 'r')
        tmstep=modelinc*tDisplay
        NumberOfElements=str(strat.ny)+' ' +str(strat.nx) + ' ' +str(strat.nz)
        Dimensions= str(len(strat.x))+' '+ str(strat.nz)
        print(Dimensions)
        xmf_file='xmf/sed.time'+str(modelinc)+'.xmf'
        print ('output will be '+xmf_file)
        
        texfile=open(xmf_file,'w')
        texfile.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
 <Domain>
  <Grid GridType="Collection" CollectionType="Spatial">
  <Time Type="Single" Value="{tmstep}"/>
    <Grid Name="Block.0" GridType="Uniform">
     <Topology TopologyType="3DSMesh" Format="HDF" NumberOfElements="{NumberOfElements}" BaseOffset="1">
     </Topology>
      <Geometry GeometryType="X_Y_Z">
        <DataItem Name="X" Dimensions="{Dimensions}" NumberType="Float" Precision="4" Format="HDF"> {hdf_file}:/X </DataItem>
        <DataItem Name="Y" Dimensions="{Dimensions}" NumberType="Float" Precision="4" Format="HDF"> {hdf_file}:/Y </DataItem>
        <DataItem Name="Z" Dimensions="{Dimensions}" NumberType="Float" Precision="4" Format="HDF"> {hdf_file}:/layDepth </DataItem>
      </Geometry>
""")
        for g in (f.keys()): ## if there are any other attributes add them to the xmf for paraview.
            if g !='X' and g !='Y' and g !='coords':
                texfile.write(f"""
      <Attribute Type="Scalar" Center="Node" Name="{str(g)}">
        <DataItem Dimensions="{NumberOfElements}" Format="HDF" NumberType="Float" Precision="4"> {hdf_file}:/{str(g)} </DataItem>
      </Attribute>""")
        texfile.write(f""" 
    
    </Grid>
   </Grid>
 </Domain>
</Xdmf>
""")
        texfile.close()


#write all of the xmf files and then the XDMF for a time series. 
def sedfile_xdmf(modelfile):
    model=xmlParser.xmlParser(modelfile,'False') #just read in the xml 'False' stops it writing output directories'
    modelh5dir=model.outDir+'/h5'
    maxSteps=int(model.tEnd/model.tDisplay)
    modelXmfdir=model.outDir+'/xmf'
    SeaLvl=model.seapos
    
    ## send to xmf function
    for i in range(1,maxSteps+1):
        hdf_file=modelh5dir+'/sed.time'+str(i)+'.hdf5'
        sedfile_xmf(hdf_file,i)
  
##XMDF part
#Write the xdmf time series file to tie this all together
    xdmf_file=model.outDir+'/sedtime_series.xdmf'
    texfile=open(xdmf_file,'w')
    texfile.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
 <Domain>
    <Grid GridType="Collection" CollectionType="Temporal">
""")
    for i in range(1,maxSteps+1):
        xmf_file=str(modelXmfdir)+'/sed.time'+str(i)+'.xmf'
        texfile.write(f"""     <xi:include href="xmf/sed.time{i}.xmf" xpointer="xpointer(//Xdmf/Domain/Grid)"/>
""")
    
    texfile.write(f"""  </Grid>
 </Domain>
</Xdmf>
""")
    texfile.close()   
    
    

    
    