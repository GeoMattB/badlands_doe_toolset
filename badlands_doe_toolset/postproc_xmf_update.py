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
        tmstep=float(i*model.tDisplay)
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
        tmstep=float(i*model.tDisplay)
        cells=len(f['cells'])
        geom=len(f['coords'])
#        row=pd.DataFrame()
#        row['cells']=cells,
#        row['geom']=geom,
        texfile=open(xmf_file,'w')
        texfile.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
 <Domain>
    <Grid GridType="Collection" CollectionType="Spatial">
      <Time Type="Single" Value="{tmstep}"/>
      <Grid Name="Block.0">
         <Topology Type="Triangle" NumberOfElements="{(int(row['cells']))}" BaseOffset="1">
          <DataItem Format="HDF" DataType="Int" Dimensions="{int(cells)} 2">h5/tin.time{i}.hdf5:/cells</DataItem>
         </Topology>
         <Geometry Type="XYZ">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 3">h5/tin.time{i}.hdf5:/coords</DataItem>
         </Geometry>
         <Attribute Type="Scalar" Center="Node" Name="lake">
          <DataItem Format="HDF" NumberType="Integer" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/lake</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="Discharge [m3/s]">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/discharge</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="cumdiff">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumdiff</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="cumhill [m3/s]">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumhill</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="cumfail adim">
          <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time{i}.hdf5:/cumfail</DataItem>
         </Attribute>
         <Attribute Type="Scalar" Center="Node" Name="Sealevel">
          <DataItem ItemType="Function" Function="$0 * 0.00000000001 + {SeaLvl}" Dimensions="{int(geom)} 1">
           <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="{int(geom)} 1">h5/tin.time1.hdf5:/cumdiff</DataItem>
          </DataItem>
         </Attribute>
         """)
        for g in (f.keys()): ## if there are any new parameters add them to the xmf for paraview.
            if g !='coords' and g !='cells' and g !='lake' and g !='discharge' and g !='cumdiff' and g !='cumhill' and g !='cumhill' and g !='cumfail' and g !='lake':
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
    
