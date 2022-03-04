"""
Ingests DoEgen .csv output to produce a series of Badlands .xml experiment configuration files.
"""

# TODO: These imports seem a bit excessive; go back through and remove unused ones
import argparse
import xml.etree.ElementTree as etree
from xml.etree import ElementTree
import pandas as pd 
from pathlib import Path
import os
import sys

"""Constants and Environment"""

ap = argparse.ArgumentParser()
ap.add_argument("-f","--file", default = "Designtable_optimal.csv", type=Path)
ap.add_argument("-i","--input", default = "./design_results", type=Path)
ap.add_argument("-o","--output", default = "./badlands_xml", type=Path)
ap.add_argument("-s","--scenarios", default = "./scenarios", type=Path)
args = ap.parse_args()

"""
Comment friendly XML builder
From stack: https://stackoverflow.com/a/34324359
TODO: Switch lxml for xml parsing? comments are supported natively!
"""
class CommentedTreeBuilder(ElementTree.TreeBuilder):
    def comment(self, data):
        self.start(ElementTree.Comment, {})
        self.data(data)
        self.end(ElementTree.Comment)

nsmap = {
    'xsi':"http://www.w3.org/2001/XMLSchema-instance"
}

"""
Define scenario dicts
TODO: build functions to generate these xml snippets based on user supplied inputs.
"""

def read_scenarios(scenario_type):
    scenarios = Path(args.scenarios).glob(scenario_type)
    dict = {}
    for xml in scenarios:
        with open(xml) as f:
            scenario = f.read()
        dict[str(os.path.basename(xml))] = scenario
    return dict

def badlands_encode(row, *args):
    
    """
    Function for writing an xml file in badlands format from a pandas df row.
    Current logic depends on df column names to populate the xml file
    The first set of if statements checks to see if the csv (generated from the DoEgen outputs) has specified any complex
    structures. That is, multiple, variable precipitation, tectonic or erodibility sub- xmls. If so, it loads them as the sting for that section
    of the configuration.
    """
# make the dictionaries where there are options.
    erodability_dict = read_scenarios('erodability/*.xml')
    precipitation_dict = read_scenarios('precipitation/*.xml')
    tectonic_dict = read_scenarios('tectonic/*.xml')
    carbonate_dict = read_scenarios('carbonate/*.xml')

# Load some optional default values
    if "mindt" in row.keys():
        mindt=row['mindt']
    else: 
        mindt=1

    if "maxdt" in row.keys():
        maxdt=row['maxdt']
    else:
        maxdt=row['display']
    
    if "position" in row.keys():
        position=row['position']
    else:
        position=0    

# Sea level curve
    if "curve" in row.keys():
        curve=f"""<!--  Optional sea-level curve file  -->
        <curve>data/{row['curve']}</curve>"""
    else:
        curve= '<!-- no optional sea level curve defined -->'


# Tectonic structure.
    if 'Tectonic' in row.keys():
        try:
            tectonic_structure = tectonic_dict[row['Tectonic']]
        except(KeyError): # this is the error returned if the specified xml file is not found.
            print('No Tectonic structure files found')
    elif 'dstart' in row.keys():
        tectonic_structure=f"""    <tectonic>
        <!-- Number of tectonic events -->
		<events>1</events>
		<!-- Displacement definition -->
		<disp>
			<!-- Displacement start time [a] -->
			<dstart>{row['dstart']}</dstart>
			<!-- Displacement end time [a] -->
			<dend>{row['dend']}</dend>
			<!-- Displacement map [m] -->
			<dfile>{row['dfile']}</dfile>
        </disp>
    </tectonic>"""
    else:
        tectonic_structure = '<!-- no tectonics defined -->'

# Precipitation structure
    if '.xml' in str(row['Precipitation']):
        try:
            precipitation_structure = precipitation_dict[row['Precipitation']]
        except(KeyError): # this is the error returned if the specified xml file is not found.
            print('No precipitation structure files found using single event climate values')
    else:
        precipitation_structure = f"""<!-- Precipitation structure -->
    <precipitation>
        <!-- Number of precipitation events -->
        <climates>1</climates>
        <!-- Precipitation definition -->
        <rain>
            <!-- Rain start time [a] -->
            <rstart>{row['rstart']}</rstart>
            <!-- Rain end time [a] -->
            <rend>{row['rend']}</rend>
            <!-- Precipitation value [m/a] - (optional) -->
            <rval>{row['Precipitation']}</rval>
        </rain>
    </precipitation>"""

## Stream power law optional parameters:
    if "diffnb" in row.keys():
        diffnb = f"""<diffnb>{row['diffnb']}</diffnb>"""
    else:
        diffnb='<!-- optional diffnb not defined -->'
    
    if "diffprop" in row.keys():
        diffprop = f"""<diffprop>{row['diffprop']}</diffprop>"""
    else:
        diffprop='<!-- optional diffprop not defined -->'
    
    if "propa" in row.keys():
        propa = f"""<propa>{row['propa']}</propa>"""
    else:
        propa='<!-- optional propa not defined -->'
    
    if "propb" in row.keys():
        propb = f"""<propb>{row['propb']}</propb>"""
    else:
        propb='<!-- optional propb not defined -->'
    
    if "dens_cr" in row.keys():
        dens_cr = f"""<dens_cr>{row['dens_cr']}</dens_cr>"""
    else:
        dens_cr='<!-- optional dens_cr not defined -->'
    
    if "deepbasin" in row.keys():
        deepbasin = f"""<deepbasin>{row['deepbasin']}</deepbasin>"""
    else:
        deepbasin='<!-- optional deepbasin not defined -->'

#Hillslope diffuction optional / default parameters    
    if "sfail" in row.keys():
        sfail = f"""<sfail>{row['sfail']}</sfail>"""
    else:
        sfail='<!-- sfail default set -->'
    
    if "cfail" in row.keys():
        cfail = f"""<cfail>{row['cfail']}</cfail>"""
    else:
        cfail='<!-- cfail default set -->'   
    
    if "cslp" in row.keys():
        cslp = f"""<cslp>{row['cslp']}</cslp>"""
    else:
        cslp='<!-- cslp default set -->'     
        
    if "cmarine" in row.keys():
        cmarine = f"""<cmarine>{row['cmarine']}</cmarine>"""
    else:
        cmarine='<!-- cmarine default set -->'   
    
    if "criver" in row.keys():
        criver = f"""<criver>{row['criver']}</criver>"""
    else:
        criver='<!-- criver default set -->'   
        
 
#Flexure structure   
# Optional parameters are created only if they are present in the design table.     
    if "elasticGrid" in row.keys():
        elasticGrid = f"""<elasticGrid>data/{row['elasticGrid']}</elasticGrid>"""
    else:
        elasticGrid = '<!-- no optional elastic grid defined -->'
    
    if "fnx" in row.keys():
        fnx = f"""<fnx>{row['fnx']}</fnx>"""
    else:
        fnx = '<!-- no optional fnx defined -->'
    
    if "fny" in row.keys():
        fny = f"""<fny>{row['fny']}</fny>"""
    else:
        fny = '<!-- no optional fny defined -->' 
    
    if "poisson" in row.keys():
        poisson = f"""<poisson>{row['poisson']}</poisson>"""
    else:
        poisson = '<!-- no optional poisson defined -->'   
        
    if "elasticA1" in row.keys():
        elasticA1 = f"""<elasticA1>{row['elasticA1']}</elasticA1>"""
    else:
        elasticA1 = '<!-- no optional elasticA1 defined -->' 
    
    if "elasticA2" in row.keys():
        elasticA2 = f"""<elasticA2>{row['elasticA2']}</elasticA2>"""
    else:
        elasticA2 = '<!-- no optional elasticA2 defined -->'

# here we build a complete flexure structure on the assumption that if a dmantle is specified, you will have defined the rest of the required parameters for the xml, if not, errors will occur.         
    if "dmantle" in row.keys():
        flexure=f"""<flexure>
        <ftime>{row['ftime']}</ftime>
        
        {fnx}
        {fny}
        
        <dmantle>{row['dmantle']}</dmantle>
        <dsediment>{row['dsediment']}</dsediment>
        <youngMod>{row['youngMod']}</youngMod>
        <elasticH>{row['elasticH']}</elasticH>
       
        <boundary_W>{row['boundary_W']}</boundary_W>
        <boundary_E>{row['boundary_E']}</boundary_E>
        <boundary_S>{row['boundary_S']}</boundary_S>
        <boundary_N>{row['boundary_N']}</boundary_N>
        
        {poisson}
        {elasticGrid}
        {elasticA1}
        {elasticA2}
    </flexure>"""
    
    else:
        flexure= '<!-- flexure not enabled -->'     
    
##################################################################################################    
    if "stratdx" in row.keys():
        strata = f"""<strata>
        <!-- Stratal grid resolution [m] -->
        <stratdx>{row['stratdx']}.</stratdx>
        <!-- Stratal layer interval [a] -->
        <laytime>{row['laytime']}</laytime>
        <!-- Surface porosity  -->
        <poro0>{row['poro0']}</poro0>
        <!-- characteristic constant for Athy's porosity law  [/km] -->
        <poroC>{row['poroC']}</poroC>
        <!-- Only write the final sed.timeXX.hdf5 strat file-->
        <laststrat>True</laststrat>
    </strata>"""
    else:
        strata='<!-- stratal grid not enabled-->'

    if "fillmax" in row.keys():
        fillmax=f"""<fillmax>{row['fillmax']}.</fillmax>"""
    else:
        fillmax='<fillmax>100</fillmax>'
        

# Erodability
    if '.xml' in str(row['Erodibility']):
        try:
            erodibility_structure = erodability_dict[row['Erodibility']]
        except(KeyError): # this is the error returned if the specified xml file is not found.
            print('No Erodibility structure files found using single Erodibility values')
    else:
        erodability_structure = f"""<erodibility>{row['Erodibility']}</erodibility>"""         

# sedfluxfunction structure
    if "modeltype" in row.keys():
        sedfluxfunction=f"""<sedfluxfunction>
        <!-- Incision model type is defined with an integer between 0 and 4:
          + 0 - detachment limited (default) does not required to set additional parameters.
          + 1 - generalised undercapacity model (linear sedflux dependency) [cover effect]
          + 2 - parabolic sedflux dependency [tool & cover effect]
          + 3 - Turowski sedflux dependency [tool & cover effect]
          + 4 - saltation abrasion incision model
         See Hobley et al. (2011), JGR, 116 for more information.
    -->
    <modeltype>{row['modeltype']}</modeltype>
    <!-- Volumetric sediment transport capacity formulation is built with a stream power law
         and requires the definition of 2 exponents for water discharge (mt) and slope (nt). -->
    <mt>{row['mt']}</mt>
    <nt>{row['nt']}.</nt>
    <!-- Transportability of channel sediment (erodibility coefficient) -->
    <kt>{row['kt']}</kt>
    <!-- Power law relation between channel width and discharge -->
    <kw>{row['kw']}</kw>
    <b>{row['b']}</b>
    <!-- Erodibility dependence to the precipitation is defined with an exponent.
       Default value is set to 0. See Murphy et al. (2016), Nature, 532. -->
    <mp>{row['mp']}.</mp>
    <!-- Bedload versus slope dependency. This option changes the amount of incision based on
         the proportion of bedload material (i.e. gravels) present in stream. For any point in
         the landscape the amount of bedload material is assumed slope-dependent. The user can
         choose between the following options:
          + 0 - no dependency (default)
          + 1 - linear dependency
          + 2 - exponential growth
          + 3 - logarithmic growth
    -->
    <bedslp>{row['bedslp']}</bedslp>
</sedfluxfunction>"""
    else:
        sedfluxfunction='<!-- no optional sedfluxfunction defined -->'


# waveglobal structure
    if "wmodel" in row.keys():
        waveglobal=f"""<waveglobal>
    <!-- Wave model to consider either SWAN or WaveSed.
         Default is WaveSed (wmodel = 0). -->
    <wmodel>{row['wmodel']}</wmodel>
    <!-- Wave interval [a] -->
    <twave>{row['twave']}.</twave>
    <!-- Wave grid resolution [m] -->
    <wres>{row['wres']}.</wres>
    <!-- Maximum depth for wave influence [m] -->
    <wbase>{row['wbase']}</wbase>
    <!-- Number of wave climate temporal events. -->
    <events>1</events>
    <!-- Mean grain size diameter [m] -->
    <d50>{row['d50']}</d50>
    <!-- Wave sediment diffusion coefficient. Default is 50. -->
    <wCd>{row['wCd']}.</wCd>
    <!-- Wave sediment entrainment coefficient. Value needs to be
         set between ]0,1]. Default is 0.5 -->
    <wCe>{row['wCe']}</wCe>
    <!-- Maximum wave-induced erosion rate [m/yr] -->
    <wEro>{row['wEro']}</wEro>
    <!-- Maximum depth for wave influence [m] -->
    <wbase>{row['wbase']}</wbase>
    <!--  Steps used to perform sediment transport.
          Default is 1000. -->
    <tsteps>{row['tsteps']}</tsteps>
    <!--  Steps used to perform sediment diffusion.
          Default is 1000. -->
    <dsteps>{row['dsteps']}</dsteps>
</waveglobal>"""

    else:
        waveglobal='<!-- optional Wave global parameters structure not defined-->'

# Carbonate structure.
    if 'Carbonate' in row.keys():
        try:
            carbonate_structure = carbonate_dict[row['Carbonate']]
        except(KeyError): # this is the error returned if the specified xml file is not found.
            print('No carbonate structure files found')
    elif 'tcarb' in row.keys():
        carbonate_structure=f"""<carb>
    <!-- Specify initial basement structure (0) for hard rock and (1) for loose sediment. -->
    <baseMap>{row['baseMap']}</baseMap>
    <!-- Carbonate growth time interval [a] -->
    <tcarb>{row['tcarb']}.</tcarb>
    <!-- Specify the number of reef growth events -->
    <growth_events>{row['growth_events']}</growth_events>
    <!-- Specify Species 1 and 2 growth rates for specific reef growth events-->
    <event>
        <!-- Reef growth event start time [a] -->
        <gstart>{row['gstart']}.</gstart>
        <!-- Reef growth event end time [a] -->
        <gend>{row['gend']}.</gend>
        <!-- Species 1 growth rate during event [m/yr]. -->
        <growth_sp1>{row['growth_sp1']}</growth_sp1>
        <!-- Species 2 growth rate during event [m/yr]. -->
        <growth_sp2>{row['growth_sp2']}</growth_sp2>
    </event>
</carb>"""
    else:
        carbonate_structure='<!-- optional Carbonate growth parameters structure not defined-->'

       
#####################################################    
#Generate the main badlands experiment configurations by putting together all the elements from here on.
 
    badlands = f"""<?xml version="1.0" encoding="UTF-8"?>
    <badlands xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <!-- Regular grid structure -->
    <grid>
        <demfile>data/{row['Initial Topography']}</demfile>
        <boundary>fixed</boundary>
        <resfactor>{row['resfactor']}</resfactor>
    </grid>
    <!-- Simulation time structure -->
    <time>
        <!-- Simulation start time [a] -->
        <start>{row['start']}</start>
        <!-- Simulation end time [a] -->
        <end>{row['end']}</end>
        <mindt>{mindt}</mindt>
        <maxdt>{maxdt}</maxdt>
        <!-- Display interval [a] -->
        <display>{row['display']}</display>
    </time>
    <!-- Simulation stratigraphic structure -->
    {strata}
    <!--  Sea-level structure  -->
    <sea>
        <!--  Relative sea-level position [m]  -->
        <position>{row['position']}.</position>
        {curve}
        <!-- Limit flow network computation based on
        water depth [m] -->
        <limit>{row['limit']}</limit>
    </sea>
    
    {tectonic_structure}
    {precipitation_structure}

    <!--
    Stream power law parameters:
    -->
    <sp_law>
        <!-- Make the distinction between purely erosive models (0) and erosion /
        deposition ones (1). Default value is 1 -->
        <dep>{row['dep']}</dep>
        <!-- Critical slope used to force aerial deposition for alluvial plain,
        in [m/m] (optional). -->
        <slp_cr>{row['slp_cr']}</slp_cr>
        <!-- Maximum percentage of deposition at any given time interval from rivers
         sedimentary load in alluvial plain. Value ranges between [0,1] (optional). -->
        <perc_dep>{row['perc_dep']}</perc_dep>
        <!--
         Planchon & Darboux filling thickness limit [m]. This parameter is used
         to define maximum accumulation thickness in depression area per time
         step.
         -->
        {fillmax}
        <!--
         Values of m and n indicate how the incision rate scales
         with bed shear stress for constant value of sediment flux
         and sediment transport capacity.
         -->
        <m>{row['M']}</m>
        <n>{row['N']}</n>
        <!--
         The erodibility coefficient is scale-dependent and its value depend
         on lithology and mean precipitation rate, channel width, flood
         frequency, channel hydraulics.
        -->
        {erodability_structure}
        {diffnb}
        <!-- Proportion of marine sediment deposited on downstream nodes. It needs
         to be set between ]0,1[. Default value is 0.9 (optional). -->
        {diffprop}
        {propa}
        {propb}
        {dens_cr}
        {deepbasin}
    </sp_law>
    
    <creep>
        <!--  Surface diffusion coefficient [m2/a]  -->
        <caerial>{row['caerial']}</caerial>
        <!--  Marine diffusion coefficient [m2/a]  -->
        {cmarine}
        {cslp}
        <!-- River transported sediment diffusion
         coefficient in marine realm [m2/a] -->
        {criver}
        {sfail}
        {cfail}
    </creep>
    
    {sedfluxfunction}
    {flexure}
    {carbonate_structure}
    
    <!-- Output folder path -->
    <outfolder>Experiments/experiment_{row['Nexp']}</outfolder>
</badlands>
    """


    output_path = args[0]

    filename = os.path.join(output_path, "experiment_" + str(row['Nexp']) + ".xml")

    parser = ElementTree.XMLParser(target=CommentedTreeBuilder())

    #badlands_xml = etree.ElementTree(etree.fromstring(badlands))

    badlands_xml = etree.ElementTree(etree.fromstring(badlands, parser))

    with open(filename, 'wb') as f:
        badlands_xml.write(f, encoding='utf-8', xml_declaration=True)
        
"""
DoEgen to badlands converter
Takes a pandas dataframe of the DoEgen .csv output and writes a badlands .xml for every row.
"""

def doegen_to_badlands(input_doegen_df, output_path = None):
    if output_path is None:
        output_path = Path('badlands_test_xml/')

    os.makedirs(Path(output_path), exist_ok = True)

    args = [output_path]
    
    input_doegen_df.apply(badlands_encode, axis=1, args = args)


def main(input_doegen_csv):
#    input_doegen_csv = args.input / args.file
    input_doegen_df = pd.read_csv(input_doegen_csv)
    # set up calculated / dependent parameters
    M=round((input_doegen_df['N']*input_doegen_df['MNrat']),4)
    input_doegen_df['M']=M
    doegen_to_badlands(input_doegen_df)


#if __name__ == "__main__":
#    main()
