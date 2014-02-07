import pandas as pd
import numpy as np
import xlrd
import re


##########################################################################
def sample_dict(dataframe, station_names_col):
    """Generates dictionary of stations, with an array per each station containing relevant sample
    numbers.
    INPUTS:
        dataframe: a master Pandas dataframe, with index corresponding to sample #
        station_names_col: the column of dataframe (ex. dataframe['Station Names']) with station names
    OUTPUT:
        dictionary of stations, with corresponding sample numbers"""
    st_indices=dict((name,index) for index, name in enumerate(station_names_col))
    end_indices=[st_indices[item] for item in pd.unique(station_names_col)]
    start=0
    s_dict=dict()
    for i in range(0, np.size(end_indices)):
        st_sample_nums=[k for k in dataframe.index[start:end_indices[i]+1]]
        s_dict['%s' % pd.unique(station_names_col)[i]]=st_sample_nums
        start=end_indices[i]+1
    return s_dict


##########################################################################
def gsdatain_ange2(xls_name, sheet_name=None, load_duplicates=False):
    """Load g.s. data from excel file from Ange's sediment lab
    INPUTS:
        xls_name = 'excel file name with extension'
        sheet_name = 'sheet name containing grain size data'  (usually SDSZ or Grain Size)
        load_duplicates = False (default). Set True if want to load duplicates
    OUTPUTS:
        master: all data from grain size and sample info tabs, indexed by sample #
        s_dict: dictionary of names of stations, and the sample #s associated with those stations

    Look up data for any particular station or sample like so if you want to look up a station:
        master.loc[s_dict['NAME OF STATION IN DICTIONARY']]
        or if you know the sample # (ex. sample 15), just:
        master.loc[15] """

    #"""CAN MANUALLY SET FILENAMES AND OTHER SETTINGS HERE"""
    filename = 'L712AK.xls'
    sheet = 'SDSZ'
    #Set to True if you want to load all stations, including duplicates.
    if load_duplicates is None:
        load_duplicates = False
    #"""***************************************************"""


    '''Load file'''

    if xls_name is None:
        xls_name = filename
        if sheet == 'SDSZ' or sheet == 'Grain Size':
            try:
                data=pd.read_excel(xls_name, sheet, header=2, index_col=0)
            except xlrd.XLRDError:
                if sheet == 'SDSZ':
                    try:
                        data=pd.read_excel(xls_name, 'Grain Size', header=2, index_col=0)
                    except xlrd.XLRDError:
                        print('Manually input sheet name, try again')
                elif sheet == 'Grain Size':
                    try:
                        data=pd.read_excel(xls_name, 'SDSZ', header=2, index_col=0)
                    except xlrd.XLRDError:
                        print('Manually input sheet name, try again')
        else:
            data=pd.read_excel(xls_name, sheet, header=2, index_col=0)
    else:
        try:
            data=pd.read_excel(xls_name, sheet_name, header=2, index_col=0)
        except xlrd.XLRDError:
            print('Failed to read excel spreadsheet')


    '''Load sample locality data'''
    l_data=pd.read_excel(xls_name, 'Sample Info', header=2, index_col=None)

    '''Load Sample Names'''

    names=data.index[pd.notnull(data.index)]    # remove blank rows

    try:    # remove weird table at bottom
        names=names[0:np.where(names == 'Lab Fraction Names')[0][0]]
    except IndexError:
        print("could not find 'Lab Fractions' table at bottom of sheet")
        pass

    if load_duplicates is False:    # removes samples listed as "duplicates" from the lab
        names=[name for name in names if 'duplicate' not in name]

    '''Get phi and mm bins'''   # MAY BE UNNECESSARY?
    #phibins = data.columns[1:np.where(data.columns == r" % Gravel")[0][0]]
    #mmbins = [2**(-phi) for phi in phibins]

    data=data.T[names]
    data=data.T

    '''Make master table, combining grain size and station data, index by sample number'''
    master=l_data.set_index(l_data['Sample #']).join(data.set_index(l_data['Sample #']))
    master['Station Name']=l_data.index     # add Station Name to master data

    '''Make lists of stations and depths'''
    s_dict=sample_dict(master, master['Station Name'])

    '''Make array of station names in the master dataframe'''
    s_names=pd.unique(master['Station Name'])



    return master, s_dict, s_names


##########################################################################
def split_layers(master, s_dict, s_name):
    '''Read dataframe generated from gsdatain_****.py for a given locality and for a particular station, split into plottable layers
    INPUT:
        master: master dataframe
        s_dict: dictionary of stations and corresponding sample numbers
        s_name: station name in the master dataframe that you want data for
    OUTPUT:
        layer_data, list of dataframes for each layer, each dataframe contains info from master
    USAGE:
        If you have an array of unique station names (s_names):
            for c in s_names:
                layer=split_layers(master,s_dict,c)
                ...
        And remember, output (layer_data) is a list, containing dataframes, not a dataframe itself
    '''
    #DETERMINE WHERE LAYER SPLITS ARE
    st_data = master.loc[s_dict[s_name]]

    try:
        top_bot_1=[[float(re.split('-',i)[0]),float(re.split('-',i)[1])] for i in st_data['Interval']]
    except ValueError:
        top_bot_1 = np.array([[0,3]])

    layer_splits=[]
    if np.shape(top_bot_1)[0] == 1:
        for i in range(0,1):
            layer_splits.append(1)
    else:
        for i in range(0, np.shape(top_bot_1)[0]-1):
            if np.abs(top_bot_1[i+1][0] - top_bot_1[i][1]) >= 0.5:
                layer_splits.append(st_data.index[i])   # split after each of these sample numbers

    if layer_splits[0] == 1 and np.size(layer_splits) == 1:
        num_layers = 1
    else:
        num_layers = np.size(layer_splits)+1

    layer_data=[]
    for i in range(0, num_layers):
        if i == 0 and num_layers ==1:
            layer_data.append(st_data)
        elif i == 0:
            layer_data.append(st_data.loc[:layer_splits[i]])
        elif i == np.max(num_layers)-1:
            layer_data.append(st_data.loc[layer_splits[i-1]+1:])
        else:
            layer_data.append(st_data.loc[layer_splits[i-1]+1:layer_splits[i]])

    return layer_data


###########################################################################
def get_modes(layer_data):
    '''Returns a list of the modes (in phi) for a layer from split_layers
    Inputs:
        layer_data -> a pandas DataFrame for a single sample
    Output:
        mode_list -> list of modes for the sample in phi
    '''
            modes=[i for i in layer_data[' Mode in Phi']]
            mode_list=[]
            for i in modes:
                if type(i) == float:
                    mode_list.append([i])
                elif i == ' None':
                    pass
                else:
                    mode_list.append([float(j) for j in re.split('&', i)])
            return mode_list

###########################################################################
def get_graph_D(layer_data,percentile):
    '''For a given layer, returns the graphically interpolated grain size in phi
    for a given percentile (ex. d50, d84).
    Percentiles usually available are:
        [d5,d10,d16,d25,d50,d75,d84,d90,d95]
    Inputs:
        layer_data -> a pandas DataFrame for a single sample
        percentile -> integer for the percentile desired (ex. 50 for d50)
    Output:
        size -> grain size in 50 representing the 'd' percentile of the sample
    '''
    try:
        percentile_string = ' D%i Phi' % percentile   # finds the column with the appropriate percentile
        d = layer_data[percentile_string]  # an object containing sample# and phi
        size = float(d)   # just the grain size in ph
    except KeyError:
        print('No percentile calculated for d%i.\nUse an integer from:\n[d5,d10,d16,d25,d50,d75,d84,d90,d95]' %percentile)
        size = None

    return size
