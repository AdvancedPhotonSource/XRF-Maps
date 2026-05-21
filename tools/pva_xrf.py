#!/usr/bin/env python

import sys
import json
import time
import random
import numpy as np
import pvaccess as pva
import redis 
import pyxrfmaps as px

element_csv_filename = "../reference/xrf_library.csv"
element_henke_filename = "../reference/henke.xdr"

class XRF_Stream_Source():
    def __init__(self, name, param_override, CHANNEL):
        self.node_name = name
        self.bnode_name = bytes(name, 'utf-8')

        self.dataset_directory = "/"
        self.dataset_name = "Stream"

        self.row = 0
        self.col = 0
        self.width = 50
        self.height = 50

        # Select fitting routine
        #self.fit_rout = px.fitting.routines.nnls()
        self.fit_rout = px.fitting.routines.roi()

        # Use Gausian Model
        self.model = px.fitting.models.GaussModel()

        # Load fit parameters 
        self.po = param_override

        self.channel = pva.Channel(CHANNEL)
        # monitor
        self.channel.monitor(self.monitor)

        self.sprectra_streamer = px.workflow.SpectraNetStreamer("43434")
        self.sprectra_streamer.set_send_counts(True)
        self.sprectra_streamer.set_send_spectra(False)

    def init_fit_routine(self, int_spec):
        # Initialize model and fit routine with fit parameters
        self.energy_range = px.get_energy_range(int_spec.size, po.fit_params)
        self.model.update_fit_params_values(po.fit_params)
        self.fit_rout.initialize(model, po.elements_to_fit, energy_range)

    def config_change_handler(self, msg):
        if msg['data'] == self.bnode_name:
            print(msg)

    def fit_spec(self, spec, id):
        if spec is not None:
            # create a stream block
            sb = px.StreamBlock(-1, self.row, self.col, self.height, self.width, self.dataset_directory, self.dataset_name)
            sb.init_fitting_blocks( {px.FittingRoutines.ROI: self.fit_rout}, self.po.elements_to_fit)
            #sb.spectra = spec
            sb.fitting_blocks[px.FittingRoutines.ROI] = px.StreamFittingBlock()
            # fit the spectra
            sb.fitting_blocks[px.FittingRoutines.ROI].fit_counts = self.fit_rout.fit_counts(self.model, spec, self.po.elements_to_fit)
            print (self.col, id, sb.fitting_blocks[px.FittingRoutines.ROI].fit_counts)
            self.sprectra_streamer.stream(sb)
            self.col += 1
            if self.col >= self.width:
                self.col = 0
                self.row += 1
                if self.row >= self.height:
                    self.row = 0
            if self.row >= self.height:
                    self.row = 0
            # push to redis
            ##

    def monitor(self, pv):
        xdim = 0
        yxim = 0
        data_arr = None
        print(pv['attribute'][3])
        if len(pv['dimension']) == 2:
            xdim = pv['dimension'][0]['size']
            ydim = pv['dimension'][1]['size']
            if 'shortValue' in pv['value'][0]:
                data_arr = pv['value'][0]['shortValue'].reshape(xdim * ydim, 1).astype(np.float32) # TODO: redo so we get proper dims
            elif 'ushortValue' in pv['value'][0]:
                data_arr = pv['value'][0]['ushortValue'].reshape(xdim * ydim, 1).astype(np.float32) # TODO: redo so we get proper dims
            elif 'intValue' in pv['value'][0]:
                data_arr = pv['value'][0]['intValue'].reshape(xdim * ydim, 1).astype(np.float32) # TODO: redo so we get proper dims
            elif 'uintValue' in pv['value'][0]:
                data_arr = pv['value'][0]['uintValue'].reshape(xdim * ydim, 1).astype(np.float32) # TODO: redo so we get proper dims
            elif 'floatValue' in pv['value'][0]:
                data_arr = pv['value'][0]['floatValue'].reshape(xdim * ydim, 1) # TODO: redo so we get proper dims
            else:
                print(pv['value'][0].keys())
            #print('Got image: %d' % pv['uniqueId'], xdim, ydim, data_arr[:] )
            self.fit_spec(data_arr, pv['uniqueId'])
        #'value', 'codec', 'compressedSize', 'uncompressedSize', 'dimension', 'uniqueId', 'dataTimeStamp', 'attribute', 'descriptor', 'alarm', 'timeStamp', 'display'

def main():
    CHANNEL = None
    config_dict = None
    param_override = None
    p = None
    if len(sys.argv) < 2:
        print("Please provide config file as second arg")
        exit(1)
    
    print ('Loading json file ', sys.argv[1])
    with open(sys.argv[1]) as json_file:
        config_dict = json.load(json_file)
        print (config_dict)

    if config_dict == None:
        print('Error loading json file ', sys.argv[1])

    if 'redis_host' in config_dict:
        r = redis.Redis(config_dict['redis_host'])
        # get initial config
        doc = r.json().get('xrf_workers', '$')
        if not config_dict['node_name'] in doc[0]:
            print ('Config not found for ', config_dict['node_name'])
            print (doc)
            exit(1)
            
        # subscribe to config_change for live changes
        #p = r.pubsub()
        #p.subscribe(**{'config_change': config_change_handler})
        # connect to PVA channel
        CHANNEL = str(doc[0][config_dict['node_name']]['PVA']) # 'bdpSimDetector:Pva1:Image'
        print(config_dict['node_name'], CHANNEL)

        # get param override from redis

    if CHANNEL == None:
        if 'PVA' in config_dict:
            CHANNEL = config_dict['PVA']
        else:
            print("Please add 'PVA' to config.json or add a 'redis_host' to get config data from")
            exit(1)

    if param_override == None:
        if "xrf_param_override" in config_dict:
            config_dict["xrf_param_override"]["dataset_dir"] = config_dict["xrf_param_override"]["dataset_dir"] + '/'
            param_override = px.load_override_params(config_dict["xrf_param_override"]["dataset_dir"], config_dict["xrf_param_override"]["detector_num"], True)

    # load ref info
    px.load_element_info(element_henke_filename, element_csv_filename)

    if 'node_name' not in config_dict:
        config_dict['node_name'] = 'node_' + str(random.randint(100,1000))
        print('Could not find "node_name" in config.json. Using random generated name = ', config_dict['node_name'] )
    source = XRF_Stream_Source(config_dict['node_name'], param_override, CHANNEL)

    # loop and check config changes 
    for i in range(1000):
        #if p is not None:
            #p.get_message()
        time.sleep(1)
    #time.sleep(1000)


if __name__ == '__main__':
    main()


