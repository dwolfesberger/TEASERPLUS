# -*- coding: utf-8 -*-

"""This module contains functions to save lca_data and lca_data_fallback objects."""
import warnings
import json
import teaser.logic.utilities as utilities
import collections


def save_lca_data(lca_data, data_class):
    """lca data saver.

    Appends given LCA-dataset and its properties to the LCA-data-JSON file.

    Parameters
    ----------
    lca_data : En15804LcaData()
        instance of TEASERS En15804LcaData class to be saved

    data_class : DataClass()
        DataClass containing the bindings for LCA-Data, TypeBuildingElement and
        Material (typically this is the data class stored in prj.data,
        but the user can individually change that.

    """
    data_class.lca_data_bind["version"] = "0.7"
    add_to_json = True

    warning_text = ("LCA data with same name and id already "
                    "exists in JSON")

    for id, check in data_class.lca_data_bind.items():
        if id != "version":
            if check["name"] == lca_data.name and check == lca_data.lca_data_id:

                warnings.warn(warning_text)
                
                add_to_json = False
                break

    if add_to_json is True:
        data_class.lca_data_bind[lca_data.lca_data_id] = collections.OrderedDict()
        data_class.lca_data_bind[lca_data.lca_data_id]["name"] = lca_data.name
        
        data_class.lca_data_bind[lca_data.lca_data_id]["pere"] = lca_data.pere.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["pert"] = lca_data.pert.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["penre"] = lca_data.penre.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["penrm"] = lca_data.penrm.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["penrt"] = lca_data.penrt.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["sm"] = lca_data.sm.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["rsf"] = lca_data.rsf.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["nrsf"] = lca_data.nrsf.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["fw"] = lca_data.fw.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["hwd"] = lca_data.hwd.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["nhwd"] = lca_data.nhwd.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["rwd"] = lca_data.rwd.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["cru"] = lca_data.cru.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["mfr"] = lca_data.mfr.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["mer"] = lca_data.mer.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["eee"] = lca_data.eee.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["eet"] = lca_data.eet.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["gwp"] = lca_data.gwp.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["odp"] = lca_data.odp.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["pocp"] = lca_data.pocp.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["ap"] = lca_data.ap.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["ep"] = lca_data.ep.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["adpe"] = lca_data.adpe.get_values_as_dict()
        data_class.lca_data_bind[lca_data.lca_data_id]["adpf"] = lca_data.adpf.get_values_as_dict()
        
        ref_flow_dict = {}
        ref_flow_dict["value"] = lca_data.ref_flow_value
        ref_flow_dict["unit"] = lca_data.ref_flow_unit
        
        data_class.lca_data_bind[lca_data.lca_data_id]["ref_flow"] = ref_flow_dict
        data_class.lca_data_bind[lca_data.lca_data_id]["fallback"] = lca_data.fallback
        
        

    with open(utilities.get_full_path(data_class.path_lcad), 'w') as file:
        

        
        file.write("test")
        
        json.dump(
            data_class.lca_data_bind,
            file,
            indent=4,
            separators=(',', ': '))
        

def save_lca_data_fallback(lca_data, data_class):
    """lca data fallback saver.

    Appends given LCA-dataset and its properties to the LCA-data-fallback-JSON file.
    
    Parameters
    ----------
    lca_data : En15804LcaData()
        instance of TEASERS En15804LcaData class to be saved

    data_class : DataClass()
        DataClass containing the bindings for LCA-Data, TypeBuildingElement and
        Material (typically this is the data class stored in prj.data,
        but the user can individually change that.

    """
    data_class.lca_data_fallback_bind["version"] = "0.7"
    add_to_json = True

    warning_text = ("LCA data with same name and id already "
                    "exists in JSON")

    for id, check in data_class.lca_data_fallback_bind.items():
        if id != "version":
            if check["name"] == lca_data.name and check == lca_data.lca_data_id:

                warnings.warn(warning_text)
                
                add_to_json = False
                break

    if add_to_json is True:
        data_class.lca_data_fallback_bind[lca_data.lca_data_id] = collections.OrderedDict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["name"] = lca_data.name
        
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["pere"] = lca_data.pere.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["pert"] = lca_data.pert.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["penre"] = lca_data.penre.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["penrm"] = lca_data.penrm.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["penrt"] = lca_data.penrt.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["sm"] = lca_data.sm.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["rsf"] = lca_data.rsf.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["nrsf"] = lca_data.nrsf.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["fw"] = lca_data.fw.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["hwd"] = lca_data.hwd.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["nhwd"] = lca_data.nhwd.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["rwd"] = lca_data.rwd.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["cru"] = lca_data.cru.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["mfr"] = lca_data.mfr.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["mer"] = lca_data.mer.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["eee"] = lca_data.eee.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["eet"] = lca_data.eet.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["gwp"] = lca_data.gwp.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["odp"] = lca_data.odp.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["pocp"] = lca_data.pocp.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["ap"] = lca_data.ap.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["ep"] = lca_data.ep.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["adpe"] = lca_data.adpe.get_values_as_dict()
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["adpf"] = lca_data.adpf.get_values_as_dict()
        
        ref_flow_dict = {}
        ref_flow_dict["value"] = lca_data.ref_flow_value
        ref_flow_dict["unit"] = lca_data.ref_flow_unit
        
        data_class.lca_data_fallback_bind[lca_data.lca_data_id]["ref_flow"] = ref_flow_dict
                
        

    with open(utilities.get_full_path(data_class.path_lcad_fallback), 'w') as file:
        file.write(json.dumps(
            data_class.lca_data_fallback_bind,
            indent=4,
            separators=(',', ': ')))
