# -*- coding: utf-8 -*-
"""This module contains functions to load LCA JSON-data from dataclass into En15804LcaData objects"""



from teaser.logic.buildingobjects.buildingphysics.en15804indicatorvalue import En15804IndicatorValue


def load_en15804_lca_data_id(lca_data, lca_id, data_class):
    """LCA-data loader with id as identification.

    Loads LCA-data specified in the JSON by given LCA-ID.

    Parameters
    ----------
    lca_data : En15804LcaData()
        instance of TEASERS En15804LcaData class to store selected LCA-dataset

    lca_id : str
        uuid of LCA-data from JSON

    data_class : DataClass()
        DataClass containing the bindings for LCA-data
        (typically this is the data class stored in prj.data,
        but the user can individually change that.

    """
    
    
    binding = data_class.lca_data_bind

    for id, data in binding.items():
        if id != "version":
            if id == lca_id:
                
                lca_data.lca_data_id = id
                lca_data.name = data["name"]

                ref_flow_value = data["ref_flow"]["value"]
                lca_data.ref_flow_unit = data["ref_flow"]["unit"]
                
                pere = En15804IndicatorValue()
                pert = En15804IndicatorValue()
                penre = En15804IndicatorValue()
                penrm = En15804IndicatorValue()
                penrt = En15804IndicatorValue()
                sm = En15804IndicatorValue()
                rsf = En15804IndicatorValue()
                nrsf = En15804IndicatorValue()
                fw = En15804IndicatorValue()
                hwd = En15804IndicatorValue()
                nhwd = En15804IndicatorValue()
                rwd = En15804IndicatorValue()
                cru = En15804IndicatorValue()
                mfr = En15804IndicatorValue()
                mer = En15804IndicatorValue()
                eee = En15804IndicatorValue()
                eet = En15804IndicatorValue()
                gwp = En15804IndicatorValue()
                odp = En15804IndicatorValue()
                pocp = En15804IndicatorValue()
                ap = En15804IndicatorValue()
                ep = En15804IndicatorValue()
                adpe = En15804IndicatorValue()
                adpf = En15804IndicatorValue()
                
                pere.set_values(**data["pere"])
                pert.set_values(**data["pert"])
                penre.set_values(**data["penre"])
                penrm.set_values(**data["penrm"])
                penrt.set_values(**data["penrt"])
                sm.set_values(**data["sm"])
                rsf.set_values(**data["rsf"])
                nrsf.set_values(**data["nrsf"])
                fw.set_values(**data["fw"])
                hwd.set_values(**data["hwd"])
                nhwd.set_values(**data["nhwd"])
                rwd.set_values(**data["rwd"])
                cru.set_values(**data["cru"])
                mfr.set_values(**data["mfr"])
                mer.set_values(**data["mer"])
                eee.set_values(**data["eee"])
                eet.set_values(**data["eet"])
                gwp.set_values(**data["gwp"])
                odp.set_values(**data["odp"])
                pocp.set_values(**data["pocp"])
                ap.set_values(**data["ap"])
                ep.set_values(**data["ep"])
                adpe.set_values(**data["adpe"])
                adpf.set_values(**data["adpf"]) 
                
             
                lca_data.pere = pere * (1/ref_flow_value)
                lca_data.pert = pert * (1/ref_flow_value)
                lca_data.penre = penre * (1/ref_flow_value)
                lca_data.penrm = penrm * (1/ref_flow_value)
                lca_data.penrt = penrt * (1/ref_flow_value)
                lca_data.sm = sm * (1/ref_flow_value)
                lca_data.rsf = rsf * (1/ref_flow_value)
                lca_data.nrsf = nrsf * (1/ref_flow_value)
                lca_data.fw = fw * (1/ref_flow_value)
                lca_data.hwd = hwd * (1/ref_flow_value)
                lca_data.nhwd = nhwd * (1/ref_flow_value)
                lca_data.rwd = rwd * (1/ref_flow_value)
                lca_data.cru = cru * (1/ref_flow_value)
                lca_data.mfr = mfr * (1/ref_flow_value)
                lca_data.mer = mer * (1/ref_flow_value)
                lca_data.eee = eee * (1/ref_flow_value)
                lca_data.eet = eet * (1/ref_flow_value)
                lca_data.gwp = gwp * (1/ref_flow_value)
                lca_data.odp = odp * (1/ref_flow_value)
                lca_data.pocp = pocp * (1/ref_flow_value)
                lca_data.ap = ap * (1/ref_flow_value)
                lca_data.ep = ep * (1/ref_flow_value)
                lca_data.adpe = adpe * (1/ref_flow_value)
                lca_data.adpf = adpf * (1/ref_flow_value)
                
                lca_data.ref_flow_value = 1
                
                
                
                
                if data["fallback"]:
                    lca_data.load_fallbacks(data["fallback"], data_class)
                    
                    lca_data.add_fallbacks()
                else:
                    lca_data.fallback = []
                    

def load_en15804_lca_data_fallback_id(lca_data, lca_id, data_class):
    """LCA-data-fallback loader with id as identification.

    Loads LCA-data-fallbacks specified in the JSON by given LCA-ID.
    LCA-fallbacks are specified in an seperated JSON-file to clarify they are
    just partial defined.

    Parameters
    ----------
    lca_data : En15804LcaData()
        instance of TEASERS En15804nLcaData class to store selected LCA-fallback

    lca_id : str
        id of LCA-data from JSON

    data_class : DataClass()
        DataClass containing the bindings for LCA-data-fallbacks
        (typically this is the data class stored in prj.data,
        but the user can individually change that.

    """
    

    binding = data_class.lca_data_fallback_bind
    
    for id, data in binding.items():
        if id != "version":
            if id == lca_id:
                
                lca_data.lca_data_id = id
                lca_data.name = data["name"]

                lca_data.ref_flow_value = data["ref_flow"]["value"]
                lca_data.ref_flow_unit = data["ref_flow"]["unit"]
                
                pere = En15804IndicatorValue()
                pert = En15804IndicatorValue()
                penre = En15804IndicatorValue()
                penrm = En15804IndicatorValue()
                penrt = En15804IndicatorValue()
                sm = En15804IndicatorValue()
                rsf = En15804IndicatorValue()
                nrsf = En15804IndicatorValue()
                fw = En15804IndicatorValue()
                hwd = En15804IndicatorValue()
                nhwd = En15804IndicatorValue()
                rwd = En15804IndicatorValue()
                cru = En15804IndicatorValue()
                mfr = En15804IndicatorValue()
                mer = En15804IndicatorValue()
                eee = En15804IndicatorValue()
                eet = En15804IndicatorValue()
                gwp = En15804IndicatorValue()
                odp = En15804IndicatorValue()
                pocp = En15804IndicatorValue()
                ap = En15804IndicatorValue()
                ep = En15804IndicatorValue()
                adpe = En15804IndicatorValue()
                adpf = En15804IndicatorValue()
                
                pere.set_values(**data["pere"])
                pert.set_values(**data["pert"])
                penre.set_values(**data["penre"])
                penrm.set_values(**data["penrm"])
                penrt.set_values(**data["penrt"])
                sm.set_values(**data["sm"])
                rsf.set_values(**data["rsf"])
                nrsf.set_values(**data["nrsf"])
                fw.set_values(**data["fw"])
                hwd.set_values(**data["hwd"])
                nhwd.set_values(**data["nhwd"])
                rwd.set_values(**data["rwd"])
                cru.set_values(**data["cru"])
                mfr.set_values(**data["mfr"])
                mer.set_values(**data["mer"])
                eee.set_values(**data["eee"])
                eet.set_values(**data["eet"])
                gwp.set_values(**data["gwp"])
                odp.set_values(**data["odp"])
                pocp.set_values(**data["pocp"])
                ap.set_values(**data["ap"])
                ep.set_values(**data["ep"])
                adpe.set_values(**data["adpe"])
                adpf.set_values(**data["adpf"]) 
                
                lca_data.pere = pere 
                lca_data.pert = pert
                lca_data.penre = penre
                lca_data.penrm = penrm
                lca_data.penrt = penrt
                lca_data.sm = sm
                lca_data.rsf = rsf
                lca_data.nrsf = nrsf
                lca_data.fw = fw
                lca_data.hwd = hwd
                lca_data.nhwd = nhwd
                lca_data.rwd = rwd
                lca_data.cru = cru
                lca_data.mfr = mfr
                lca_data.mer = mer
                lca_data.eee = eee
                lca_data.eet = eet
                lca_data.gwp = gwp
                lca_data.odp = odp
                lca_data.pocp = pocp
                lca_data.ap = ap
                lca_data.ep = ep
                lca_data.adpe = adpe
                lca_data.adpf = adpf
                
                lca_data.fallback = None