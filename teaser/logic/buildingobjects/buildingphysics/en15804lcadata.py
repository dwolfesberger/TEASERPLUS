# -*- coding: utf-8 -"*-
"""This module includes the En15804LcaData-class"""

import uuid

from teaser.logic.buildingobjects.buildingphysics.en15804indicatorvalue import En15804IndicatorValue

import teaser.data.input.lca_data_input as lca_data_input
import teaser.data.output.lca_data_output as lca_data_output

class En15804LcaData(object):
    """En15804LcaData class
    
    This class can be used to store environmental product declarations
    according to EN 15804. All 25 environmental indicators of EN 15804 
    can be specified.
    
    Parameters
    ----------
    
    parent : object, optional
        parent of the object. Default is None
    
    Attributes
    ----------
        lca_data_id : str
            universally unique identifier of the LCA-dataset
        name : str
            name of the EPD (e.g. name of the material)
        ref_flow_value : float [-]
            value of the reference flow
        ref_flow_unit : str
            unit of the reference flow
        pere : En15804IndicatorValue [MJ]
            Use of renewable primary energy
        perm : En15804IndicatorValue [MJ]
            Use of renewable primary energy resources used as raw materials
        pert : En15804IndicatorValue [MJ]
            Total use of renewable primary energy resources
        penre : En15804IndicatorValue [MJ]
            Use of non renewable primary energy 
        penrm : En15804IndicatorValue [MJ]
            Use of non renewable primary energy resources used as raw materials
        penrt : En15804IndicatorValue [MJ]
            Total use of non renewable primary energy resource 
        sm : En15804IndicatorValue [kg]
            Use of secondary material 
        rsf : En15804IndicatorValue [MJ]
            Use of renewable secondary fuels 
        nrsf : En15804IndicatorValue [MJ]
            Use of non renewable secondary fuels 
        fw : En15804IndicatorValue [m^3]
            Use of net fresh water 
        hwd : En15804IndicatorValue [kg]
            Hazardous waste disposed
        nhwd : En15804IndicatorValue [kg]
            Non hazardous waste dispose 
        rwd : En15804IndicatorValue [kg]
            Radioactive waste disposed
        cru : En15804IndicatorValue [kg]
            Components for re-use 
        mfr : En15804IndicatorValue [kg]
            Materials for recycling 
        mer : En15804IndicatorValue [kg]
            Materials for energy recovery
        eee : En15804IndicatorValue [MJ]
            Exported electrical energy 
        eet : En15804IndicatorValue [MJ]
            Exported thermal energy
        gwp : En15804IndicatorValue [kg CO2 eq.]
            Global warming potential
        odp : En15804IndicatorValue [kg R11 eq.]
            Depletion potential of the stratospheric ozone layer 
        pocp : En15804IndicatorValue [kg Ethene eq.]
            Formation potential of tropospheric ozone 
        ap : En15804IndicatorValue [kg SO2 eq.]
            Acidification potential of soil and water 
        ep : En15804IndicatorValue [kg Phosphate eq.]
            Eutrophication potential 
        adpe : En15804IndicatorValue [kg Sb eq.]
            Abiotic depletion potential for non fossil resources
        adpf : En15804IndicatorValue [MJ]
            Abiotic depletion potential for fossil resources 
        fallback : dict
            Dictonarie with stages as key and LCA-datasets as values. 
            The LCA dataset serves as the fallback for the specified stage
        _fallback_added : boolean
            is True when fallbacks allready added to "lca_data"
    """
    def __init__(self, parent = None):
      
        self.lca_data_id = str(uuid.uuid1())
        
        self._name = None
        self.parent = parent
        
        self._ref_flow_value = None
        self._ref_flow_unit = None
        
        self._pere = En15804IndicatorValue()
        self._perm = En15804IndicatorValue()
        self._pert = En15804IndicatorValue()
        self._penre = En15804IndicatorValue()
        self._penrm = En15804IndicatorValue()
        self._penrt = En15804IndicatorValue()
        self._sm = En15804IndicatorValue()
        self._rsf = En15804IndicatorValue()
        self._nrsf = En15804IndicatorValue()
        self._fw = En15804IndicatorValue()
        self._hwd = En15804IndicatorValue()
        self._nhwd = En15804IndicatorValue()
        self._rwd = En15804IndicatorValue()
        self._cru = En15804IndicatorValue()
        self._mfr = En15804IndicatorValue()
        self._mer = En15804IndicatorValue()
        self._eee = En15804IndicatorValue()
        self._eet = En15804IndicatorValue()
        self._gwp = En15804IndicatorValue()
        self._odp = En15804IndicatorValue()
        self._pocp = En15804IndicatorValue()
        self._ap = En15804IndicatorValue()
        self._ep = En15804IndicatorValue()
        self._adpe = En15804IndicatorValue()
        self._adpf = En15804IndicatorValue()
        
        self._fallback = None
        self._fallback_added = False
        
        
    def _check_unit(self, unit, unit_expected, var_name = None):
        """function to check if unit equals the excpected unit. The  
        name of the checked variable can be passed for the error-message
        
        Parameters
        ----------
        unit : str
            Unit to be checked
        unit_expected : str
            unit to be expected
        var_name : str, optional
            Name of the Variable used for the error-message. The default is 
            None.

        Returns
        -------
        Boolean

        """
        
        if unit == unit_expected:
            return(True)
        else:
            if var_name:
                print("Variable '{}' should be specified in {}".format(var_name, unit_expected))
            else:
                print("Please insert value in {}".format(unit_expected))
            return(False)
        
    def _check_en15804indicatorvalue_class(self, value, var_name = None):
        """fuction to check if value is En15804IndicatorValue-Object. The  
        namf of the checked variable can be passed for the error-message
        

        Parameters
        ----------
        value : En15804IndicatorValue
            value to be checked
        var_name : str, optional
            Name of the Variable used for the error-message. The default is 
            None.

        Returns
        -------
        Boolean

        """

        if value:
            if isinstance(value, En15804IndicatorValue): 
                return(True)
            else:
                if var_name:
                    print("'{}' has to be an En15804IndicatorValue-Object!".format(var_name))
                else:
                    print("Please insert a En15804IndicatorValue-Object!")
                return(False)
        else:
            return(False)
        
        
        
        

    @property
    def ref_flow_value(self):
        return self._ref_flow_value
    
    @ref_flow_value.setter
    def ref_flow_value(self, value):
        self._ref_flow_value = value
        
    @property
    def ref_flow_unit(self):
        return self._ref_flow_unit
    
    @ref_flow_unit.setter
    def ref_flow_unit(self, value):
        self._ref_flow_unit = value
        
    @property
    def fallback(self):
        return self._fallback 
    
    @fallback.setter
    def fallback(self, value):
        self._fallback = value
        
    @property
    def pere(self):
        return self._pere
    
    @pere.setter
    def pere(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "pere"):
            if self._check_unit(value.unit, "MJ", "pere"):
                self._pere = value
                
    @property
    def perm(self):
        return self._perm
    
    @perm.setter
    def perm(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "perm"):
            if self._check_unit(value.unit, "MJ", "perm"):
                self._perm = value
                         
    @property
    def pert(self):
        return self._pert 
    
    @pert.setter
    def pert(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "pert"):
            if self._check_unit(value.unit, "MJ", "pert"):
                self._pert = value
                
    @property
    def penre(self):
        return self._penre
    
    @penre.setter
    def penre(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "penre"):
            if self._check_unit(value.unit, "MJ", "penre"):
                self._penre = value
                
    @property
    def penrm(self):
        return self._penrm
    
    @penrm.setter
    def penrm(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "penrm"):
            if self._check_unit(value.unit, "MJ", "penrm"):
                self._penrm = value
                
    @property
    def penrt(self):
        return self._penrt
    
    @penrt.setter
    def penrt(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "penrt"):
            if self._check_unit(value.unit, "MJ", "penrt"):
                self._penrt = value
    
    @property
    def sm(self):
        return self._sm
    
    @sm.setter
    def sm(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "sm"):
            if self._check_unit(value.unit, "kg", "sm"):
                self._sm = value
    
    @property
    def rsf(self):
        return self._rsf
    
    @rsf.setter
    def rsf(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "rsf"):
            if self._check_unit(value.unit, "MJ", "rsf"):
                self._rsf = value
    
    @property
    def nrsf(self):
        return self._nrsf
    
    @nrsf.setter
    def nrsf(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "nrsf"):
            if self._check_unit(value.unit, "MJ", "nrsf"):
                self._nrsf = value
    
    @property
    def fw(self):
        return self._fw
    
    @fw.setter
    def fw(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "fw"):
            if self._check_unit(value.unit, "m^3", "fw"):
                self._fw = value
    
    @property
    def hwd(self):
        return self._hwd
    
    @hwd.setter
    def hwd(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "hwd"):
            if self._check_unit(value.unit, "kg", "hwd"):
                self._hwd = value
                
    @property
    def nhwd(self):
        return self._nhwd
    
    @nhwd.setter
    def nhwd(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "nhwd"):
            if self._check_unit(value.unit, "kg", "nhwd"):
                self._nhwd = value
                
    @property
    def rwd(self):
        return self._rwd
    
    @rwd.setter
    def rwd(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "rwd"):
            if self._check_unit(value.unit, "kg", "rwd"):
                self._rwd = value
                
    @property
    def cru(self):
        return self._cru
    
    @cru.setter
    def cru(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "cru"):
            if self._check_unit(value.unit, "kg", "cru"):
                self._cru = value
                
    @property
    def mfr(self):
        return self._mfr
    
    @mfr.setter
    def mfr(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "mfr"):
            if self._check_unit(value.unit, "kg", "mfr"):
                self._mfr = value
                
    @property
    def mer(self):
        return self._mer
    
    @mer.setter
    def mer(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "mer"):
            if self._check_unit(value.unit, "kg", "mer"):
                self._mer = value
                
    @property
    def eee(self):
        return self._eee
    
    @eee.setter
    def eee(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "eee"):
            if self._check_unit(value.unit, "MJ", "eee"):
                self._eee = value
                
    @property
    def eet(self):
        return self._eet
    
    @eet.setter
    def eet(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "eet"):
            if self._check_unit(value.unit, "MJ", "eet"):
                self._eet = value
                
    @property
    def gwp(self):
        return self._gwp
    
    @gwp.setter
    def gwp(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "gwp"):
            if self._check_unit(value.unit, "kg CO2 eq.", "gwp"):
                self._gwp = value
    
    @property
    def odp(self):
        return self._odp
    
    @odp.setter
    def odp(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "odp"):
            if self._check_unit(value.unit, "kg R11 eq.", "odp"):
                self._odp = value
    
    @property
    def pocp(self):
        return self._pocp
    
    @pocp.setter
    def pocp(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "pocp"):
            if self._check_unit(value.unit, "kg Ethene eq.", "pocp"):
                self._pocp = value
    
    @property
    def ap(self):
        return self._ap
    
    @ap.setter
    def ap(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "ap"):
            if self._check_unit(value.unit, "kg SO2 eq.", "ap"):
                self._ap = value
    
    @property
    def ep(self):
        return self._ep
    
    @ep.setter
    def ep(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "ep"):
            if self._check_unit(value.unit, "kg Phosphate eq.", "ep"):
                self._ep = value
    
    @property
    def adpe(self):
        return self._adpe
    
    @adpe.setter
    def adpe(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "adpe"):
            if self._check_unit(value.unit, "kg Sb eq.", "adpe"):
                self._adpe = value
    
    @property
    def adpf(self):
        return self._adpf
    
    @adpf.setter
    def adpf(self, value):
        if self._check_en15804indicatorvalue_class(value, var_name = "adpf"):
            if self._check_unit(value.unit, "MJ", "adpf"):
                self._adpf = value
    
    def set_values(self, **values):
        """Procedure to set all object-attributes at once

        Parameters
        ----------
        **values : dict
            Dictionary with attribute-name as key and attribute-value as value 

        """
        
        for attr, value in values.items():
            setattr(self, attr, value)
            
    def _ignore_none_mul(self, lca_indicator, scalar):
        if lca_indicator:
            return(lca_indicator * scalar)
        else:
            return None
            

    def __mul__(self, scalar):
        """Multiplies every indicator with a scalar
        

        Parameters
        ----------
        scalar : int, float
            scalar to be mutliplied

        Returns
        -------
        new : En15804LcaData
            product of "self" and scalar
        """       
        
        
        values = {"pere": self._ignore_none_mul(self.pere , scalar),
                "pert": self._ignore_none_mul(self.pert , scalar),
                "penre": self._ignore_none_mul(self.penre , scalar),
                "penrm": self._ignore_none_mul(self.penrm , scalar),
                "penrt": self._ignore_none_mul(self.penrt , scalar),
                "sm": self._ignore_none_mul(self.sm , scalar),
                "rsf": self._ignore_none_mul(self.rsf , scalar),
                "nrsf": self._ignore_none_mul(self.nrsf , scalar),
                "fw": self._ignore_none_mul(self.fw , scalar),
                "hwd": self._ignore_none_mul(self.hwd , scalar),
                "nhwd": self._ignore_none_mul(self.nhwd , scalar),
                "rwd": self._ignore_none_mul(self.rwd , scalar),
                "cru": self._ignore_none_mul(self.cru , scalar),
                "mfr": self._ignore_none_mul(self.mfr , scalar),
                "mer": self._ignore_none_mul(self.mer , scalar),
                "eee": self._ignore_none_mul(self.eee , scalar),
                "eet": self._ignore_none_mul(self.eet , scalar),
                "gwp": self._ignore_none_mul(self.gwp , scalar),
                "odp": self._ignore_none_mul(self.odp , scalar),
                "pocp": self._ignore_none_mul(self.pocp , scalar),
                "ap": self._ignore_none_mul(self.ap , scalar),
                "ep": self._ignore_none_mul(self.ep , scalar),
                "adpe": self._ignore_none_mul(self.adpe , scalar),
                "adpf": self._ignore_none_mul(self.adpf , scalar),
                "unit": self.ref_flow_unit
                }
        
        new = En15804LcaData()
        new.set_values(**values)
        
        return(new)

        
    def __add__(self, other):
        """function to add two En15804LcaData-Objects.

        Parameters
        ----------
        other : En15804LcaData
            other En15804LcaData-Object to be added

        Returns
        -------
        new : En15804LcaData
            sum of other En15804LcaData-Object and self.

        """
        
        if isinstance(other, En15804LcaData):

            values = {"pere": self.pere + other.pere,
                        "pert": self.pert + other.pert,
                        "penre": self.penre + other.penre,
                        "penrm": self.penrm + other.penrm,
                        "penrt": self.penrt + other.penrt,
                        "sm": self.sm + other.sm,
                        "rsf": self.rsf + other.rsf,
                        "nrsf": self.nrsf + other.nrsf,
                        "fw": self.fw + other.fw,
                        "hwd": self.hwd + other.hwd,
                        "nhwd": self.nhwd + other.nhwd,
                        "rwd": self.rwd + other.rwd,
                        "cru": self.cru + other.cru,
                        "mfr": self.mfr + other.mfr,
                        "mer": self.mer + other.mer,
                        "eee": self.eee + other.eee,
                        "eet": self.eet + other.eet,
                        "gwp": self.gwp + other.gwp,
                        "odp": self.odp + other.odp,
                        "pocp": self.pocp + other.pocp,
                        "ap": self.ap + other.ap,
                        "ep": self.ep + other.ep,
                        "adpe": self.adpe + other.adpe,
                        "adpf": self.adpf + other.adpf    
                    }
                        
            new = En15804LcaData()
            new.set_values(**values)
            return new
        else:      
            print("Addend must be an 'En15804LcaData'-Object!")
            
    def load_lca_data_template(self, lca_id, data_class=None):
        """LCA-data loader.

        Loads LCA-data specified in the json.

        Parameters
        ----------

        lca_id : str
            LCA-data Identifier

        data_class : DataClass
            DataClass containing the bindings for LCA-data and LCA-data
            -fallbacks (typically this is the data class stored in prj.data,
            but the user can individually change that.)

        """
        
        if data_class == None:
            data_class = self.parent.parent.parent.parent.data

        lca_data_input.load_en15804_lca_data_id(lca_data=self,
                                     lca_id=lca_id,
                                     data_class=data_class)
        
    def save_lca_data_template(self, data_class=None):
        """LCA-data loader.

        Saves LCA-data specified in the json.

        Parameters
        ----------

        data_class : DataClass
            DataClass containing the bindings for LCA-data and LCA-data
            -fallbacks (typically this is the data class stored in prj.data,
            but the user can individually change that.)

        """
        
        if data_class is None:
            data_class = self.parent.parent.parent.parent.data
        

        lca_data_output.save_lca_data(lca_data=self,
                                     data_class=data_class)
        
    def load_fallbacks(self, fallback_dictonarie, data_class):
        """procedure to load fallbacks specified in the fallback-dictonarie
        into the self.fallback attribute

        Parameters
        ----------
        fallback_dictonarie : dict
            Dictonarie with stages as key and LCA-dataset-uuid as value
        data_class : DataClass
            DataClass containing the bindings for LCA-data and LCA-data
            -fallbacks (typically this is the data class stored in prj.data,
            but the user can individually change that.)

        Returns
        -------
        None.

        """
        if fallback_dictonarie:
            self.fallback = {}
            for stage in fallback_dictonarie:
                
                fallback_id = fallback_dictonarie[stage]
                
                fallback_object = En15804LcaData(self.parent)
                
                self.fallback[stage] = fallback_object
                
                self.fallback[stage].load_lca_data_fallback_template(fallback_id, data_class)
    
    def load_lca_data_fallback_template(self, lca_id, data_class=None):
        """LCA-data-fallback loader.

        Loads LCA-data-fallbacks specified in the json.

        Parameters
        ----------

        lca_id : str
            LCA-Data Identifier

        data_class : DataClass
            DataClass containing the bindings for LCA-data and LCA-data
            -fallbacks (typically this is the data class stored in prj.data,
            but the user can individually change that.)

        """

        if data_class == None:
            data_class = self.parent.parent.parent.parent.data

        lca_data_input.load_en15804_lca_data_fallback_id(lca_data=self,
                                     lca_id=lca_id,
                                     data_class=data_class)
        
    def convert_ref_unit(self, target_unit, area = None, thickness = None, density = None):
        """converts the values of the environmental indicators to a new 
        reference unit. All parameters are optional and are only used if they 
        are necessary for the conversion.
        

        Parameters
        ----------
        target_unit : str
            target unit of the conversion
        area : float, optional
            area (e.g. area of a buildingelement)
        thickness : float, optional
            thickness (e.g. thickness of a layer). The default is None.
        density : float, optional
            density (e.g. the density of a material). The default is None.

        Returns
        -------
        result : En15804LcaData
            LCA-Data with new reference flow

        """
        result = En15804LcaData()
        
        if target_unit == "pcs":

            if self.ref_flow_unit == "m^3":
                scalar = (area * thickness)
                
            elif self.ref_flow_unit == "kg":
                scalar = (area * thickness * density)

            elif self.ref_flow_unit == "m^2":
                scalar = area
            
            elif self.ref_flow_unit == "pcs":
                scalar = 1
                
            else:
                scalar = 1
                target_unit = self.ref_flow_unit
                raise ValueError("Unable to convert {} into {}!".format(self.ref_flow_unit, target_unit))
        
        elif target_unit == "kg":
            
            if self.ref_flow_unit == "m^3":
                scalar = 1 / density
            elif self.ref_flow_unit == "m^2":
                scalar = thickness / density
            
            else:
                scalar = 1
                target_unit = self.unit
                raise ValueError("Unable to convert unit into target unit!")
        
        elif target_unit == "m^2":
            if self.ref_flow_unit == "m^3":
                scalar = 1/thickness
            
        scalar = scalar / self.ref_flow_value
        result = self * scalar
        result.ref_flow_unit = target_unit
        result.ref_flow_value = 1
        
        return result
            
    def sum_to_b4(self):
        """function to sum up all every stage indicators to stage
        B4 'replacement'
        

        Returns
        -------
        result : En15804LcaData()
            LCA-data with all stages sumed up to B4

        """
        
        result = En15804LcaData()
        
        values = {"pere": En15804IndicatorValue() ,
                "pert": En15804IndicatorValue(),
                "penre": En15804IndicatorValue(),
                "penrm": En15804IndicatorValue(),
                "penrt": En15804IndicatorValue(),
                "sm": En15804IndicatorValue(),
                "rsf": En15804IndicatorValue(),
                "nrsf": En15804IndicatorValue(),
                "fw": En15804IndicatorValue(),
                "hwd": En15804IndicatorValue(),
                "nhwd": En15804IndicatorValue(),
                "rwd": En15804IndicatorValue(),
                "cru": En15804IndicatorValue(),
                "mfr": En15804IndicatorValue(),
                "mer": En15804IndicatorValue(),
                "eee": En15804IndicatorValue(),
                "eet": En15804IndicatorValue(),
                "gwp": En15804IndicatorValue(),
                "odp": En15804IndicatorValue(),
                "pocp": En15804IndicatorValue(),
                "ap": En15804IndicatorValue(),
                "ep": En15804IndicatorValue(),
                "adpe": En15804IndicatorValue(),
                "adpf": En15804IndicatorValue(),
                "ref_flow_value": self.ref_flow_value,
                "ref_flow_unit": self.ref_flow_unit
            }
        
        values["pere"].unit = self.pere.unit
        values["pert"].unit = self.pert.unit
        values["penre"].unit = self.penre.unit
        values["penrm"].unit = self.penrm.unit
        values["penrt"].unit = self.penrt.unit
        values["sm"].unit = self.sm.unit
        values["rsf"].unit = self.rsf.unit
        values["nrsf"].unit = self.nrsf.unit
        values["fw"].unit = self.fw.unit
        values["hwd"].unit = self.hwd.unit
        values["nhwd"].unit = self.nhwd.unit
        values["rwd"].unit = self.rwd.unit
        values["cru"].unit = self.cru.unit
        values["mfr"].unit = self.mfr.unit
        values["mer"].unit = self.mer.unit
        values["eee"].unit = self.eee.unit
        values["eet"].unit = self.eet.unit
        values["gwp"].unit = self.gwp.unit
        values["odp"].unit = self.odp.unit
        values["pocp"].unit = self.pocp.unit
        values["ap"].unit = self.ap.unit
        values["ep"].unit = self.ep.unit
        values["adpe"].unit = self.adpe.unit
        values["adpf"].unit = self.adpf.unit

        
        values["pere"].b4 = self.pere.sum_stages()
        values["pert"].b4 = self.pert.sum_stages()
        values["penre"].b4 = self.penre.sum_stages()
        values["penrm"].b4 = self.penrm.sum_stages()
        values["penrt"].b4 = self.penrt.sum_stages()
        values["sm"].b4 = self.sm.sum_stages()
        values["rsf"].b4 = self.rsf.sum_stages()
        values["nrsf"].b4 = self.nrsf.sum_stages()
        values["fw"].b4 = self.fw.sum_stages()
        values["hwd"].b4 = self.hwd.sum_stages()
        values["nhwd"].b4 = self.nhwd.sum_stages()
        values["rwd"].b4 = self.rwd.sum_stages()
        values["cru"].b4 = self.cru.sum_stages()
        values["mfr"].b4 = self.mfr.sum_stages()
        values["mer"].b4 = self.mer.sum_stages()
        values["eee"].b4 = self.eee.sum_stages()
        values["eet"].b4 = self.eet.sum_stages()
        values["gwp"].b4 = self.gwp.sum_stages()
        values["odp"].b4 = self.odp.sum_stages()
        values["pocp"].b4 = self.pocp.sum_stages()
        values["ap"].b4 = self.ap.sum_stages()
        values["ep"].b4 = self.ep.sum_stages()
        values["adpe"].b4 = self.adpe.sum_stages()
        values["adpf"].b4 = self.adpf.sum_stages()
        
        result.set_values(**values)
        return result
    
    def add_fallbacks(self):
        """adds the indicators from the lca-data-fallbacks specified in 
        self.fallback to the matching stages of self.The attribute 
        "_fallback_added" is True, if the fallbacks are allready added

        Returns
        -------
        None.

        """
        
        if self._fallback_added is False:
            
            pere_backup = self.pere
            pert_backup = self.pert
            penre_backup = self.penre
            penrm_backup = self.penrm
            penrt_backup = self.penrt
            sm_backup = self.sm
            rsf_backup = self.rsf
            nrsf_backup = self.nrsf
            fw_backup = self.fw
            hwd_backup = self.hwd
            nhwd_backup = self.nhwd
            rwd_backup = self.rwd
            cru_backup = self.cru
            mfr_backup = self.mfr
            mer_backup = self.mer
            eee_backup = self.eee
            eet_backup = self.eet
            gwp_backup = self.gwp
            odp_backup = self.odp
            pocp_backup = self.pocp
            ap_backup = self.ap
            ep_backup = self.ep
            adpe_backup = self.adpe
            adpf_backup = self.adpf


            self._fallback_added = True
            
            for stage in self.fallback:
                
                if self.fallback[stage].ref_flow_unit != self.ref_flow_unit:
                    try:
                    
                        self.fallback[stage].convert_ref_unit(
                            target_unit = self.ref_flow_unit,
                            density = self.parent.density,
                            thickness = self.parent.parent.thickness                           
                            )
                    
                    except:
                        
                        print("Error while trying to convert fallback {} of {} reference unit".format(self.fallback[stage].lca_data_id, self.lca_data_id),
                             "{} to {}".format( self.fallback[stage].ref_flow_unit, self.ref_flow_unit))
                        
                        self.pere = pere_backup
                        self.pert = pert_backup
                        self.penre = penre_backup
                        self.penrm = penrm_backup
                        self.penrt = penrt_backup
                        self.sm = sm_backup
                        self.rsf = rsf_backup
                        self.nrsf = nrsf_backup
                        self.fw = fw_backup
                        self.hwd = hwd_backup
                        self.nhwd = nhwd_backup
                        self.rwd = rwd_backup
                        self.cru = cru_backup
                        self.mfr = mfr_backup
                        self.mer = mer_backup
                        self.eee = eee_backup
                        self.eet = eet_backup
                        self.gwp = gwp_backup
                        self.odp = odp_backup
                        self.pocp = pocp_backup
                        self.ap = ap_backup
                        self.ep = ep_backup
                        self.adpe = adpe_backup
                        self.adpf = adpf_backup
                        
                        self.fallback_added = False
                        break
                      

                    
                
                self.pere = self.pere.add_stage(stage, self.fallback[stage].pere)
                self.pert = self.pert.add_stage(stage, self.fallback[stage].pert)
                self.penre = self.penre.add_stage(stage, self.fallback[stage].penre)
                self.penrm = self.penrm.add_stage(stage, self.fallback[stage].penrm)
                self.penrt = self.penrt.add_stage(stage, self.fallback[stage].penrt)
                self.sm = self.sm.add_stage(stage, self.fallback[stage].sm)
                self.rsf = self.rsf.add_stage(stage, self.fallback[stage].rsf)
                self.nrsf = self.nrsf.add_stage(stage, self.fallback[stage].nrsf)
                self.fw = self.fw.add_stage(stage, self.fallback[stage].fw)
                self.hwd = self.hwd.add_stage(stage, self.fallback[stage].hwd)
                self.nhwd = self.nhwd.add_stage(stage, self.fallback[stage].nhwd)
                self.rwd = self.rwd.add_stage(stage, self.fallback[stage].rwd)
                self.cru = self.cru.add_stage(stage, self.fallback[stage].cru)
                self.mfr = self.mfr.add_stage(stage, self.fallback[stage].mfr)
                self.mer = self.mer.add_stage(stage, self.fallback[stage].mer)
                self.eee = self.eee.add_stage(stage, self.fallback[stage].eee)
                self.eet = self.eet.add_stage(stage, self.fallback[stage].eet)
                self.gwp = self.gwp.add_stage(stage, self.fallback[stage].gwp)
                self.odp = self.odp.add_stage(stage, self.fallback[stage].odp)
                self.pocp = self.pocp.add_stage(stage, self.fallback[stage].pocp)
                self.ap = self.ap.add_stage(stage, self.fallback[stage].ap)
                self.ep = self.ep.add_stage(stage, self.fallback[stage].ep)
                self.adpe = self.adpe.add_stage(stage, self.fallback[stage].adpe)
                self.adpf = self.adpf.add_stage(stage, self.fallback[stage].adpf)
                

            