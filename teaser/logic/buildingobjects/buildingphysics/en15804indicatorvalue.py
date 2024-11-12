# -*- coding: utf-8 -*-

# created September 2021
# by Johannes Linus Cuypers

class En15804IndicatorValue(object):
    """En15804IndicatorValue class
    
    This class holds one value for every lifecycle-stage according to EN 15804.
    It can be used to set the value of an environmental indicator.    
    
    Parameters
    ----------
    
    parent : XXX, optional
        Default is None
    
    Attributes
    ----------
    unit : str
        unit of the values
    a1 : float
        Product stage: Raw material supply
    a2 : float
        Product stage: Transport
    a3 : float
        Product stage: Manufacturing
    a1_a3 : float
        sum of stage a1, a2 and a3
    a4 : float
        Construction on process stage: Transport from the gate to the site
    a5 : float
        Construction on process stage: Assembly
    b1 : float
        Use stage: Use
    b2 : float
        Use Stage: Maintenance
    b3 : float
        Use Stage: Repair
    b4 : float
        Use Stage: Replacement
    b5 : float
        Use Stage: Refurbishment
    b6 : float
        Use Stage: Operational energy use
    b7 : float
        Use Stage: Operational water use
    c1 : float
        End of Life Stage: De-Construction demolition
    c2 : float
        End of Life Stage: Transport
    c3 : float
        End of Life Stage: Waste processing
    c4 : float
        End of Life stage: Disposal
    d : float
        Benefits and loads beyond the system boundaries: 
        Reuse-Recovery-Recycling-Potential
    
    
    """
    def __init__(self, parent = None):
        """Constructor of En15804IndicatorValue
        """
        
        self.parent = parent
        self._unit = None
        
        self._a1 = None
        self._a2 = None
        self._a3 = None
        self._a1_a3 = None
        
        self._a4 = None
        self._a5 = None
        
        self._b1 = None
        self._b2 = None
        self._b3 = None
        self._b4 = None
        self._b5 = None
        self._b6 = None
        self._b7 = None
        
        self._c1 = None
        self._c2 = None
        self._c3 = None
        self._c4 = None
        
        self._d = None
        
    
        
    def _validate_stage_value(self, value, stage_name):
        """Function to validate the value of an stage.

        Parameters
        ----------
        value : float
            value to be validated 
        stage_name : str
            printable name of a stage. Used for Error-Messages.

        Returns
        -------
        validated value: float, nonetype

        """
        if isinstance(value, float):
            return(value)
        else:
            if not value:
                return(None)
            try:
                value = float(value)
                return(value)
            except ValueError:
                print("Can't convert value of '{}' to float".format(stage_name))       
            except TypeError:
                print("Can´t convert {} into float".format(type(value)))

                
    @property
    def unit(self):
        return(self._unit)
    
    @unit.setter
    def unit(self, value):
        if value != None:
            if isinstance(value, str):
                self._unit = value
            else:
                print("En15804IndicatorValue.unit must be string!")
        else:
            self._unit = None

    @property
    def a1(self):
        return(self._a1)
    
    @a1.setter
    def a1(self, value):
       self._a1 = self._validate_stage_value(value, "A1")
       
    @property
    def a2(self):
        return(self._a2)
    
    @a2.setter
    def a2(self, value):
       self._a2 = self._validate_stage_value(value, "A2")
       
    @property
    def a3(self):
        return(self._a3)
    
    @a3.setter
    def a3(self, value):
       self._a3 = self._validate_stage_value(value, "A3")
    
    @property
    def a1_a3(self):
        return(self._a1_a3)
    
    @a1_a3.setter
    def a1_a3(self, value):
       self._a1_a3 = self._validate_stage_value(value, "A1-A3")
       
    @property
    def a4(self):
        return(self._a4)
    
    @a4.setter
    def a4(self, value):
       self._a4 = self._validate_stage_value(value, "A4")
       
    @property
    def a5(self):
        return(self._a5)
    
    @a5.setter
    def a5(self, value):
       self._a5 = self._validate_stage_value(value, "A5")
       
    @property
    def b1(self):
        return(self._b1)
    
    @b1.setter
    def b1(self, value):
       self._b1 = self._validate_stage_value(value, "B1")
       
    @property
    def b2(self):
        return(self._b2)
    
    @b2.setter
    def b2(self, value):
       self._b2 = self._validate_stage_value(value, "B2")
       
    @property
    def b3(self):
        return(self._b3)
    
    @b3.setter
    def b3(self, value):
       self._b3 = self._validate_stage_value(value, "B3")
       
    @property
    def b4(self):
        return(self._b4)
    
    @b4.setter
    def b4(self, value):
       self._b4 = self._validate_stage_value(value, "B4")
       
    @property
    def b5(self):
        return(self._b5)
    
    @b5.setter
    def b5(self, value):
       self._b5 = self._validate_stage_value(value, "B5")
       
    @property
    def b6(self):
        return(self._b6)
    
    @b6.setter
    def b6(self, value):
       self._b6 = self._validate_stage_value(value, "B6")
       
    @property
    def b7(self):
        return(self._b7)
    
    @b7.setter
    def b7(self, value):
       self._b7 = self._validate_stage_value(value, "B7")
       
    @property
    def c1(self):
        return(self._c1)
    
    @c1.setter
    def c1(self, value):
       self._c1 = self._validate_stage_value(value, "C1")
       
    @property
    def c2(self):
        return(self._c2)
    
    @c2.setter
    def c2(self, value):
       self._c2 = self._validate_stage_value(value, "C2")
       
    @property
    def c3(self):
        return(self._c3)
    
    @c3.setter
    def c3(self, value):
       self._c3 = self._validate_stage_value(value, "C3")
       
    @property
    def c4(self):
        return(self._c4)
    
    @c4.setter
    def c4(self, value):
       self._c4 = self._validate_stage_value(value, "C4")
       
    @property
    def d(self):
        return(self._d)
    
    @d.setter
    def d(self, value):
       self._d = self._validate_stage_value(value, "D")
       
    def set_values(self, **values):
        """Procedure to set all object-attributes at once

        Parameters
        ----------
        **values : dict
            Dictionary with attribute-name as key and attribute-value as value 

        """
        
        if "a1_a3" not in values or values["a1_a3"] is None:
            if "a1" in values and "a2" in values and "a3" in values:
                if values["a1"] is not None or values["a2"] is not None or values["a3"] is not None:
                    a1_a3 = self._ignore_none_sum(self._ignore_none_sum(values["a1"], values["a2"]), values["a3"])
                    
                    values["a1_a3"] = a1_a3

        for attr, value in values.items():
            setattr(self, attr, value)

    
    def _ignore_none_sum(self, addend1, addend2):
        """sums up two numbers, but interprets None addends as 0. When both
        input-addends are None, it returns None
        
        Parameters
        ----------
        addend1 : int, float, none
            addend 1
        addend2 : int, float, none
            addend 2

        Returns
        -------
        sum of addends : int, float
        or None
        """
        if addend1 and addend2:
            return(addend1 + addend2)
        else:
            if not addend1 and not addend2:
                return(None)
            else:
                if addend1:
                    return(addend1)
                else:
                    return(addend2)
                
    def __add__(self, other):
        """Adds two En15804IndicatorValue-Objects. Every stage is summed up 
        separately. Both objects must have the same unit
        

        Parameters
        ----------
        other : En15804IndicatorValue
            The addend En15804IndicatorValue-Object

        Returns
        -------
        sum: En15804IndicatorValue
           sum of both En15804IndicatorValue-Objects. 

        """
        if isinstance(other, En15804IndicatorValue):
            
            if self.unit == other.unit or not self.unit or not other.unit:

                values = {"a1": self._ignore_none_sum(self.a1, other.a1),
                    "a2": self._ignore_none_sum(self.a2, other.a2),
                    "a3": self._ignore_none_sum(self.a3, other.a3),
                    "a1_a3": self._ignore_none_sum(self.a1_a3, other.a1_a3),
                    "a4": self._ignore_none_sum(self.a4, other.a4),
                    "a5": self._ignore_none_sum(self.a5, other.a5),
                    "b1": self._ignore_none_sum(self.b1, other.b1),
                    "b2": self._ignore_none_sum(self.b2, other.b2),
                    "b3": self._ignore_none_sum(self.b3, other.b3),
                    "b4": self._ignore_none_sum(self.b4, other.b4),
                    "b5": self._ignore_none_sum(self.b5, other.b5),
                    "b6": self._ignore_none_sum(self.b6, other.b6),
                    "b7": self._ignore_none_sum(self.b7, other.b7),
                    "c1": self._ignore_none_sum(self.c1, other.c1),
                    "c2": self._ignore_none_sum(self.c2, other.c2),
                    "c3": self._ignore_none_sum(self.c3, other.c3),
                    "c4": self._ignore_none_sum(self.c4, other.c4),
                    "d": self._ignore_none_sum(self.d, other.d),
                    }
                
                if self.unit:
                    values["unit"] = self.unit
                else:
                    values["unit"] = other.unit
                
                new = En15804IndicatorValue()
                new.set_values(**values)
                return(new)
            
            else:
                print("Addends must have the same unit!")
            
        else:
            
            print("Addend must be an 'En15804IndicatorValue'-Object!")
            
    def _ignore_none_mul(self, factor1, factor2):
        """Multiplies two factors. If one is None, it returns None

        Parameters
        ----------
        factor1 : int, float
            Factor 1
        factor2 : TYPE
            Factor 2

        Returns
        -------
            Product of factors or None

        """
        if not factor1 or not factor2:
            return(None)
        else:
            return(factor1 * factor2)
        
    def __mul__(self, scalar):
        """Multiplies every stage-value with a scalar

        Parameters
        ----------
        scalar : int, float
            scalar to be mutliplied

        Returns
        -------
        new : En15804IndicatorValue
            product of "self" and scalar

        """
        try:
            scalar = float(scalar)
            
            values = {"a1": self._ignore_none_mul(self.a1, scalar),
                "a2": self._ignore_none_mul(self.a2, scalar),
                "a3": self._ignore_none_mul(self.a3, scalar),
                "a1_a3": self._ignore_none_mul(self.a1_a3, scalar),
                "a4": self._ignore_none_mul(self.a4, scalar),
                "a5": self._ignore_none_mul(self.a5, scalar),
                "b1": self._ignore_none_mul(self.b1, scalar),
                "b2": self._ignore_none_mul(self.b2, scalar),
                "b3": self._ignore_none_mul(self.b3, scalar),
                "b4": self._ignore_none_mul(self.b4, scalar),
                "b5": self._ignore_none_mul(self.b5, scalar),
                "b6": self._ignore_none_mul(self.b6, scalar),
                "b7": self._ignore_none_mul(self.b7, scalar),
                "c1": self._ignore_none_mul(self.c1, scalar),
                "c2": self._ignore_none_mul(self.c2, scalar),
                "c3": self._ignore_none_mul(self.c3, scalar),
                "c4": self._ignore_none_mul(self.c4, scalar),
                "d": self._ignore_none_mul(self.d, scalar),
                "unit": self.unit}
            
            new = En15804IndicatorValue()
            new.set_values(**values)
            return(new)
            
        except ValueError:
            print("Can't convert value of '{}' to float. Please insert scalar!".format(scalar))
        
        except TypeError:
            print("Can´t convert {} into float. Please insert scalar!".format(type(scalar)))
    
    def sum_stages(self, add_stage_d = False):
        """sums up all values of the object

        Parameters
        ----------
        add_stage_d : Boolean
            Is true, when stage d should be added. The default is False.

        Returns
        -------
            sum : float
                overall value of all stages. With or without stage d 

        """

        addends = []
        
        
        
        if self.a1_a3:
            addends.append(self.a1_a3)
        else:
            if self.a1: addends.append(self.a1)
            if self.a2: addends.append(self.a2)
            if self.a3: addends.append(self.a3)
        if self.a4: addends.append(self.a4)
        if self.a5: addends.append(self.a5)
        if self.b1: addends.append(self.b1)
        if self.b2: addends.append(self.b2)
        if self.b3: addends.append(self.b3)
        if self.b4: addends.append(self.b4)
        if self.b5: addends.append(self.b5)
        if self.b6: addends.append(self.b6)
        if self.b7: addends.append(self.b7)
        if self.c1: addends.append(self.c1)
        if self.c2: addends.append(self.c2)
        if self.c3: addends.append(self.c3)
        if self.c4: addends.append(self.c4)
        if self.d and add_stage_d: addends.append(self.d)
        
        return(sum(addends))

    def stages(self):
        """returns a list with all stages that are not None


        Returns
        -------
            list of stages

        """
        stages = []
        
        if self.a1 is not None: stages.append("a1")
        if self.a2 is not None: stages.append("a2")
        if self.a3 is not None: stages.append("a3")
        if self.a1_a3 is not None: stages.append("a1_a3")
        if self.a4 is not None: stages.append("a4")
        if self.a5 is not None: stages.append("a5")
        if self.b1 is not None: stages.append("b1")
        if self.b2 is not None: stages.append("b2")
        if self.b3 is not None: stages.append("b3")
        if self.b4 is not None: stages.append("b4")
        if self.b5 is not None: stages.append("b5")
        if self.b6 is not None: stages.append("b6")
        if self.b7 is not None: stages.append("b7")
        if self.c1 is not None: stages.append("c1")
        if self.c2 is not None: stages.append("c2")
        if self.c3 is not None: stages.append("c3")
        if self.c4 is not None: stages.append("c4")
        if self.d is not None: stages.append("d")
        
        return stages

    def add_stage(self, stage, other):
        """Function which adds only a specific stage from an En15804Indicator-
        Object to self
        

        Parameters
        ----------
        stage : str
            Stage to add (e.g. 'b4')
        other : En15804IndicatorValue
            the other En15804IndicatorValue

        Returns
        -------
        result En15804IndicatorValue

        """
        result = En15804IndicatorValue()
        addend = En15804IndicatorValue()
        addend.unit = other.unit
        if isinstance(other, En15804IndicatorValue):
            other_stages = other.stages()
            if stage in other_stages:
                if len(other_stages) != 1:
                    if stage == "a1":
                        addend.a1 = other.a1
                    elif stage == "a2":
                        addend.a2 = other.a2
                    elif stage == "a3":
                        addend.a3 = other.a3
                    elif stage == "a1_a3":
                        addend.a1_a3 = other.a1_a3
                    elif stage == "a4":
                        addend.a4 = other.a4
                    elif stage == "a5":
                        addend.a5 = other.a5
                    elif stage == "b1":
                        addend.b1 = other.b1
                    elif stage == "b2":
                        addend.b2 = other.b2
                    elif stage == "b3":
                        addend.b3 = other.b3
                    elif stage == "b4":
                        addend.b4 = other.b4
                    elif stage == "b5":
                        addend.b5 = other.b5
                    elif stage == "b6":
                        addend.b6 = other.b6
                    elif stage == "b7":
                        addend.b7 = other.b7
                    elif stage == "c1":
                        addend.c1 = other.c1
                    elif stage == "c2":
                        addend.c2 = other.c2
                    elif stage == "c3":
                        addend.c3 = other.c3
                    elif stage == "c4":
                        addend.c4 = other.c4
                    elif stage == "d":
                        addend.d = other.d

                else:
                    addend = other

                    
                    #print(vars(other))
                    #other = other.set_values(**{stage: vars(other)[stage]})
                result = self + addend
            else:
                print("No value for Stage {}".format(stage))
                result = self
        
        else:
            print("Addend must be an 'En15804IndicatorValue'-Object!")
                
        
        return(result)
    
    def get_values_as_dict(self):
        return_dict = {}
        
        if self.a1 is not None:
            return_dict["a1"] = self.a1
        if self.a2 is not None:
            return_dict["a2"] = self.a2
        if self.a3 is not None:
            return_dict["a3"] = self.a3
        if self.a1_a3 is not None:
            return_dict["a1_a3"] = self.a1_a3
        if self.a4 is not None:
            return_dict["a4"] = self.a4
        if self.a5 is not None:
            return_dict["a5"] = self.a5
            
        if self.b1 is not None:
            return_dict["b1"] = self.b1
        if self.b2 is not None:
            return_dict["b2"] = self.b2
        if self.b3 is not None:
            return_dict["b3"] = self.b3
        if self.b4 is not None:
            return_dict["b4"] = self.b4
        if self.b5 is not None:
            return_dict["b5"] = self.b5
        if self.b6 is not None:
            return_dict["b6"] = self.b6
        if self.b7 is not None:
            return_dict["b7"] = self.b7
            
        if self.c1 is not None:
            return_dict["c1"] = self.c1
        if self.c2 is not None:
            return_dict["c2"] = self.c2
        if self.c3 is not None:
            return_dict["c3"] = self.c3
        if self.c4 is not None:
            return_dict["c4"] = self.c4
        if self.d is not None:
            return_dict["d"] = self.d

        return_dict["unit"] = self.unit
        
        return(return_dict)
                    

           
                
                    
            
            
       
                
                
    
        