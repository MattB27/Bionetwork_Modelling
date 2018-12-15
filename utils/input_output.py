"""
modelling.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A module to read sbml models and vonstruct gekko models.
"""

import pandas as pd
import numpy as np
from libsbml import *

class bionetwork_model:
    def __init__(self, filepath):
        doc = readSBML(filepath)
        model = doc.getModel()
        species = model.getListOfSpecies()
        params = model.getListOfParameters()
        reactions = model.getListOfReactions()
        rules = model.getListOfRules()
        self.CV = []
        self.FV = []
        self.MV = []
        self.SV = []
        self.model = model
        
        species_id = [x.getId() for x in species]
        species_name = [x.getName() for x in species]
        species_value = [x.getInitialConcentration() if x.getInitialAmount() == 0\
                         else x.getInitialAmount() for x in species]
        self.species = pd.DataFrame({"ID": species_id, "Name": species_name, "Value": species_value})
        variables = pd.DataFrame({"ID": species_id, "Name": species_name, "Value": species_value})
        variables = variables.set_index("ID")
        self.variables = variables
 
        param_id = [x.getId() for x in params]
        param_name = [x.getName() for x in params]
        param_type = [x.getConstant() for x in params]
        param_value = [x.getValue() if x.getConstant() else np.nan for x in params]
        for rx in reactions:
            if rx.getKineticLaw().getListOfLocalParameters is not None:
                local_params = rx.getKineticLaw().getListOfLocalParameters()
            for param in local_params:
                param_id.append(param.getId())
                param_name.append(rx.getId()+'_'+param.getId())
                param_type.append(True)
                param_value.append(param.getValue())
        self.parameters = pd.DataFrame({"ID": param_id, "Name": param_name, "Constant":param_type, "Value": param_value})
        
        rules_data = [[rule.getId(),rule.getElementName(),rule.getFormula()] for rule in rules]
        rules_df = pd.DataFrame(rules_data, columns=['ID','Type','Eqn'])
        self.rules = rules_df
        
        rate_laws = []
        for rxn in reactions:
            rate_law_string = rxn.getKineticLaw().getFormula()
            local_params = rxn.getKineticLaw().getListOfLocalParameters()
            if local_params is not None:
                for param in local_params:
                    new_param_string = rxn.getId() + '_' + param.getId()
                    string_length = len(list(param.getId()))
                    new_string_length = len(new_param_string)
                    i = 0
                    while (i < len(rate_law_string)):
                        if ''.join(rate_law_string[i:(i+string_length)]) == param.getId():
                            join_list = [rate_law_string[:i], new_param_string, rate_law_string[(i+string_length):]]
                            rate_law_string = ''.join(join_list)
                            i += new_string_length
                        else:
                            i += 1
            if '* pow' in rate_law_string:
                num_pow = len(rate_law_string.split(sep='* pow')) - 1
                for i in range(num_pow):
                    manip_string = rate_law_string.split(sep='* pow')[1]
                    base, power = manip_string.split(sep=',')
                    new_string = '* ' + base[1:] + '*m.exp(' + power[1:]
                    replace_string = '* pow' + manip_string
                    rate_law_string = rate_law_string.replace(replace_string, new_string)
            rate_laws.append(rate_law_string)
        kinetic_laws = [rxn.getId() for rxn in reactions]
        kinetic_df = pd.DataFrame({"ID": kinetic_laws, "Eqn": rate_laws})
        self.kinetic = kinetic_df
        
        storage = np.zeros([len(reactions),len(species)])
        species_names = [x.getId() for x in species]
        reaction_names = [rxn.getId() for rxn in reactions]
        df = pd.DataFrame(storage,columns=species_names, index=reaction_names)
        for reaction in reactions:
            reaction_name = reaction.getId()
            reaction_series = pd.Series(np.zeros(len(species)),name=reaction_name,index=species_names)
            for reactant in reaction.getListOfReactants():
                reactant_series = pd.Series(-reactant.getStoichiometry(),index=[reactant.getSpecies()])
                reaction_series.update(reactant_series)
            for product in reaction.getListOfProducts():
                product_series = pd.Series(product.getStoichiometry(),index=[product.getSpecies()])
                reaction_series.update(product_series)
            df.loc[reaction_name].update(reaction_series)
        self.stoich_table = df
        
        
    def change_var_type(self,var_id,var_type,**kwargs):
        # grab value from variable df
        if var_type == "CV":
            self.CV.append(CV(self.variables.loc[var_id]["Value"], **kwargs))
        if var_type == "SV":
            self.SV.append(SV(self.variables.loc[var_id]["Value"], **kwargs))
        if var_type == "MV":
            self.MV.append(MV(self.variables.loc[var_id]["Value"], **kwargs))
        if var_type == "FV":
            self.FV.append(FV(self.variables.loc[var_id]["Value"], **kwargs))
        new_variables = self.variables
        new_variables = new_variables.drop(var_id)
        self.variables = new_variables
        return
        
    def write_global_options(self, option_dict):
        string = "# Global Options\n"
        for option, setting in option_dict.items():
            string += "m.options." + option + " = " + str(setting) + "\n"
        string += "m.GUI = True\n"
        string += "\n"
        
        return string
    
    def write_gekko(self, time):
        string = "''' SBML file converted to Gekko Format\n"
        string += 'This model assumes a time basis of seconds, events occuring in sequential order and ignores units.\n'
        string += "'''\n" + 'import numpy as np\nfrom gekko import GEKKO\nimport matplotlib.pyplot as plt\n'
        string += '\nm = GEKKO()\nm.time = np.linspace(0,' + str(time) + ',' + str(time+1) + ')\n'
        
        return string
    
    def write_var(self):
        string = "# Variables\n"
        for var in self.variables.index.values:
            string += var + ' = m.Var(value = ' + str(self.variables.loc[var]['Value']) + ')\n'
        string += '\n'
        return string
        
    def write_parameters(self):
        string = "# Parameters\n"
        parameter_bool = (self.parameters['Constant'] == True)
        parameter_df = self.parameters[parameter_bool]
        for i in range(len(parameter_df)):
            string += parameter_df['ID'].iloc[i] + ' = m.Param(value = ' + str(parameter_df['Value'].iloc[i]) + ')\n'
        string += '\n'
        return string
        
    def write_intermediates(self):
        string = '# Intermediates\n'
        for i in range(len(self.rules)):
            if self.rules['Type'].iloc[i] == 'assignmentRule':
                rule_string = self.rules['Eqn'].iloc[i]
                if '* pow' in rule_string:
                    num_pow = len(rule_string.split(sep='* pow')) - 1
                    for i in range(num_pow):
                        manip_string = rule_string.split(sep='* pow')[1]
                        base, power = manip_string.split(sep=',')
                        new_string = '* ' + base[1:] + '*m.exp(' + power[1:]
                        replace_string = '* pow' + manip_string
                        rule_string = rule_string.replace(replace_string, new_string)
                string += self.rules['ID'].iloc[i] + ' = m.Intermediate(' + rule_string + ')\n'
        for i in range(len(self.kinetic)):
            string += self.kinetic['ID'].iloc[i] + ' = m.Intermediate(' + self.kinetic['Eqn'].iloc[i] + ')\n'
        string += '\n'
        return string
    
    def write_equations(self):
        string1 = '# Equations\n'
        for column in self.stoich_table:
            string = 'm.Equation(' + column + '.dt() == ('
            for index in self.stoich_table.index:
                if self.stoich_table[column][index] != 0:
                    string += '(' + str(self.stoich_table[column][index]) + ' * ' + index + ') + '
            if string == 'm.Equation(' + column + '.dt() == (':
                string += '0))\n'
            else:
                string = ''.join(list(string[:-3]))
                string += '))\n'
            string1 += string
        for i in range(len(self.rules)):
            if self.rules['Type'].iloc[i] == 'rateRule':
                rule_string = self.rules['Eqn'].iloc[i]
                if '* pow' in rule_string:
                    num_pow = len(rule_string.split(sep='* pow')) - 1
                    for i in range(num_pow):
                        manip_string = rule_string.split(sep='* pow')[1]
                        base, power = manip_string.split(sep=',')
                        new_string = '* ' + base[1:] + '*m.exp(' + power[1:]
                        replace_string = '* pow' + manip_string
                        rule_string = rule_string.replace(replace_string, new_string)
                string1 += 'm.Equation(' + self.rules['ID'].iloc[i] + '.dt() == ' + rule_string + ')\n'
        string1 += '\n'
        return string1
    
    def write_model(self,file_name,time,option_dict=None):
        with open (file_name,'w+') as f:
            gekko_string = self.write_gekko(time)
            variable_string = self.write_var()
            param_string = self.write_parameters()
            intermediates_string = self.write_intermediates()
            reaction_string = self.write_equations()
            f.write(gekko_string)
            f.write(variable_string)
            f.write(param_string)
            f.write(intermediates_string)
            f.write(reaction_string)
            if option_dict is not None:
                options_string = self.write_global_options(option_dict)
                f.write(options_string)
            f.write('\nm.solve()')
            
        return
    
class CV:
    def __init__(self,value,bias=0,cost=0,critical=0,fdelay=0,fstatus=0,lower=0,lstval=1,
                 meas=1,meas_gap=1e-3,model=1,pred=1,pstatus=1,sp=0,sphi=1e20,splo=-1e20,
                 status=0,tau=60,tier=1,tr_init=0,tr_open=1,upper=1e20,vdvl=1e20,vlaction=0,
                 vlhi=1e20,vllo=-1e20,wmeas=20,wmodel=2,wsp=20,wsphi=20,wsplo=20):
        self.value = value
        self.bias = bias
        self.cost = cost
        self.critical = critical
        self.fdelay = fdelay
        self.fstatus = fstatus
        self.lower = lower
        self.lstval = lstval
        self.meas = meas
        self.meas_gap = meas_gap
        self.model = model
        self.pred = pred
        self.pstatus = pstatus
        self.sp = sp
        self.sphi = sphi
        self.splo = splo
        self.status = status
        self.tau = tau
        self.tier = tier
        self.tr_init = tr_init
        self.tr_open = tr_open
        self.upper = upper
        self.vdvl=vdvl
        self.vlaction = vlaction
        self.vlhi = vlhi
        self.vllo = vllo
        self.wmeas = wmeas
        self.wmodel = wmodel
        self.wsp = wsp
        self.wsphi = wsphi
        self.wsplo = wsplo

class MV:
    def __init__(self,value,aws=0,cost=0,critical=0,dcost=0.00001,dmax=1e20,dmaxhi=1e20,
                 dmaxlo=-1e20,dpred=0,fstatus=1,lower=-1e20,lstval=1,meas=1,mv_step_hor=0,
                 newval=1,nxtval=1,pred=1,pstatus=1,reqonctrl=0,status=1,tier=1,upper=1e20,vdvl=1e20,
                 vlaction=0,vlhi=1e20,vllo=-1e20):
        self.value = value
        self.aws = aws
        self.cost = cost
        self.critical = critical
        self.dcost = dcost
        self.dmax = dmax
        self.dmaxhi = dmaxhi
        self.dmaxlo = dmaxlo
        self.dpred = dpred
        self.fstatus = fstatus
        self.lower = lower
        self.lstval = lstval
        self.meas = meas
        self.mv_step_hor = mv_step_hor
        self.newval = newval
        self.nxtval = nxtval
        self.pred = pred
        self.pstatus = pstatus
        self.reqonctrl = reqonctrl
        self.status = status
        self.tier = tier
        self.upper = upper
        self.vdvl = vdvl
        self.vlaction = vlaction
        self.vlhi = vlhi
        self.vllo = vllo
        
class FV:
    def __init__(self, value, critical=0,dmax=1e20,dmaxhi=1e20,dmaxlo=1e20,fstatus=1,
                 lower=-1e20,lstval=1,meas=1,newval=1,pstatus=1,status=1,upper=1e20,vdvl=1e20,
                 vlaction=0,vlhi=1e20,vllo=-1e20):
        self.value = value
        self.critical = critical
        self.dmax = dmax
        self.dmaxhi = dmaxhi
        self.dmaxlo = dmaxlo
        self.fstatus = fstatus
        self.lower = lower
        self.lstval = lstval
        self.meas = meas
        self.newval = newval
        self.pstatus = pstatus
        self.status = status
        self.upper = upper
        self.vdvl = vdvl
        self.vlaction = vlaction
        self.vlhi = vlhi
        self.vllo = vllo      

class SV:
    def __init__(self, value, fstatus=0,lower=-1e20,meas=1,model=1,pred=1,upper=1e20):
        self.value = value
        self.fstatus = 0
        self.lower = lower
        self.meas = meas
        self.model = model
        self.pred = pred
        self.upper = upper
        
def read_model(filepath):
    m = bionetwork_model(filepath)
    
    return m