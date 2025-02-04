import sbmltoodepy.modelclasses
from scipy.integrate import odeint
import numpy as np
import operator
import math

class SBMLmodel(sbmltoodepy.modelclasses.Model):

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['parameter_1'] = sbmltoodepy.modelclasses.Parameter(1.2e-07, 'parameter_1', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T12tot"))
		self.p['parameter_2'] = sbmltoodepy.modelclasses.Parameter(2.5e-07, 'parameter_2', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T21tot"))
		self.p['parameter_3'] = sbmltoodepy.modelclasses.Parameter(2.5e-07, 'parameter_3', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("A1tot"))
		self.p['parameter_4'] = sbmltoodepy.modelclasses.Parameter(1e-06, 'parameter_4', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("dI1tot"))
		self.p['parameter_5'] = sbmltoodepy.modelclasses.Parameter(5e-07, 'parameter_5', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("A2tot"))
		self.p['parameter_6'] = sbmltoodepy.modelclasses.Parameter(1.25e-07, 'parameter_6', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPtot"))
		self.p['parameter_7'] = sbmltoodepy.modelclasses.Parameter(1.5e-08, 'parameter_7', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHtot"))
		self.p['parameter_8'] = sbmltoodepy.modelclasses.Parameter(74000.0, 'parameter_8', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kTA21"))
		self.p['parameter_9'] = sbmltoodepy.modelclasses.Parameter(14000.0, 'parameter_9', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kTA12"))
		self.p['parameter_10'] = sbmltoodepy.modelclasses.Parameter(53000.0, 'parameter_10', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kAI1"))
		self.p['parameter_11'] = sbmltoodepy.modelclasses.Parameter(24000.0, 'parameter_11', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("krAI1"))
		self.p['parameter_12'] = sbmltoodepy.modelclasses.Parameter(31000.0, 'parameter_12', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kAI2"))
		self.p['parameter_13'] = sbmltoodepy.modelclasses.Parameter(28000.0, 'parameter_13', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kTAI21"))
		self.p['parameter_14'] = sbmltoodepy.modelclasses.Parameter(140000.0, 'parameter_14', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kTAI12"))
		self.p['parameter_15'] = sbmltoodepy.modelclasses.Parameter(28000.0, 'parameter_15', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kAIrA1"))
		self.p['parameter_16'] = sbmltoodepy.modelclasses.Parameter(100000000.0, 'parameter_16', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kplus"))
		self.p['parameter_17'] = sbmltoodepy.modelclasses.Parameter(100000000.0, 'parameter_17', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kplusH"))
		self.p['parameter_18'] = sbmltoodepy.modelclasses.Parameter(6.8, 'parameter_18', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kminusON12"))
		self.p['parameter_19'] = sbmltoodepy.modelclasses.Parameter(262.0, 'parameter_19', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kminusOFF12"))
		self.p['parameter_20'] = sbmltoodepy.modelclasses.Parameter(24.7, 'parameter_20', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kminusON21"))
		self.p['parameter_21'] = sbmltoodepy.modelclasses.Parameter(267.0, 'parameter_21', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kminusOFF21"))
		self.p['parameter_22'] = sbmltoodepy.modelclasses.Parameter(7.6, 'parameter_22', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kminusH1"))
		self.p['parameter_23'] = sbmltoodepy.modelclasses.Parameter(1.6, 'parameter_23', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kminusH2"))
		self.p['parameter_24'] = sbmltoodepy.modelclasses.Parameter(0.05, 'parameter_24', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kcatON12"))
		self.p['parameter_25'] = sbmltoodepy.modelclasses.Parameter(0.002, 'parameter_25', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kcatOFF12"))
		self.p['parameter_26'] = sbmltoodepy.modelclasses.Parameter(0.08, 'parameter_26', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kcatON21"))
		self.p['parameter_27'] = sbmltoodepy.modelclasses.Parameter(0.02, 'parameter_27', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kcatOFF21"))
		self.p['parameter_28'] = sbmltoodepy.modelclasses.Parameter(0.05, 'parameter_28', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kcatH1"))
		self.p['parameter_29'] = sbmltoodepy.modelclasses.Parameter(0.24, 'parameter_29', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kcatH2"))

		self.c = {} #Dictionary of compartments
		self.c['compartment_1'] = sbmltoodepy.modelclasses.Compartment(1.0, 3, True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("compartment"))

		self.s = {} #Dictionary of chemical species
		self.s['species_1'] = sbmltoodepy.modelclasses.Species(1.2e-07, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T12"))
		self.s['species_2'] = sbmltoodepy.modelclasses.Species(5e-07, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("A2"))
		self.s['species_3'] = sbmltoodepy.modelclasses.Species(2.5e-07, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T21"))
		self.s['species_4'] = sbmltoodepy.modelclasses.Species(2.5e-07, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("A1"))
		self.s['species_5'] = sbmltoodepy.modelclasses.Species(1e-06, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("dI1"))
		self.s['species_6'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rA1"))
		self.s['species_7'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rI2"))
		self.s['species_8'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T12A2"))
		self.s['species_9'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T21A1"))
		self.s['species_10'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("A1dI1"))
		self.s['species_11'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("rA1dI1"))
		self.s['species_12'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("A2rI2"))
		self.s['species_13'] = sbmltoodepy.modelclasses.Species(1.25e-07, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAP"))
		self.s['species_14'] = sbmltoodepy.modelclasses.Species(1.5e-08, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseH"))
		self.s['species_15'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPT12A2"))
		self.s['species_16'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPT12"))
		self.s['species_17'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPT21A1"))
		self.s['species_18'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPT21"))
		self.s['species_19'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHrA1dI1"))
		self.s['species_20'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHA2rI2"))

		self.r = {} #Dictionary of reactions
		self.r['reaction_1'] = reaction_1(self)
		self.r['reaction_2'] = reaction_2(self)
		self.r['reaction_3'] = reaction_3(self)
		self.r['reaction_4'] = reaction_4(self)
		self.r['reaction_5'] = reaction_5(self)
		self.r['reaction_6'] = reaction_6(self)
		self.r['reaction_7'] = reaction_7(self)
		self.r['reaction_8'] = reaction_8(self)
		self.r['reaction_9'] = reaction_9(self)
		self.r['reaction_10'] = reaction_10(self)
		self.r['reaction_11'] = reaction_11(self)
		self.r['reaction_12'] = reaction_12(self)
		self.r['reaction_13'] = reaction_13(self)
		self.r['reaction_14'] = reaction_14(self)
		self.r['reaction_15'] = reaction_15(self)
		self.r['reaction_16'] = reaction_16(self)
		self.r['reaction_17'] = reaction_17(self)
		self.r['reaction_18'] = reaction_18(self)
		self.r['reaction_19'] = reaction_19(self)
		self.r['reaction_20'] = reaction_20(self)

		self.f = {} #Dictionary of function definitions
		self.time = 0

		self.AssignmentRules()



	def AssignmentRules(self):

		if self.time <= 0 :
			isConstantValue = self.s['species_1']._constant
			self.s['species_1']._constant = False
			self.s['species_1'].concentration = self.p['parameter_1'].value
			self.s['species_1']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_2']._constant
			self.s['species_2']._constant = False
			self.s['species_2'].concentration = self.p['parameter_5'].value
			self.s['species_2']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_3']._constant
			self.s['species_3']._constant = False
			self.s['species_3'].concentration = self.p['parameter_2'].value
			self.s['species_3']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_4']._constant
			self.s['species_4']._constant = False
			self.s['species_4'].concentration = self.p['parameter_3'].value
			self.s['species_4']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_5']._constant
			self.s['species_5']._constant = False
			self.s['species_5'].concentration = self.p['parameter_4'].value
			self.s['species_5']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_13']._constant
			self.s['species_13']._constant = False
			self.s['species_13'].concentration = self.p['parameter_6'].value
			self.s['species_13']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_14']._constant
			self.s['species_14']._constant = False
			self.s['species_14'].concentration = self.p['parameter_7'].value
			self.s['species_14']._constant = isConstantValue

		return

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['species_1'].amount, self.s['species_2'].amount, self.s['species_3'].amount, self.s['species_4'].amount, self.s['species_5'].amount, self.s['species_6'].amount, self.s['species_7'].amount, self.s['species_8'].amount, self.s['species_9'].amount, self.s['species_10'].amount, self.s['species_11'].amount, self.s['species_12'].amount, self.s['species_13'].amount, self.s['species_14'].amount, self.s['species_15'].amount, self.s['species_16'].amount, self.s['species_17'].amount, self.s['species_18'].amount, self.s['species_19'].amount, self.s['species_20'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1,0,0,0.,0,0.],[-1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1.],[ 0,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,1,0,0.,0,0.],[ 0,-1,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0.,0,0.],[ 0,0,-1,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0.,1,0.],[ 0,0,0,-1,0,0,0,-1,0,0,0,0,1,0,1,0,0,0.,0,0.],[ 0,0,0,0,-1,-1,0,0,0,0,0,0,0,1,0,1,0,0.,0,0.],[ 1,0,0,0,0,-1,0,0,-1,0,0,0,1,0,0,0,0,0.,0,0.],[ 0,1,0,0,0,0,-1,0,0,-1,0,0,0,1,0,0,0,0.,0,0.],[ 0,0,1,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0.,0,0.],[ 0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,-1,0.,0,0.],[ 0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,-1.,0,0.],[ 0,0,0,0,0,0,0,0,-1,-1,-1,-1,1,1,1,1,0,0.,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1.,1,1.],[ 0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0.,0,0.],[ 0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0.,0,0.],[ 0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0.,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0.,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0.,-1,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,-1.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['reaction_1'](), self.r['reaction_2'](), self.r['reaction_3'](), self.r['reaction_4'](), self.r['reaction_5'](), self.r['reaction_6'](), self.r['reaction_7'](), self.r['reaction_8'](), self.r['reaction_9'](), self.r['reaction_10'](), self.r['reaction_11'](), self.r['reaction_12'](), self.r['reaction_13'](), self.r['reaction_14'](), self.r['reaction_15'](), self.r['reaction_16'](), self.r['reaction_17'](), self.r['reaction_18'](), self.r['reaction_19'](), self.r['reaction_20']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['species_1'].amount, self.s['species_2'].amount, self.s['species_3'].amount, self.s['species_4'].amount, self.s['species_5'].amount, self.s['species_6'].amount, self.s['species_7'].amount, self.s['species_8'].amount, self.s['species_9'].amount, self.s['species_10'].amount, self.s['species_11'].amount, self.s['species_12'].amount, self.s['species_13'].amount, self.s['species_14'].amount, self.s['species_15'].amount, self.s['species_16'].amount, self.s['species_17'].amount, self.s['species_18'].amount, self.s['species_19'].amount, self.s['species_20'].amount], dtype = np.float64)
		self.s['species_1'].amount, self.s['species_2'].amount, self.s['species_3'].amount, self.s['species_4'].amount, self.s['species_5'].amount, self.s['species_6'].amount, self.s['species_7'].amount, self.s['species_8'].amount, self.s['species_9'].amount, self.s['species_10'].amount, self.s['species_11'].amount, self.s['species_12'].amount, self.s['species_13'].amount, self.s['species_14'].amount, self.s['species_15'].amount, self.s['species_16'].amount, self.s['species_17'].amount, self.s['species_18'].amount, self.s['species_19'].amount, self.s['species_20'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class reaction_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Activation1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_9'].value * self.parent.s['species_1'].concentration * self.parent.s['species_2'].concentration

class reaction_2:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Activation2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_8'].value * self.parent.s['species_3'].concentration * self.parent.s['species_4'].concentration

class reaction_3:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Annihilation1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_10'].value * self.parent.s['species_4'].concentration * self.parent.s['species_5'].concentration

class reaction_4:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Annihilation2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_11'].value * self.parent.s['species_6'].concentration * self.parent.s['species_5'].concentration

class reaction_5:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Annihilation3")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_12'].value * self.parent.s['species_2'].concentration * self.parent.s['species_7'].concentration

class reaction_6:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Inhibition1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_14'].value * self.parent.s['species_8'].concentration * self.parent.s['species_7'].concentration

class reaction_7:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Inhibition2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_13'].value * self.parent.s['species_9'].concentration * self.parent.s['species_5'].concentration

class reaction_8:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Release")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_15'].value * self.parent.s['species_6'].concentration * self.parent.s['species_10'].concentration

class reaction_9:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPbinding1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * (self.parent.p['parameter_16'].value * self.parent.s['species_13'].concentration * self.parent.s['species_8'].concentration - self.parent.p['parameter_18'].value * self.parent.s['species_15'].concentration)

class reaction_10:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPbinding2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * (self.parent.p['parameter_16'].value * self.parent.s['species_13'].concentration * self.parent.s['species_9'].concentration - self.parent.p['parameter_20'].value * self.parent.s['species_17'].concentration)

class reaction_11:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPbinding3")

	def __call__(self):
		return self.parent.c['compartment_1'].size * (self.parent.p['parameter_16'].value * self.parent.s['species_13'].concentration * self.parent.s['species_1'].concentration - self.parent.p['parameter_19'].value * self.parent.s['species_16'].concentration)

class reaction_12:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPbinding4")

	def __call__(self):
		return self.parent.c['compartment_1'].size * (self.parent.p['parameter_16'].value * self.parent.s['species_13'].concentration * self.parent.s['species_3'].concentration - self.parent.p['parameter_21'].value * self.parent.s['species_18'].concentration)

class reaction_13:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPcat1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_24'].value * self.parent.s['species_15'].concentration

class reaction_14:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPcat2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_26'].value * self.parent.s['species_17'].concentration

class reaction_15:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPcat3")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_25'].value * self.parent.s['species_16'].concentration

class reaction_16:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNAPcat4")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_27'].value * self.parent.s['species_18'].concentration

class reaction_17:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHbinding1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * (self.parent.p['parameter_17'].value * self.parent.s['species_14'].concentration * self.parent.s['species_11'].concentration - self.parent.p['parameter_22'].value * self.parent.s['species_19'].concentration)

class reaction_18:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHbinding2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * (self.parent.p['parameter_17'].value * self.parent.s['species_14'].concentration * self.parent.s['species_12'].concentration - self.parent.p['parameter_23'].value * self.parent.s['species_20'].concentration)

class reaction_19:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHcat1")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_28'].value * self.parent.s['species_19'].concentration

class reaction_20:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("RNaseHcat2")

	def __call__(self):
		return self.parent.c['compartment_1'].size * self.parent.p['parameter_29'].value * self.parent.s['species_20'].concentration

