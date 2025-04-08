import pickle
import gzip
import copy
import joblib
from coffea.hist import Hist
import numpy as np
import json
import argparse
import os

parser = argparse.ArgumentParser(description="Scale the QCD background.")
parser.add_argument("-p", "--path", required=True, help="Path to the input files.")
args = parser.parse_args()

# Set the path variable
path = args.path

# Function to load the .pkl.gz file
def load_histograms(file_path):
    with gzip.open(file_path, 'rb') as f:
        return pickle.load(f)

# Function to save the modified histograms to a new .pkl.gz file
def save_histograms(histograms, output_path):
    with gzip.open(output_path, 'wb') as f:
        pickle.dump(histograms, f)
        f.close()
        
if os.path.exists(path + 'QCD.pkl.gz'):
  print('WARNING: QCD file already exists in path... moving to ', path + 'old/QCD.pkl.gz.old')
  if not os.path.exists(path + 'old/'):
    os.makedirs(path + 'old/')
  os.rename(path + 'QCD.pkl.gz', path + 'old/QCD.pkl.gz.old')
        

need_to_cut=True
########################################## Initialization of factors: driven from Save_QCD. To be incorporated automatized

'''
fake_rates={('e_minus','2j1b'):0.6,
	('e_plus','2j1b'):0.46,
	('m_minus','2j1b'):0.315,
	('m_plus','2j1b'):0.324,
	('e_minus','3j1b'):0.31,
	('e_plus','3j1b'):1.02,
	('m_minus','3j1b'):0.68,
	('m_plus','3j1b'):0.42,																								#FAKE RATES IN THE ht_atlas<170 CASE. NOT TO BE USED UNLESS CHECKING THAT QUESTION
	('e_minus','3j2b'):0.31,
	('e_plus','3j2b'):1.02,
	('m_minus','3j2b'):0.68,								#These are the Fake rates (parte denominador)
	('m_plus','3j2b'):0.42,
		}

cr_rates={('e_minus','2j1b'):18.77,
	('e_plus','2j1b'):28.73,
	('m_minus','2j1b'):40.67,
	('m_plus','2j1b'):59.33,
	('e_minus','3j1b'):1.84,
	('e_plus','3j1b'):0.65,
	('m_minus','3j1b'):2.56,
	('m_plus','3j1b'):3.39,								#These are the normalizations in each category (N_dat^{CR}-N_MC^CR). In yields or counts
	('e_minus','3j2b'):0.96,
	('e_plus','3j2b'):0,
	('m_minus','3j2b'):0,
	('m_plus','3j2b'):0.95,
		}
'''
# Load fake_rates and cr_rates from the JSON file
with open(path+"rates_QCD.json", "r") as f:
    rates = json.load(f)

# Convert string keys back to tuples
fake_rates = {eval(k): v for k, v in rates["fake_rates"].items()}
cr_rates = {eval(k): v for k, v in rates["cr_rates"].items()}


joint_rates={key: fake_rates[key] * cr_rates[key]/302 for key in fake_rates}   #operate both factors  AND divide by lumi. This could be done (to be more rigid in the ordering of operations) when the 'CR factors' are calculated



##############################################################################			Up to here the factors


input_file= path+"QCD_shape/QCD_shapes.pkl.gz" #point to input file (driven from QCD mc after passing the same process as other samples, only difference the relaxation of cuts) 
#input_file="prueba/QCD_prueba.pkl.gz"
output_file = path+"QCD.pkl.gz"

# Load the histograms
hists = load_histograms(input_file)

# Make a copy of the histograms to avoid modifying the original
modified_hists = copy.deepcopy(hists)


##########################################################################				First part: integrate all channels and levels to get a single shape for all the categories. 
ch=['e_fake_minus','e_fake_plus','m_fake_minus','m_fake_plus']						 # we are integrating in fake channels since those are the enriched by the selection
ch_el=['e_fake_minus','e_fake_plus']
ch_mu=['m_fake_minus','m_fake_plus']


lv=['2j1b','3j1b','3j2b']

for name, hist in modified_hists.items():
	if isinstance(hist, Hist):  # Ensure it's a coffea Hist object
		if not hist.values():
			continue          #Skip if we find not filled histograms (like dummy)
			
		print(name)	
		if name in ['MVAtth','MVAwp','pttrig','etatrig','abseeta','absmeta','counts_metl20','coste_2']:#'ept','eeta','mpt','meta',
			continue         # for the moment we skip the 'e' only or 'mu' only (lep 0 will serve for this)

		
		identifiers=hist[name].axes()[3].identifiers()
		for systematic in identifiers:
			systematic=systematic.name
			if name not in ['ept','eeta']: h0=np.zeros(len(hist.values()[('QCD', 'm_minus', '2j1b', 'norm')])) #initialize to 0 the array
			else:h0=np.zeros(len(hist.values()[('QCD', 'e_minus', '2j1b', 'norm')])) 
			for c in ch:
				for l in lv:
					key=('QCD', c, l, systematic)          #loop through channels and levels (and define the key meanwhile) to sum the entries in each category
					if key not in hist.values():
						continue                       #check if the key wasn't filled while constructing the QCD_original, skip that one and keep summing the filled ones. If this is not added, an error will raise
					
					h0+=hist.values()[key]
					
		#once this loop is finished, we have the UNIQUE (integrated over the 12 cats) SHAPE that has to be set as starting point in all the 12 categories

############################################# exception for e only and mu only		
		
			if name in ['mpt','meta']:
				for c1 in ch_el:
					for l in lv:
						key=('QCD', c1, l, systematic)
					#if key not in hist.values():
					#	continue
						if name=='mpt':h0+=modified_hists['ept'].values()[key];
						else: h0+=modified_hists['eeta'].values()[key]
					
			if name in ['ept','eeta']:
				for c1 in ch_mu:
					for l in lv:
						key=('QCD', c1, l, systematic)
					#if key not in hist.values():
					#	continue
						if name=='ept':h0+=modified_hists['mpt'].values()[key]
						else: h0+=modified_hists['meta'].values()[key]
################################################################################ end of exception for e only and mu only	
		
		#In the coming loop, we will go through the signal categories and assign to all of them the UNIQUE SHAPE just calculated
			for (channel, level), fake_rate in joint_rates.items():
 				key=('QCD', channel, level, systematic)
					
 				if key not in hist.values():
 					fill_args = {}
 					for axis in hist.axes():
						 if axis.name == "sample":
							 fill_args["sample"] = "QCD"
						 elif axis.name == "channel":
							 fill_args["channel"] = channel
						 elif axis.name == "level":
							 fill_args["level"] = level
						 elif axis.name == "syst":
							 fill_args["syst"] = systematic
						 elif axis.name == name:
							 fill_args[axis.name] = np.zeros(1)
					
 					hist.fill(**fill_args)                              #To here we are just tricking if in the original signal categories (e.g. e_plus 2j1b) the histogram wasn't filled
																    # in that case, we create the histogram for that category and fill it with dummy values
 				hist.values()['QCD',channel,level,systematic][:]=h0       #here is where we assign the Unique Shape to all the categories (both the ones which were already filled with very low values and the empty ones that now contain dummy values)
																  #VERY IMPORTANT is the [:]. If not in coffea the values of a histogram dont get updated
		

################################################################################## end of first part. At this point we have a Hist with a unique shape per varaible in all the categories


################################################################################## Second part: the scalings by factors 

for name, hist in modified_hists.items():
	if isinstance(hist, Hist):  # Ensure it's a coffea Hist object
		
		if 'channel' in hist.axes() and 'level' in hist.axes():
			for (channel, level), fake_rate in joint_rates.items():
				
				valores=hist.integrate('channel',channel).integrate('level',level).values() #Do the integration of the shape...
				if valores !={}:
					integral=np.sum(valores[('QCD','norm')]) #... if we have smth, if not...
				else:
					integral=1     #   ... a naive 1
					
				hist.scale({(channel,level):fake_rate/integral},axis=('channel','level'))  #The actual scaling by the rates and division by integral. NOTE that the integral has to be calculated before the scaling

				#if name=='MVAscore_relaxed_b10_cut':hist.values()['QCD',channel,level,'norm'][:10]=np.zeros(10);hist.values()['QCD',channel,level,'norm'][-10:]=modified_hists['MVAscore_relaxed_b10'].values()['QCD',channel,level,'norm'][-10:]

				
				

# Save the modified histograms to a new file
save_histograms(modified_hists, output_file)

print(f"Modified histograms saved to {output_file}")
