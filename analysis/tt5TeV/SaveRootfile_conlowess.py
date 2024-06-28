# python ljets/test.py --path histos/2jun2022_btagM_jet25/
import copy
from config import *
from cafea.modules.fileReader import *

import statsmodels.api as sm
import argparse

def get_smoothed_scale_factor( smoothed_ratio_diff, nom_hist, var_hist, nom_var, var_var):
    '''determine an overall systematic template smoothing which minimizes the chi-squared between
        the smoothed systematic and the unsmoothed template it is derived from.
        Taken from:
             https://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2018_077_v4.pdf'''

    total_var = np.sqrt( nom_var + var_var)
    unsmoothed_diff = var_hist - nom_hist
    ratio_factor = smoothed_ratio_diff * nom_hist
    numerator = np.sum( ratio_factor * unsmoothed_diff / total_var )
    denominator = np.sum( (ratio_factor / np.sqrt(total_var))**2 )
    return numerator/denominator


if var is None: var = 'absu0eta'

plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
name = GetFileNameFromPath(path)

htop  = GetHisto(path+ 'tchannel_analysis' +   '.pkl.gz', var, group=None)
htbar=GetHisto(path+ 'tbarchannel_analysis' +   '.pkl.gz', var, group=None) #These two lines load the nominal



tuneup_top , tunedo_top = GetModSystHistos(path, 'tchannel_UE',    'UE', var=var)
tuneup_tbar , tunedo_tbar = GetModSystHistos(path, 'tbarchannel_UE',    'UE', var=var)
tuneup=tuneup_top+tuneup_tbar
tunedo=tunedo_top+tunedo_tbar #These 4 lines load the classical UE histos


tunehistoup=copy.deepcopy(tuneup)
data_dict_up=tunehistoup.values()
tunehistodo=copy.deepcopy(tunedo)
data_dict_do=tunehistodo.values() #We make copies and build an object type that allow us making changes and saving them

#Lowess starting for UE
#xvalLowess = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65,0.75,0.85]) if var=='MVAscore_pruned' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
xvalLowess = np.array([0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425,0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875,0.925, 0.975]) if var=='MVAscore_pruned' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 

yvalUpLowess = data_dict_up[('tchan',ch[0],level,'UEUp')]
yvalDnLowess = data_dict_do[('tchan',ch[0],level,'UEDown')]
nominalLowess = htop.values()[('tchannel', ch[0], level, 'norm')]+htbar.values()[('tbarchannel', ch[0], level, 'norm')]

ratio_diff = (yvalUpLowess/nominalLowess - yvalDnLowess/nominalLowess)/2.
ratio_diff = np.nan_to_num(ratio_diff, nan=0)
diff_smooth = sm.nonparametric.lowess(ratio_diff, xvalLowess)[:,1]

# Chi2 part
up_scale = get_smoothed_scale_factor( diff_smooth, nominalLowess, yvalUpLowess, np.var(nominalLowess), np.var(yvalUpLowess))
down_scale = get_smoothed_scale_factor( diff_smooth, nominalLowess, yvalDnLowess, np.var(nominalLowess), np.var(yvalDnLowess))

up_ratio = (1 + up_scale*diff_smooth)
down_ratio = ( 1 + down_scale*diff_smooth)



tunehistoup.values()[('tchan',ch[0],level,'UEUp')][()]= nominalLowess * np.nan_to_num(up_ratio, nan=1.0)
tunehistodo.values()[('tchan',ch[0],level,'UEDown')][()] = nominalLowess * np.nan_to_num(down_ratio, nan=1.0)	#These two lines assign the lowess-returned values (right hand side of eq) to the copies of the histos we made at the beginning (left hand side of eq)
#Lowess up to here for UE

del(xvalLowess,yvalUpLowess,yvalDnLowess,nominalLowess,ratio_diff,diff_smooth,up_scale,down_scale,up_ratio,down_ratio)

#And now we repeat for hdamp
hdampup_top , hdampdo_top = GetModSystHistos(path, 'tchannel_hdamp',    'hdamp', var=var)
hdampup_tbar , hdampdo_tbar = GetModSystHistos(path, 'tbarchannel_hdamp',    'hdamp', var=var)
hdampup=hdampup_top+hdampup_tbar
hdampdo=hdampdo_top+hdampdo_tbar

hdamphistoup=copy.deepcopy(hdampup)
data_dict_up_hdamp=hdamphistoup.values()
hdamphistodo=copy.deepcopy(hdampdo)
data_dict_do_hdamp=hdamphistodo.values() #We make copies and build an object type that allow us making changes and saving them

#Lowess starting for UE
#xvalLowess = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65,0.75,0.85]) if var=='MVAscore_pruned' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
xvalLowess = np.array([0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425,0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875,0.925, 0.975]) if var=='MVAscore_pruned' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
yvalUpLowess = data_dict_up_hdamp[('tchan',ch[0],level,'hdampUp')]
yvalDnLowess = data_dict_do_hdamp[('tchan',ch[0],level,'hdampDown')]
nominalLowess = htop.values()[('tchannel', ch[0], level, 'norm')]+htbar.values()[('tbarchannel', ch[0], level, 'norm')]

ratio_diff = (yvalUpLowess/nominalLowess - yvalDnLowess/nominalLowess)/2.
ratio_diff = np.nan_to_num(ratio_diff, nan=0)
diff_smooth = sm.nonparametric.lowess(ratio_diff, xvalLowess)[:,1]

# Chi2 part
up_scale = get_smoothed_scale_factor( diff_smooth, nominalLowess, yvalUpLowess, np.var(nominalLowess), np.var(yvalUpLowess))
down_scale = get_smoothed_scale_factor( diff_smooth, nominalLowess, yvalDnLowess, np.var(nominalLowess), np.var(yvalDnLowess))

up_ratio = (1 + up_scale*diff_smooth)
down_ratio = ( 1 + down_scale*diff_smooth)



hdamphistoup.values()[('tchan',ch[0],level,'hdampUp')][()]= nominalLowess * np.nan_to_num(up_ratio, nan=1.0)
hdamphistodo.values()[('tchan',ch[0],level,'hdampDown')][()] = nominalLowess * np.nan_to_num(down_ratio, nan=1.0)	#These two lines assign the lowess-returned values (right hand side of eq) to the copies of the histos we made at the beginning (left hand side of eq)
#Lowess up to here for UE

del(xvalLowess,yvalUpLowess,yvalDnLowess,nominalLowess,ratio_diff,diff_smooth,up_scale,down_scale,up_ratio,down_ratio)

#plt.AddExtraBkgHist([hdamphistoup,hdamphistodo,tunehistoup,tunehistodo],add=True) #If lowess
plt.AddExtraBkgHist([hdampup,hdampdo,tuneup,tunedo],add=True) #If no lowess

plt.SetOutpath(path + '/combineFiles/')
plt.SetVerbose(verbose)
if not doData:
  plt.SetDataName('Asimov')
RebinVar(plt, var)

def SaveFile(level, channel):
  categories = {'level':level, 'channel':channel}
  lev  = categories['level'  ] if not isinstance(categories['level'  ], list) else categories['level'  ][0];
  chan = categories['channel'] if not isinstance(categories['channel'], list) else categories['channel'][0];
  outname = '%s_%s_%s'%(var, chan, lev)
  plt.SaveCombine(var, outname, categories=categories)

if __name__=="__main__":
  if not isinstance(ch, list): ch = [ch]
  if not isinstance(level, list): level = [level]
  for l in level:
    for c in ch:
      SaveFile(l, c)
