#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import coffea
import numpy as np
import awkward as ak
np.seterr(divide='ignore', invalid='ignore', over='ignore')
#from coffea.arrays import Initialize # Not used and gives error
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from coffea.analysis_tools import PackedSelection

from cafea.analysis.objects import *
from cafea.analysis.selection import *

#coffea.deprecations_as_errors = True

# In the future these names will be read from the nanoAOD files

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples):
        self._samples = samples

        # Create the histograms
        # In general, histograms depend on 'sample', 'channel' (final state) and 'cut' (level of selection)
        self._accumulator = processor.dict_accumulator({
        'jetpt'  : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("pt",  "Jet p_{T} (GeV) ", 25, 0, 800)),
        'jeteta' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("eta", "Jet eta", 35, -4.7, 4.7)), #antes 25, -2.5, 2.5
        'jetpteta' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("pt",  "Jet p_{T} (GeV) ", [25, 30, 60, 120,121]), hist.Bin("abseta", "Jet eta", [0, 1, 1.8,2.4,4.7])),  #antes solo hast 2.4
        'jetptetaflav' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Bin("pt",  "Jet p_{T} (GeV) ", [25, 30, 60, 120,121]), hist.Bin("abseta", "Jet eta", [0, 1, 1.8,2.4,2.5]), hist.Bin("flav", "Flavor", [0, 4, 5]) ),
        'jetptetaflav' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Bin("pt",  "Jet p_{T} (GeV) ", [25, 30, 60, 120]), hist.Bin("abseta", "Jet eta", [0, 1, 1.8,2.5]), hist.Bin("flav", "Flavor", [0, 4, 5]) ),
        'pileup'  : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("pu",  "N_{vtx}", 16, 0, 80)),

        })

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    # Main function: run on a given dataset
    def process(self, events):
        # Dataset parameters
        dataset = events.metadata['dataset']
        year   = self._samples[dataset]['year']
        xsec   = self._samples[dataset]['xsec']
        sow    = self._samples[dataset]['nSumOfWeights' ]
        isData = self._samples[dataset]['isData']
        datasets = ['SingleMuon', 'SingleElectron', 'EGamma', 'MuonEG', 'DoubleMuon', 'DoubleElectron']
        for d in datasets: 
          if d in dataset: dataset = dataset.split('_')[0] 

        # Initialize objects
        met = events.MET
        e   = events.Electron
        mu  = events.Muon
        j   = events.Jet
        
        # Call to pileup
        PV=events.PV
        pileup=PV.npvsGood #this gives pile up per event
        lengths=ak.num(j)
        expanded = np.hstack([np.full(n, v) for v, n in zip(pileup, lengths)])
        pileup_mimic_jets = ak.Array(ak.layout.ListOffsetArray64(ak.layout.Index64(np.cumsum([0] + lengths.tolist())), ak.layout.NumpyArray(expanded))) #Coding awkward arrays such that 
                                                                                                                                                        #we get pu repeated per njets in an event
        ###### end pu

        e["btagDeepB"] = ak.fill_none(e.matched_jet.btagDeepB, -99)
        mu["btagDeepB"] = ak.fill_none(mu.matched_jet.btagDeepB, -99)

        # Muon selection
        #mu['isPres'] = isPresMuon(mu.dxy, mu.dz, mu.sip3d, mu.looseId)
        #mu['isTight']= isTightMuon(mu.pt, mu.eta, mu.dxy, mu.dz, mu.pfRelIso03_all, mu.sip3d, mu.mvaTTH, mu.mediumPromptId, mu.tightCharge, mu.looseId, minpt=10)
        #mu['isGood'] = mu['isPres'] & mu['isTight']

        mu["isLoose"] = MuonLoose(mu.pt, mu.eta, mu.dxy, mu.dz, mu.sip3d, mu.mediumPromptId, mu.btagDeepB, ptCut=20, etaCut=2.4)
        mu["isMVA"]= MuonMVA(mu.miniPFRelIso_all, mu.mvaTTH)
        mu['isGood'] = mu['isLoose'] & mu['isMVA']

        leading_mu = mu[ak.argmax(mu.pt,axis=-1,keepdims=True)]
        leading_mu = leading_mu[leading_mu.isGood]
        
        mu = mu[mu.isGood]

        # Electron selection
        #e['isPres']  = isPresElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, e.lostHits, minpt=15)
        #e['isTight'] = isTightElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, e.mvaTTH, e.mvaFall17V2Iso, e.lostHits, e.convVeto, e.tightCharge, e.sieie, e.hoe, e.eInvMinusPInv, minpt=15)
        #e['isClean'] = isClean(e, mu, drmin=0.05)
        #e['isGood']  = e['isPres'] & e['isTight'] & e['isClean']

        e['isLoose'] = ElecLoose(e.pt, e.eta, e.lostHits, e.sip3d, e.dxy, e.dz, e.btagDeepB, e.convVeto, e.mvaFall17V2noIso_WPL, 20, 2.4)
        e['isMVA']   = ElecMVA(e.miniPFRelIso_all, e.mvaTTH)
        e['isGood'] = e['isLoose'] & e['isMVA']

        leading_e = e[ak.argmax(e.pt,axis=-1,keepdims=True)]
        leading_e = leading_e[leading_e.isGood]

        e  =  e[e .isGood]

        nElec = ak.num(e)
        nMuon = ak.num(mu)

        twoLeps   = (nElec+nMuon) == 2
        threeLeps = (nElec+nMuon) == 3
        twoElec   = (nElec == 2)
        twoMuon   = (nMuon == 2)
        e0 = e[ak.argmax(e.pt,axis=-1,keepdims=True)]
        m0 = mu[ak.argmax(mu.pt,axis=-1,keepdims=True)]

        # Jet selection

        jetptname = 'pt_nom' if hasattr(j, 'pt_nom') else 'pt'

        j['isClean'] = isClean(j, e, drmin=0.4)& isClean(j, mu, drmin=0.4)
        j['isGood']  = isTightJet(getattr(j, jetptname), j.eta, j.jetId, jetPtCut = 25)
        goodJets = j[(j.isClean)&(j.isGood)]
        njets = ak.num(goodJets)
        goodPileup=pileup_mimic_jets[(j.isClean)&(j.isGood)]
        ht = ak.sum(goodJets.pt,axis=-1)
        j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]

        ### We need weights for: normalization, lepSF, triggerSF, pileup, btagSF...
        weights = coffea.analysis_tools.Weights(len(events))
        weights.add('norm', np.ones_like(met.pt))

        eftweights = events['EFTfitCoefficients'] if hasattr(events, "EFTfitCoefficients") else []

        hout = self.accumulator.identity()
        normweights = weights.weight().flatten() 
        #hout['SumOfEFTweights'].fill(eftweights, sample=dataset, SumOfEFTweights=varnames['counts'], weight=normweights)

        flavSelection = {'b': (np.abs(goodJets.hadronFlavour) == 5), 'c': (np.abs(goodJets.hadronFlavour) == 4), 'l': (np.abs(goodJets.hadronFlavour) <= 3) }

        #WP = {'all' : -999., 'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100}
        #WP = {'all' : -999., 'loose': 0.0521, 'medium': 0.3033, 'tight': 0.7489}
        WP = {'all' : -999, 'loose': 0.1522 , 'medium': 0.4941, 'tight': 0.8001 } # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        btagSelection = {}
        for wp, wpvals in WP.items():
          btagSelection[wp] = (goodJets.btagDeepB>wpvals) 

        for jetype in ['b', 'c', 'l']:
          for wp in WP.keys():
            mask = (flavSelection[jetype])&(btagSelection[wp])
            selectjets = goodJets[mask]
            
            #Pileup tricks           
            selectpileup=ak.mask(goodPileup,mask)
            pileup=ak.flatten(selectpileup)
            nones_mask=~ak.is_none(pileup)
            pileup=pileup[nones_mask]
            #end of pileup tricks
            
            pts     = ak.flatten(selectjets.pt)
            etas    = ak.flatten(selectjets.eta)
            absetas = ak.flatten(np.abs(selectjets.eta))
            flavarray = np.zeros_like(pts) if jetype == 'l' else (np.ones_like(pts)*(4 if jetype=='c' else 5))
            weights =  np.ones_like(pts)
            hout['jetpt'].fill(WP=wp, Flav=jetype,  pt=pts, weight=weights)
            hout['jeteta'].fill(WP=wp, Flav=jetype,  eta=etas, weight=weights)
            hout['jetpteta'].fill(WP=wp, Flav=jetype,  pt=pts, abseta=absetas, weight=weights)
            hout['jetptetaflav'].fill(WP=wp, pt=pts, abseta=absetas, flav=flavarray, weight=weights)
            hout['pileup'].fill(WP=wp,Flav=jetype,pu=pileup,weight=weights)
    
        return hout


    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    # Load the .coffea files
    outpath= './coffeaFiles/'
    samples     = load(outpath+'samples.coffea')
    topprocessor = AnalysisProcessor(samples)

