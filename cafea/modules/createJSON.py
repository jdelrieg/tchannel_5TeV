'''
  createJSON.py

  This script looks for samples an create a json with paths and metadata
  You can use the json file to run on a dataset, providing a prefix (rxd redirector, local path...)

  You can execute this script in different modes:
  1) Rootfiles locally or through xrootd with the format:
       Directory: [Prefix]/path/to/files/
       Names:     sampleName1_0.root sampleName1_1.root ... sampleName2_0.root sampleName2_1.root ...
     (a single dir containing files for multiple datasets.. you can pass a list of sample names)
  2) Rootfiles in a local dir or accesible through xrootd with the format:
       [Prefix]/path/to/files/
       subdir1/ subdir2/...
       tree_0.root tree_1.root...
     (similar to the output from crab, where you have a structure of folders and want to collect all the rootfiles inside a parent folder)
  3) Dataset published in DAS
     
  Inputs:
  - year
  - cross section (or file to read the cross section + name in file)

  - tree name ('Events' by default)

  Example Usage:
    1) and 2)
    MC:
    >> python createJSON.py [path] --prefix "root:/REDIRECTOR/" --sampleName "ttllnunu0, ttllnunu1" --xsec cafea/cfg/xsec.cfg --xsecName TTTo2L2Nu  --year 2018
    Data:
    >> python createJSON.py [path] --prefix "root:/REDIRECTOR/" --sampleName MuonEG_2018 --year 2018

    3)
    MC:
    >> python createJSON.py [DAS_dataset] --DAS --sampleName TTTo2L2Nu --xsec cafea/cfg/xsec.cfg --prefix root://your.favorite.redirector/
    Data:
    >> python createJSON.py [DAS_dataset] --DAS --sampleName DoubleMuon_2017 --year 2017

  Note: the "--xsecName TTTo2L2Nu" argument is only needed if sampleName does not exist in cafea/cfg/xsec.cfg


'''

import os, sys
from coffea.util import save
from cafea.modules.DASsearch import GetDatasetFromDAS, RunDasGoClientCommand
from cafea.modules.paths import cafea_path
from cafea.modules.fileReader import GetFiles, GetAllInfoFromFile, GetListOfWCs, GetFilesFromJson, GetSumWeightsFromRunsTree
from cafea.modules.samples import loadxsecdic
import argparse
import json

def main():
  parser = argparse.ArgumentParser(description='Create json file with list of samples and metadata')
  parser.add_argument('path'              , default=''           , help = 'Path to directory or DAS dataset')
  parser.add_argument('--prefix','-p'     , default=''           , help = 'Prefix to add to the path (e.g. redirector)')
  parser.add_argument('--sampleName','-s' , default=''           , help = 'Sample name, used to find files and/or output name')
  parser.add_argument('--xsec','-x'       , default=1            , help = 'Cross section (number or file to read)')
  parser.add_argument('--xsecName'        , default=''           , help = 'Name in cross section .cfg (only if different from sampleName)')
  parser.add_argument('--year','-y'       , default=-1           , help = 'Year')
  parser.add_argument('--treename'        , default='Events'     , help = 'Name of the tree')
  parser.add_argument('--histAxisName'    , default=''           , help = 'Name for the samples axis of the coffea hist')
  parser.add_argument('--update', '-u'    , action='store_true'  , help = 'Name for the samples axis of the coffea hist')
  parser.add_argument('--split'           , default=0            , help = 'Split a json file into different files')

  parser.add_argument('--DAS'             , action='store_true'  , help = 'Search files from DAS dataset')
  parser.add_argument('--nFiles'          , default=None         , help = 'Number of max files (for the moment, only applies for DAS)')

  parser.add_argument('--outname','-o'    , default=''           , help = 'Out name of the json file')
  parser.add_argument('--options'         , default=''           , help = 'Sample-dependent options to pass to your analysis')
  parser.add_argument('--verbose','-v'    , action='store_true'  , help = 'Activate the verbosing')

  args, unknown = parser.parse_known_args()
  #cfgfile     = args.cfgfile
  path         = args.path
  prefix       = args.prefix
  sample       = args.sampleName
  xsec         = args.xsec
  xsecName     = args.xsecName
  year         = args.year
  options      = args.options
  treeName     = args.treename
  histAxisName = args.histAxisName
  outname      = args.outname
  isDAS        = args.DAS
  nFiles       = int(args.nFiles) if not args.nFiles is None else None
  verbose      = args.verbose
  update       = args.update
  split        = int(args.split)
  if histAxisName == '': histAxisName = args.sampleName
  if update:
    AddSumWeightToJson(path)
    return

  elif split > 0:
    SplitJson(path, split)
    return

  # Get the xsec for the dataset
  if xsecName == '': xsecName = sample
  try:
    xsec = float(xsec)
  except:
    xsecdic = loadxsecdic(xsec, verbose)
    if xsecName in xsecdic.keys():
      xsec = xsecdic[xsecName]
    else:
      xsec = 1

  sampdic = {}
  sampdic['xsec']         = xsec
  sampdic['year']         = year
  sampdic['treeName']     = treeName
  sampdic['histAxisName'] = histAxisName
  sampdic['options']      = options

  # 1) Search files with name 'sample' or 'sample1, sample2...' in path
  if not isDAS:
    paths = [path] if not ',' in path else path.replace(' ', '').split(',')
    files = []
    for p in paths:
      filesWithPrefix = GetFiles(prefix+p, sample)
      files += [(f[len(prefix):]) for f in filesWithPrefix]
      
  # 2) Get all rootfiles in a dir and all the sub dirs
    if filesWithPrefix == []:
      if not ',' in path:
        filesWithPrefix = GetFiles(prefix+path, '')
        files = [(f[len(prefix):]) for f in filesWithPrefix]
        if nFiles is not None: 
           files = files[:nFiles]
           filesWithPrefix = filesWithPrefix[:nFiles]
      else:
        paths = path.replace(' ', '').split(',')
        files = []
        for p in paths:
          filesWithPrefix = GetFiles(prefix+p, '')
          files += [(f[len(prefix):]) for f in filesWithPrefix]

  # 3) Search files in DAS dataset
  #   NOTE: For DAS searches, the isData flag is determined from DAS itself, not the files
  else:
    dataset = path
    dicFiles = GetDatasetFromDAS(dataset, nFiles, options='file', withRedirector=prefix)
    files = [f[len(prefix):] for f in dicFiles['files']]
    filesWithPrefix = dicFiles['files']

    # This DAS command for some reason returns the output doubled and will look something like this:
    #   output = " \ndata  \ndata  \n "
    # So we strip off the whitespace and spurious newlines and then only take the first of the duplicates
    dataset_part = "{name} | grep dataset.datatype".format(name=path)
    output = RunDasGoClientCommand(dataset=dataset_part,mode='')
    cleaned_output = output.strip().replace('\n',' ').split()[0]
    if cleaned_output == 'data':
      isData = True
    elif cleaned_output == 'mc':
      isData = False
    else:
      raise RuntimeError("Unknown datatype returned by DAS: ---{}---".format(output))

  # When getting data from DAS, we don't need to query every single file to get the number of events
  if isDAS and isData and nFiles is None:
    dataset_part = "{name} | sum(file.nevents)".format(name=path)
    output = RunDasGoClientCommand(dataset=dataset_part,mode='file,')
    output = float(output.split(':')[-1].strip())

    # For data this this should all be the same
    nEvents,nGenEvents,nSumOfWeights = output,output,output
  else:
    nEvents, nGenEvents, nSumOfWeights, isData = GetAllInfoFromFile(filesWithPrefix, treeName)

  #sampdic['WCnames'] = GetListOfWCs(filesWithPrefix[0])
  sampdic['files']         = files
  sampdic['nEvents']       = nEvents
  sampdic['nGenEvents']    = nGenEvents
  sampdic['nSumOfWeights'] = nSumOfWeights
  sampdic['isData']        = isData

  if outname == '':
    outname = sample
    if   isinstance(outname, list): outname = outname[0]
    elif ',' in outname:outname = sample.replace(' ', '').split(',')[0]
  if not outname.endswith('.json'): outname += '.json'
  with open(outname, 'w') as outfile:
    json.dump(sampdic, outfile, indent=2)
    print('>> New json file: %s'%outname)


def AddSumWeightToJson(pathJson):
  if isinstance(pathJson, str) and os.path.isdir(pathJson):
    folder = os.listdir(pathJson)
    jsonfiles = []
    for j in folder:
      if j.endswith('.json'): jsonfiles.append(j)
    return AddSumWeightToJson(jsonfiles)
  elif isinstance(pathJson, str) and ',' in pathJson:
    pathJson = pathJson.replace(' ', '').split(',')
    return AddSumWeightToJson(pathJson)
  elif isinstance(pathJson, list):
    for j in pathJson:
      AddSumWeightToJson(j)
    return
  else: # Ok, this is a json file
    files = GetFilesFromJson(pathJson)
    with open(pathJson) as j: 
      dic = json.load(j)
    if 'sumScaleWeights' in dic and 'sumPDFWeights' in dic: return
    weightsPDF, weightsScale = GetSumWeightsFromRunsTree(files)
    with open(pathJson, 'w') as j:
      dic['sumScaleWeights'] = list(weightsScale)
      dic['sumPDFWeights'] = list(weightsPDF)
      json.dump(dic, j, indent=2)
    print('>> Updated json: ', pathJson)

import numpy as np
def SplitJson(pathJson, nfiles):
  if isinstance(pathJson, str) and os.path.isdir(pathJson):
    folder = os.listdir(pathJson)
    jsonfiles = []
    for j in folder:
      if j.endswith('.json'): jsonfiles.append(j)
    return SplitJson(jsonfiles, nfiles)
  elif isinstance(pathJson, str) and ',' in pathJson:
    pathJson = pathJson.replace(' ', '').split(',')
    return SplitJson(pathJson, nfiles)
  elif isinstance(pathJson, list):
    for j in pathJson:
      AddSumWeightToJson(j, nfiles)
    return
  else: # Ok, this is a json file
    files = GetFilesFromJson(pathJson)
    print('Number of total files: ', len(files))
    chunks = np.array_split(files, nfiles)
    newjson = pathJson[:pathJson.find('.json')+1]
    if newjson.endswith('.'): newjson = newjson[:-1]
    with open(pathJson) as j: 
      dic = json.load(j)
    for i, chunk in enumerate(chunks):
      with open(newjson + '_part%i.json'%i, 'w') as j:
        dic['files'] = chunk.tolist()
        json.dump(dic, j, indent=2)
        print('Created json: ', (newjson + '_part%i.json'%i), ' with a total of %i files'%(len(chunk.tolist())))

if __name__ == '__main__':
  main()


