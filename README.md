## Analysis

The main analysis is in `tchannel5TeV_charge.py`. To run it, use the `run.py` scrip. You need to have json files with the samples. For example, send a job to process tt sample with 

    python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part0.json -n 64 -j -o TTPS

To execute all the samples, you can use the script `run5TeV_tchannel.sh`.

## The config script

The config script, `analysis/tt5TeV/config.py` is imported by all the plotting scripts of the repo. It contains basic information about paths, rebinnings and other stuff. You might want to modify some of the paths, especially some of the output paths, which are by default created according to the current date. It also define some labels to use in the legends, etc.

## QCD estimate

To estimate QCD, you need to have run over all the samples and have a folder with all the .pkl.gz files. The script `SaveQCD.py` takes those inputs and creates a `QCD.pkl.gz` file with the QCD estimate. Run the script as:

    python analysis/tt5TeV/SaveQCD.py -p histos5TeV/16jan2023/ -n 32

You can draw QCD plots (that is, events with fake electrons and muons) using the `DrawQCD.py` script. For example:
    
    python analysis/tt5TeV/DrawQCD.py -p histos5TeV/16jan2023/

## Plotting and tables

There are several scripts to create plots and tables. You can find a description of some of them below.

 - ControlPlots.py: Produce control plots. You can modify the variables, channels and levels. This script produces data/MC plots in general. It is also used to produce MVA and QCD control plot.
Example:

    python analysis/tt5TeV/ControlPlots.py -p histos5TeV/16jan2023/ -n 16

 <!-- - PlotSystematics.py: Produce systematic plots, including comparisions. By default, it is done for ttbar only.
Example:

    python analysis/tt5TeV/PlotSystematics.py -p histos5TeV/16jan2023/ 

 - DrawTTmodAltSamp.py: Produce plots for hdamp and UE tune uncertainties from alternative samples where the systematic uncertainty of the alternative predictions is shown, as a function of jet and b-tag multiplicities.
Example:

    python analysis/tt5TeV/DrawTTmodAltSamp.py -p histos5TeV/22jan2023/
-->
## MVA plots

MVA plots are produced with several scripts. Some of them are produced with the `ControlPlots.py` script, mentioned above. Training and performance plots are produced with `MVAplots.py` scripts as follows:

    python analysis/tt5TeV/MVAplots.py

In this case, the relevant paths are directly hardcoded in the script.

## Master Plot

The final pickle file and the master plots (yields and shapes) are produced with `MasterPlot_charge.py`script. This has to be run before creating the datacards.

    python analysis/tt5TeV/MasterPlot_charge.py -p histos5TeV/22jan2023/ -n 64

The table yields are calculated with:

    python analysis/tt5TeV/PrintYields_charge.py -p histos5TeV/22jan2023/

## Datacards

You need to create rootfiles and datacards to produce fits and extract the cross section.
Rootfiles are created with `analysis/tt5TeV/SaveRootfile_conlowess.py` and then, datacards are created with `analysis/tt5TeV/CreateDatacard_v2.py`.
A script to create all the needed rootfiles and datacards in the analysis is executed as follows:

    python analysis/tt5TeV/RunCombineProducer_v2.py -p histos5TeV/16jan2023/ -o carpetasalida

The uncertainties must be specified in the python `analysis/tt5TeV/CreateDatacard_v2.py` script. The distributions are specified first in the  `analysis/tt5TeV/SaveRootfile_conlowess.py` script. The distributions creator file has the 'surname' `conlowess` since the lowess method is applied by default. Alternatively in that script, one can choose not to apply the lowess smoothing method. In that script, also the QCD shape uncertainties are calculated.

## Fit


This is divided in several steps and is done using CMSSW software. 

### Combined card, workspace and fit result

Script `fitsHelper_tchannel.py` (combined card). The following options are needed:

    -i <path with cards> -ch <channels (typically we will use 'e_minus,e_plus,m_minus,m_plus')> -r <regions (typically '2j1b,3j1b,3j2b')> //-uD for the first, -O for the two latter if we want observed results

Additionally, one can pass:

    -V (for verbosity), -j n (number of cores to be used by the machine in which the script will be launched)

For the workspace bulding, we may use the `text2workspace.py` script from `Combine`in the following way

    text2workspace.py combcard_2j1b3j1b3j2b.txt -o workspace.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=.*/t.*chan.*:r_tch[1, 0, 6]'

If we want a single POI for t channel and t-bar channel contributions. If one wants two different POIS targetting the two individual contributions, the POI mapping (last input) sould be modified into `--PO 'map=.*/tchan:r_tch[1, 0, 6]' --PO "map=.*/tbarchan:r_tchbar[1, 0, 6]" `

To obtain the fit result, one should use the `FitDiagnostics` (in the single POI case) or the `MultiDimFit` (in the multiple POI case) as follows:
    
    combine -M FitDiagnostics  --expectSignal 1 workspace.root  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1 &> fitOutput_2j1b3j1b3j2b.txt
or

    combine -M MultiDimFit workspace.root --setParameters r_tch=1,r_tchbar=1  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --algo singles  -t -1 &> fitOutput_2j1b3j1b3j2b.txt
Remember to quit -t -1 if you want the unblinded results.

### Impact plots

Depending if one is in the single POI or multiple POI case, there are 4 commands to do:

    combineTool.py -M Impacts -m 125 -d workspace_join.root --doInitialFit --setParameters r_tch=1  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1  &> fitOutput_join_2j1b3j1b3j2b.txt
    combineTool.py -M Impacts -m 125 -d ../workspace_join.root --doFits --setParameters r_tch=1  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1  --parallel 60 &> ../fitOutput_join_2j1b3j1b3j2b.txt
    combineTool.py -M Impacts -m 125 -d ../workspace_join.root  -o impacts.json --setParameters r_tch=1  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1  --parallel 60 &> ../fitOutput_join_2j1b3j1b3j2b.txt
    combineTool.py -M Impacts -m 125 -d ../workspace_join.root  -o impacts.json --setParameters r_tch=1  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1  --parallel 60 &> ../fitOutput_join_2j1b3j1b3j2b.txt
(single POI case), or:
    
    combineTool.py -M Impacts -m 125 -d workspace.root --doInitialFit --setParameters r_tchplus=1,r_tchmninus=1 --redefineSignalPOI r_tch,r_tchbar --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 -t -1 &>fitOutput_2j1b3j1b3j2b.txt
    combineTool.py -M Impacts -m 125 -d ../workspace.root --doFits --setParameters r_tch=1,r_tchbar=1 --redefineSignalPOI r_tch,r_tchbar  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1 --parallel 60 &> ../fitOutput_2j1b3j1b3j2b.txt
    combineTool.py -M Impacts -m 125 -d ../workspace.root -o impacts.json --setParameters r_tch=1,r_tchbar=1 --redefineSignalPOI r_tch,r_tchbar  --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000   -t -1 --parallel 60 &> ../fitOutput_2j1b3j1b3j2b.txt
    plotImpacts.py -i impacts.json -o impacts_rtbar --POI r_tchbar 
(multiple POI case, the 4th commmand may be modified to obtain the plot for the impacts of the other POI)

Again, remove -t -1 if the unblinded results are the desired ones. Also, `--parallel N` controls the number of cores one wants to run over.
Finally, the customImpacts.py script can be used to produce nicer plots.


    
### Scan plots
Specifically, for the scan plots `-g -P 250` is needed
