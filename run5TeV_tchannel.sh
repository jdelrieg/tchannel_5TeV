### MC

python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part0.json -n 32 -j   -s 10000  -o TTPS_part0 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part1.json -n 32 -j   -s 10000  -o TTPS_part1 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part2.json -n 32 -j   -s 10000  -o TTPS_part2 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part3.json -n 32 -j  -s 10000   -o TTPS_part3 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part4.json -n 32 -j   -s 10000  -o TTPS_part4 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part5.json -n 32 -j  -s 10000   -o TTPS_part5 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part6.json -n 32 -j   -s 10000  -o TTPS_part6 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TTPS_part7.json -n 32 -j   -s 10000  -o TTPS_part7 -p mva6_nosyst/

python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/DYJetsToLL_M_10to50.json -n 32 -j   -s 10000  -o DY_M10to50  -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/DYJetsToLL_MLL50.json -n 32 -j   -s 10000  -o DY_M50 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/WJetsToLNu.json -n 32 -j     -o WJetsToLNu -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/W0JetsToLNu.json -n 32 -j  -s 10000 -o W0JetsToLNu -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/W1JetsToLNu.json -n 32 -j  -s 10000 -o W1JetsToLNu -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/W2JetsToLNu.json -n 32 -j  -s 10000 -o W2JetsToLNu -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/W3JetsToLNu.json -n 32 -j  -s 10000 -o W3JetsToLNu -p mva6_nosyst/


python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/tW_noFullHad_part0.json -n 32 -j  -s 10000   -o tW_part0 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/tW_noFullHad_part1.json -n 32 -j   -s 10000  -o tW_part1 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/tW_noFullHad_part2.json -n 32 -j  -s 10000   -o tW_part2 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/tW_noFullHad_part3.json -n 32 -j   -s 10000  -o tW_part3 -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/tbarW_noFullHad.json -n 32 -j   -s 10000  -o tbarW -p mva6_nosyst/



python analysis/tt5TeV/run.py cafea/json/5TeV/mva_split/tbarchannel_analysis_pdfs.json -n 32 -j -s 10000 -o tbarchannel_analysis -p mva6_nosyst/ #these for the 70% part we take for analysis
python analysis/tt5TeV/run.py cafea/json/5TeV/mva_split/tchannel_analysis_pdfs.json -n 32 -j -s 10000 -o tchannel_analysis -p mva6_nosyst/ 

### skimmed data
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/HighEGJet_skim.json -n 32 -j   -s 10000  -o HighEGJet -p mva6_nosyst/
python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/SingleMuon_skim.json -n 32 -j   -s 10000  -o SingleMuon -p mva6_nosyst/


### Uncertainties (these have to be produced for tchannel, we dont have them yet)
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TT_hdampUp.json -n 32 -j  -s 10000   -o TT_hdampUp -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TT_hdampDown.json -n 32 -j   -s 10000  -o TT_hdampDown -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TT_UEUp.json -n 32 -j   -s 10000  -o TT_UEUp -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TT_UEDown.json -n 32 -j   -s 10000  -o TT_UEDown -p mva6_nosyst/




## MVA training

#python analysis/tt5TeV/run.py cafea/json/5TeV/mva_split/tbarchannel_mva.json -n 32 -j -s 10000 -o tbarchannel_mva -p mva/lesscuts/   #these for the 30% part we take for the mva
#python analysis/tt5TeV/run.py cafea/json/5TeV/mva_split/tchannel_mva.json -n 32 -j -s 10000 -o tchannel_mva -p mva/lesscuts/ 

#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/WJetsToLNu.json -n 32 -j  -s 10000   -o Wjets_mva -p mva/lesscuts/ #These two for the bckg training (note TT has no xsec since we just use this for training, not for analysis, so there's no need to normalize) 
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/TT.json -n 32 -j  -s 10000   -o TT_mva -p mva/lesscuts/




#TCHANNEL

#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part0.json -n 32 -j  -s 10000   -o tbarchannel_part0 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part1.json -n 32 -j   -s 10000  -o tbarchannel_part1 -p mva6_nosyst/  #these if we want to run all the rootfiles of tchannel and tbarchannel
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part2.json -n 32 -j  -s 10000   -o tbarchannel_part2 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part3.json -n 32 -j   -s 10000  -o tbarchannel_part3 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part4.json -n 32 -j   -s 10000  -o tbarchannel_part4 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part5.json -n 32 -j  -s 10000   -o tbarchannel_part5 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tbarchannel_part6.json -n 32 -j  -s 10000   -o tbarchannel_part6 -p mva6_nosyst/

#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part0.json -n 32 -j   -s 10000  -o tchannel_part0 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part1.json -n 32 -j   -s 10000  -o tchannel_part1 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part2.json -n 32 -j  -s 10000   -o tchannel_part2 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part3.json -n 32 -j   -s 10000  -o tchannel_part3 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part4.json -n 32 -j   -s 10000  -o tchannel_part4 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part5.json -n 32 -j  -s 10000   -o tchannel_part5 -p mva6_nosyst/
#python analysis/tt5TeV/run.py cafea/json/5TeV/newxsecs/check_vicor/tchannel_part6.json -n 32 -j   -s 10000  -o tchannel_part6 -p mva6_nosyst/
