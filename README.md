# LightZPrimeAnalysis

Yang Bai and collaborators have a model in which a light Z' 
which only has enough mass to decay to a pair of oppositely
charged light hadrons (pi+/-, K+/-), with an additional
dark matter candidate.  The process to study is:

p p -> Z', DM -> pi+, pi-, MET

This analyzer looks for highly boosted very low invariant mass
oppositely charged hadrons recoiling from large MET.

Instructions:

cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/SridharaDasu/LightZPrimeAnalysis.git
cd LightZPrimeAnalysis/LightZPrimeAnalyzer/test
cmsRun testLightZPrimeAnalyzer.py >& test.log

Examine test.log and the root file created for your amusement.
