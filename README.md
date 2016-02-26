# LightZPrimeAnalysis

Yang Bai and collaborators have a model in which there exists a 
light Z', which only has enough mass to decay to a pair of oppositely
charged light hadrons (pi+/-, K+/-), with an additional
dark matter candidate.  The process to study is:

<bash>
p p -> Z', DM -> pi+, pi-, MET
</bash>

This analyzer looks for highly boosted very low invariant mass
oppositely charged hadrons recoiling from large MET.

Instructions:

<bash>
cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/SridharaDasu/LightZPrimeAnalysis.git
cd LightZPrimeAnalysis/LightZPrimeAnalyzer/test
cmsRun testLightZPrimeAnalyzer.py >& test.log
</bash>

Examine test.log and the root file created for your amusement.

If you are amused, work away and find it in 13-TeV LHC data we are collecting
