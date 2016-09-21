# LightZPrimeAnalysis

Yang Bai and collaborators have a model in which there exists a 
light Z', which only has enough mass to decay to a small number of
hadrons, perhaps, even a single pair of oppositely charged light 
hadrons (pi+/-, K+/-).  However, it is produced with an additional 
dark matter candidate.  For moderate masses (1 < MZprime < 10 GeV),
we get very low multiplicity "pencil" jets.  The process to study is:

```bash
p p -> Z', DM, DMBar -> u, uBar, MET -> "pencil jet", MET
```

LightZPrimeGenerator package contains a simple python file to
convert LHE files from MadGraph to CMSSW GEN formatted output.
Note that this does not do simulation step.

The first analyzer LightZPrimeGenAnalyzer makes an nTuple with
information from the generator level for testing.

The second analyzer JetAnalyzer makes simple calculations of
the energy weighted eta-, phi- width of jets.

The third analyzer LightZPrimeAnalyzer is still in development.
Its function is to impleent a real analysis running on MiniAOD.

Instructions:

```bash
cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/uhussain/LightZPrimeAnalysis.git
cd $CMSSW_BASE/src/LightZPrimeAnalysis/LightZPrimeGenerator/test
cmsRun LightZPrimeDecayAndHadronizer.py
cd $CMSSW_BASE/src/LightZPrimeAnalysis/LightZPrimeAnalyzer/test
cmsRun testLightZPrimeGenAnalyzer.py
cd $CMSSW_BASE/src/LightZPrimeAnalysis/JetAnalyzer/test
cmsRun testJetAnalyzer.py
```

Examine the root file created for your amusement.

If you are amused, work away and find it in 13-TeV LHC data we are collecting

Instructions (in order to apply all recommended MET Filters):
```bash
cmsrel CMSSW_8_0_X
cmsenv
git cms-init
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
git clone -b usamabranch https://github.com/uhussain/LightZPrimeAnalysis.git
cd CMSSW_8_0_X/src/LightZPrimeAnalysis/JetAnalyzer/test
cmsRun testJetAnalyzer.py
```
https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Details about the application of the filters:
