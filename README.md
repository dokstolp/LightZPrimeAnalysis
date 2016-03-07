# LightZPrimeAnalysis

Yang Bai and collaborators have a model in which there exists a 
light Z', which only has enough mass to decay to a pair of oppositely
charged light hadrons (pi+/-, K+/-), with an additional
dark matter candidate.  For moderate masses (1 < MZprime < 10 GeV),
we get low multiplicity jets.  The process to study is:

```bash
p p -> Z', DM, DMBar -> u, uBar, MET
```

LightZPrimeGenerator package contains a simple python file to
convert LHE files from MadGraph to CMSSW GEN formatted output.
Note that this does not do simulation step.

The first analyzer LightZPrimeGenAnalyzer prints out some event level
information from the generator level for testing.

The second analyzer LightZPrimeAnalyzer is still in development.
Its function is to impleent a real analysis running on MiniAOD.

Instructions:

```bash
cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/SridharaDasu/LightZPrimeAnalysis.git
cd $CMSSW_BASE/src/LightZPrimeAnalysis/LightZPrimeGenerator/test
cmsRun LightZPrimeDecayAndHadronizer.py
cd $CMSSW_BASE/src/LightZPrimeAnalysis/LightZPrimeAnalyzer/test
cmsRun testLightZPrimeGenAnalyzer.py
```

Examine the root file created for your amusement.

If you are amused, work away and find it in 13-TeV LHC data we are collecting
