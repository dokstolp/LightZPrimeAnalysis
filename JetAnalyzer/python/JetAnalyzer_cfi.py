import FWCore.ParameterSet.Config as cms

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

process.ak4PFchsCorrectedJets   = cms.EDProducer('CorrectedPFJetProducer',
    src         = cms.InputTag('ak4PFJetsCHS'),
    correctors  = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')
    )

process.jetSequence = cms.Sequence(process.ak4PFchsCorrectedJets*process.ak4PFCHSL1FastL2L3CorrectorChain)

JetAnalyzer = cms.EDAnalyzer('JetAnalyzer',
	jetsTag = cms.InputTag("ak4PFchsCorrectedJets"),
)
