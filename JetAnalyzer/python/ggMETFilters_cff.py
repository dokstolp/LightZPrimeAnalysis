import FWCore.ParameterSet.Config as cms

from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = cms.bool(True)

from RecoMET.METFilters.globalTightHalo2016Filter_cfi import *
globalTightHalo2016Filter.taggingMode = cms.bool(True)

from RecoMET.METFilters.primaryVertexFilter_cfi import *
primaryVertexFilter.vertexCollection = cms.InputTag('offlinePrimaryVertices')

#primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                   minimumNDOF = cms.uint32(4) ,
#                                   maxAbsZ = cms.double(24),
#                                   maxd0 = cms.double(2),
#                                   filter = cms.bool(False)
#                                   )

#goodVertices = cms.EDFilter(
#  "VertexSelector",
#  filter = cms.bool(True),
#  src = cms.InputTag("offlinePrimaryVertices"),
#  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
#)

from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

ggMETFiltersSequence = cms.Sequence(
     HBHENoiseFilterResultProducer*
     primaryVertexFilter*
     eeBadScFilter*
     globalTightHalo2016Filter*
     EcalDeadCellTriggerPrimitiveFilter
)
