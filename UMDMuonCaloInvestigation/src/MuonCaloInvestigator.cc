// -*- C++ -*-
//
// Package:    L1TMuon
// Class:      MuonCaloInvestigator
// 
/**\class MuonCaloInvestigator MuonCaloInvestigator.cc L1TriggerDPGUpgrade/L1TMuon/plugins/MuonCaloInvestigator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Belloni
//         Created:  Wed, 12 Jun 2013 16:42:12 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "DataFormats/HcalRecHit/interface/CastorRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalCalibRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"

#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace L1TMuon;

//
// class declaration
//

class MuonCaloInvestigator : public edm::EDAnalyzer {
public:
  explicit MuonCaloInvestigator(const edm::ParameterSet&);
  ~MuonCaloInvestigator();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  edm::Service<TFileService> _fileService;

  edm::ESHandle<CaloTPGTranscoder> _caloDecoder;
  
  // Surfaces to be used for extrapolation
  Cylinder::CylinderPointer _rpcCyl[4]; // 4 = number of stations
  Cylinder::CylinderPointer _hoCyl;
  
  bool _doGen;
  edm::InputTag _genInput;

  edm::InputTag _rpcInput;
  edm::InputTag _dttfInput;
  edm::InputTag _hcalInput;
  edm::InputTag _stdmuInput;
  edm::InputTag _glbmuInput;

  // Delta R values to be used for matching
  double _dRtruthToRpc, _dRrpcToDttf, _dRdttfToHcal, 
    _dRhcalToStdMu, _dRdttfToStdMu;
  edm::InputTag _truthToRpc;
  edm::InputTag _rpcToDttf;
  edm::InputTag _dttfToHcal;
  edm::InputTag _hcalToStdMu;
  edm::InputTag _dttfToStdMu; // or global: can't be too different

  // map with histograms: all deltaEta and deltaPhi plots will
  // have same boundaries (very generous), then work out useful
  // ranges with plotting macro
  std::map<std::string,TH1F*> _h1dDeltaEta;
  std::map<std::string,TH1F*> _h1dDeltaPhi;
  std::map<std::string,TH1F*> _h1dDeltaR;
  std::map<std::string,TH2F*> _h2dDeltaEtaPhi;

  std::map<std::string,TH1F*> _h1dEta;
  std::map<std::string,TH1F*> _h1dPhi;
  std::map<std::string,TH1F*> _h1dPt;

  std::map<std::string,TH2F*> _h2dEtaPhi;
  std::map<std::string,TH2F*> _h2dXY;

  std::map<std::string,TH1F*> _h1dStations;

  TH1F* _counters;
  enum { ALL=0, TRUTH, RPC, DTTF, HCAL, STDMU, GLBMU };
  
  TriggerPrimitiveRef getBestTriggerPrimitive
  (const TriggerPrimitiveList& list, unsigned subsystem) const;

  void fillDeltaEtaPhiHistograms(float eta1, float phi1,
				 float eta2, float phi2,
				 std::string key);
  void fillMapHistograms(float eta, float phi,
			 float x, float y,
			 std::string key);
  void fillKinematicHistograms(float eta, float phi, float pt,
			       std::string key);

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonCaloInvestigator::MuonCaloInvestigator(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  if( (_doGen = iConfig.getUntrackedParameter<bool>("doGen",false)) ) {
    _genInput = iConfig.getParameter<edm::InputTag>("genSrc");
  }
  _rpcInput =   iConfig.getParameter<edm::InputTag>("rpcSrc");
  _dttfInput =  iConfig.getParameter<edm::InputTag>("dttfSrc");
  _hcalInput =  iConfig.getParameter<edm::InputTag>("hcalSrc");
  _stdmuInput = iConfig.getParameter<edm::InputTag>("stdmuSrc");
  _glbmuInput = iConfig.getParameter<edm::InputTag>("glbmuSrc");

  _dRtruthToRpc =iConfig.getUntrackedParameter<double>("dRtruthToRpc" ,1.);
  _dRrpcToDttf  =iConfig.getUntrackedParameter<double>("dRrpcToDttf"  ,1.);
  _dRdttfToHcal =iConfig.getUntrackedParameter<double>("dRdttfToHcal" ,1.);
  _dRhcalToStdMu=iConfig.getUntrackedParameter<double>("dRhcalToStdMu",1.);
  _dRdttfToStdMu=iConfig.getUntrackedParameter<double>("dRdttfToStdMu",1.);

  _counters = _fileService->make<TH1F>("counter","counters",10,-0.5,9.5);

  // Build surfaces for extrapolation
  Cylinder::PositionType pos0; // mah, we do not use them, but they are
  Cylinder::RotationType rot0; // needed in Cylinder constructor
  _rpcCyl[0] = Cylinder::build(4000,pos0,rot0);
  _rpcCyl[1] = Cylinder::build(5000,pos0,rot0);
  _rpcCyl[2] = Cylinder::build(6000,pos0,rot0);
  _rpcCyl[3] = Cylinder::build(7000,pos0,rot0);

  _hoCyl = Cylinder::build(3600,pos0,rot0);

}


MuonCaloInvestigator::~MuonCaloInvestigator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonCaloInvestigator::analyze(const edm::Event& iEvent, 
			      const edm::EventSetup& iSetup)
{

  // Raise flags to indicate matching
  bool foundTruth = false;
  bool foundRpc   = false;
  bool foundDttf  = false;
  bool foundHcal  = false;
  bool foundStdMu = false;
  
  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<IdealGeometryRecord>().get(htopo);

  edm::ESHandle<CaloGeometry> hgeome;
  iSetup.get<CaloGeometryRecord>().get(hgeome);

  // Setup the decoder
  iSetup.get<CaloTPGRecord>().get(_caloDecoder);
  //_caloDecoder->setup(iSetup, CaloTPGTranscoder::HcalTPG);

  // Setup the B field
  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  // Setup the tracking geometry
  edm::ESHandle<GlobalTrackingGeometry> htrackgeo;
  iSetup.get<GlobalTrackingGeometryRecord>().get(htrackgeo);

  // Setup the track propagator
  edm::ESHandle<Propagator> shProp;
  iSetup.get<TrackingComponentsRecord>().
    get("SteppingHelixPropagatorAlong", shProp);

  // Here we open all collections!
  edm::Handle<reco::GenParticleCollection> truthParticles;
   // I will get it only if needed, later on!
  //iEvent.getByLabel(_genInput,truthParticles);

  edm::Handle<InternalTrackCollection> rpcTriggerPrimitives;
  iEvent.getByLabel(_rpcInput,rpcTriggerPrimitives);

  edm::Handle<InternalTrackCollection> dttfTriggerPrimitives;
  iEvent.getByLabel(_dttfInput,dttfTriggerPrimitives);

  edm::Handle<TriggerPrimitiveCollection> hcalTriggerPrimitives;
  iEvent.getByLabel(_hcalInput,hcalTriggerPrimitives);

  edm::Handle<reco::TrackCollection> standaloneMuons;
  iEvent.getByLabel(_stdmuInput,standaloneMuons);

  edm::Handle<reco::TrackCollection> globalMuons;
  iEvent.getByLabel(_glbmuInput,globalMuons);

  // Some extra collection, always interesting
  edm::Handle<CaloTowerCollection> caloTowers;
  iEvent.getByLabel("towerMaker",caloTowers);

  edm::Handle<HORecHitCollection> hoRecoHits;
  iEvent.getByLabel("horeco",hoRecoHits);

  // Weird idea: if running on data, let us use the global muons
  // instead of truth muons in the outer loop
  // May have to use good old unsigned int loop on collections
  // instead of nice loop with iterators
  //auto btruth = _doGen ? truthParticles->cbegin() : globalMuons->cbegin();
  //auto etruth = _doGen ? truthParticles->cend() : globalMuons->cend();
  //if (_doGen) {
  iEvent.getByLabel(_genInput,truthParticles);
  auto btruth = truthParticles->cbegin();
  auto etruth = truthParticles->cend();
  //}

  for( ; btruth != etruth; ++btruth ) {
    // Initial quality cuts
    if (_doGen) {
      if (std::abs(btruth->pdgId()) != 13 ||
	  btruth->pt()<20.)
	continue;
    }
    else {
      if (btruth->pt()<20.)
	continue;
    }
    foundTruth = true;
    
    fillKinematicHistograms(btruth->eta(),btruth->phi(),btruth->pt(),"truth");

    // end of quality cuts, let's start with RPC
    auto brpc = rpcTriggerPrimitives->cbegin();
    auto erpc = rpcTriggerPrimitives->cend();
    for( ; brpc != erpc; ++brpc ) {
      double rpcEta=0.,rpcPhi=0.;
      int rpcStat=0; // number of RPC stations with a trigger primitive
      TriggerPrimitiveStationMap stubs = brpc->getStubs();
      // Chris' function will return configuration of RPC trigger
      // primitives: how many and which stations
      //int rpctype = getType(stubs);
 
      // Loop on RPC stations, for each of which you can get a TP
      unsigned station;
      for( station = 1; station <= 4; ++station ) {
	const unsigned idx = 4*1+station-1; // RPCb=1
	if( !stubs.count(idx) ) continue;
	TriggerPrimitiveList tpRpcB = stubs[idx];
	TriggerPrimitiveRef bestRpcB = getBestTriggerPrimitive(tpRpcB,1);
	rpcEta+=bestRpcB->getCMSGlobalEta();
	rpcPhi+=bestRpcB->getCMSGlobalPhi();
	++rpcStat;
      }
      if (rpcStat>0) {
	rpcEta/=rpcStat;
	rpcPhi/=rpcStat;
      }
      // Let's try to extrapolate to RPC region 
      // the truth tracks
      FreeTrajectoryState initial(GlobalPoint(btruth->vx(),
					      btruth->vy(),
					      btruth->vz()),
				  GlobalVector(btruth->px(),
					       btruth->py(),
					       btruth->pz()),
				  btruth->charge(),
				  &*bField);
      
      TrajectoryStateOnSurface final = 
	shProp->propagate(initial,*_rpcCyl[station-1]);

      // Let us fill some histograms!
      fillDeltaEtaPhiHistograms(btruth->eta(),btruth->phi(),
				rpcEta,rpcPhi,
				"truth-rpc");
      if (final.isValid())
	fillDeltaEtaPhiHistograms(final.globalPosition().eta(),
				  final.globalPosition().phi(),
				  rpcEta,rpcPhi,
				  "extratruth-rpc");
      if (btruth == truthParticles->cbegin()) {
      	fillMapHistograms(rpcEta,rpcPhi,0.,0.,
      			  "rpc");
      }
      ///////////////////////////////////////
      
      // Continue with matching only if we did find a match
      if (sqrt(reco::deltaR2(btruth->eta(),btruth->phi(),
			     rpcEta,rpcPhi))>_dRtruthToRpc)
	continue;
      foundRpc = true;

      // Now loop on DTTF tracks...
      auto bdttf = dttfTriggerPrimitives->cbegin();
      auto edttf = dttfTriggerPrimitives->cend();
      for( ; bdttf != edttf; ++bdttf ) {
	double dttfEta=0.,dttfPhi=0.;
	int dttfStat=0; // number of DTTF stations with a trigger primitive
	TriggerPrimitiveStationMap stubs = bdttf->getStubs();
	// Chris' function will return configuration of DTTF trigger
	// primitives: how many and which stations
	//int dttftype = getType(stubs);
	
	// Loop on DTTF stations, for each of which you can get a TP
	unsigned station;
	for( station = 1; station <= 4; ++station ) {
	  const unsigned idx = 4*0+station-1; // DTTF=0
	  if( !stubs.count(idx) ) continue;
	  TriggerPrimitiveList tpDttf = stubs[idx];
	  TriggerPrimitiveRef bestDttf = getBestTriggerPrimitive(tpDttf,0);
	  dttfEta+=bestDttf->getCMSGlobalEta();
	  dttfPhi+=bestDttf->getCMSGlobalPhi();
	  ++dttfStat;
	}
	if (dttfStat>0) {
	  dttfEta/=dttfStat;
	  dttfPhi/=dttfStat;
	}
	// Let us fill some histograms here too
	// Carefull that we are inside a double-loop
	// I want to fill the RPC-DTTF plot only once per truth muon
	// and, similarly, fill the TRUTH-DTTF only once per RPC stub 
	if (brpc == rpcTriggerPrimitives->cbegin()) {
	  fillDeltaEtaPhiHistograms(btruth->eta(),btruth->phi(),
				    dttfEta,dttfPhi,
				    "truth-dttf");
	}
	if (btruth == truthParticles->cbegin()) {
	  fillDeltaEtaPhiHistograms(rpcEta,rpcPhi,
				    dttfEta,dttfPhi,
				    "rpc-dttf");
	}
	if (btruth == truthParticles->cbegin() &&
	    brpc == rpcTriggerPrimitives->cbegin()) {
	  fillMapHistograms(dttfEta,dttfPhi,0.,0.,
			    "dttf");
	}
	

	// Check: what is the eta of the internal track vs. the
	// one of the trigger primitive?
	//if (btruth == truthParticles->cbegin() &&
	//    brpc == rpcTriggerPrimitives->cbegin()) {
	//  fillDeltaEtaPhiHistograms(bdttf->parent()->etaValue(),
	//			    bdttf->parent()->phiValue(),
	//			    dttfEta,dttfPhi,
	//			    "dttfTK-dttfTP");
	//}
	///////////////////////////////////////

	// Continue with matching only if we did find a match
	if (sqrt(reco::deltaR2(rpcEta,rpcPhi,
			       dttfEta,dttfPhi))>_dRrpcToDttf)
	  continue;
	foundDttf = true;

	// One other layer of complication: loop on HCAL TP
	auto bhcal = hcalTriggerPrimitives->cbegin();
	auto ehcal = hcalTriggerPrimitives->cend();
	for( ; bhcal != ehcal; ++bhcal ) {
	  //TriggerPrimitiveStationMap stubs = bhcal->getStubs();
	  //const unsigned idx = 4*4+1-1; // HCAL=4, station=1
	  //if( !stubs.count(idx) ) continue;
	  //TriggerPrimitiveList tpHcal = stubs[idx];
	  //TriggerPrimitiveRef bestHcal = 
	  //  getBestTriggerPrimitive(tpHcal,4);
	  if (brpc == rpcTriggerPrimitives->cbegin() &&
	      bdttf == dttfTriggerPrimitives->cbegin() ) {
	    fillDeltaEtaPhiHistograms(btruth->eta(),btruth->phi(),
				      bhcal->getCMSGlobalEta(),
				      bhcal->getCMSGlobalPhi(),
				      "truth-hcal");
	  }
	  if (btruth == truthParticles->cbegin() &&
	      bdttf == dttfTriggerPrimitives->cbegin() ) {
	    if (fabs(bhcal->getCMSGlobalEta())<1.3)
	      fillDeltaEtaPhiHistograms(rpcEta,rpcPhi,
					bhcal->getCMSGlobalEta(),
					bhcal->getCMSGlobalPhi(),
					"rpc-hcal");
	  }
	  if (btruth == truthParticles->cbegin() &&
	      brpc == rpcTriggerPrimitives->cbegin() ) {
	    if (fabs(bhcal->getCMSGlobalEta())<1.3)
	      fillDeltaEtaPhiHistograms(dttfEta,dttfPhi,
					bhcal->getCMSGlobalEta(),
					bhcal->getCMSGlobalPhi(),
					"dttf-hcal");
	  }
	  if (btruth == truthParticles->cbegin() &&
	      brpc == rpcTriggerPrimitives->cbegin() &&
	      bdttf == dttfTriggerPrimitives->cbegin() ) {
	    fillMapHistograms(bhcal->getCMSGlobalEta(),
	  		      bhcal->getCMSGlobalPhi(),
	  		      0.,0.,
	  		      "hcal");
	  }
	  

	  // Notice that here I have a truth muon, an RPC TP,
	  // a DTTF track and an HCAL TP, not necessarily matching

	  // Continue with matching only if we did find a match
	  if (sqrt(reco::deltaR2(dttfEta,dttfPhi,
				 bhcal->getCMSGlobalEta(),
				 bhcal->getCMSGlobalPhi()))>_dRdttfToHcal)
	    continue;
	  foundHcal = true;
	  
	  //std::cout << _caloDecoder->hcaletValue
	  //  (bhcal->detId<HcalTrigTowerDetId>().ieta(),
	  //   bhcal->detId<HcalTrigTowerDetId>().iphi(),
	  //   bhcal->getHCALData().SOI_compressedEt)
	  //	    << std::endl;

	  // Let's loop on muons
	  // I will be sneaky and write this code only once
	  // then use it for standalone and global muons by just
	  // modifying the python file...
	  auto bstdmu = standaloneMuons->cbegin();
	  auto estdmu = standaloneMuons->cend();
	  for( ; bstdmu != estdmu; ++bstdmu ) {
	    if (brpc == rpcTriggerPrimitives->cbegin() &&
		bdttf == dttfTriggerPrimitives->cbegin() &&
		bhcal == hcalTriggerPrimitives->cbegin() ) {
	      fillDeltaEtaPhiHistograms(btruth->eta(),btruth->phi(),
					bstdmu->eta(),bstdmu->phi(),
					"truth-standalone");
	    }

	    // Continue with matching only if we did find a match
	    if (sqrt(reco::deltaR2(bhcal->getCMSGlobalEta(),
				   bhcal->getCMSGlobalPhi(),
				   bstdmu->eta(),
				   bstdmu->phi()))>_dRhcalToStdMu)
	      continue;
	    foundStdMu = true;
  
	    // Technically, here I have a truth muon matched to
	    // an RPC TP, matched to a DTTF TP, matched to an HCAL TP,
	    // matched to a standalone (or global, if I change the
	    // collection name) muon

	    // Let's loop on Calo Towers
	    // I want to check if there is an HORecHit in front of the
	    // global muon, and check its energy. How does that
	    // hit match with HcalTrigPrimitives? And RPC and DTTF?
	    // Notice that there will be tons of noisy towers
	    auto bcalo = caloTowers->begin();
	    auto ecalo = caloTowers->end();
	    for( ; bcalo != ecalo; ++bcalo ) {
	      if (btruth == truthParticles->cbegin() &&
		  brpc == rpcTriggerPrimitives->cbegin() &&
		  bdttf == dttfTriggerPrimitives->cbegin() &&
		  bhcal == hcalTriggerPrimitives->cbegin() ) {
		fillDeltaEtaPhiHistograms(bstdmu->eta(),bstdmu->phi(),
					  bcalo->hadPosition().eta(),
					  bcalo->hadPosition().phi(),
					  "standalone-calotower");
	      }
	      if (brpc == rpcTriggerPrimitives->cbegin() &&
		  bdttf == dttfTriggerPrimitives->cbegin() &&
		  bhcal == hcalTriggerPrimitives->cbegin() &&
		  bstdmu == standaloneMuons->cbegin() ) {
		fillDeltaEtaPhiHistograms(btruth->eta(),btruth->phi(),
					  bcalo->hadPosition().eta(),
					  bcalo->hadPosition().phi(),
					  "truth-calotower");
	      }
	      if (btruth == truthParticles->cbegin() &&
		  bdttf == dttfTriggerPrimitives->cbegin() &&
		  bhcal == hcalTriggerPrimitives->cbegin() &&
		  bstdmu == standaloneMuons->cbegin() ) {
		fillDeltaEtaPhiHistograms(rpcEta,rpcPhi,
					  bcalo->hadPosition().eta(),
					  bcalo->hadPosition().phi(),
					  "rpc-calotower");
	      }
	      if (btruth == truthParticles->cbegin() &&
		  brpc == rpcTriggerPrimitives->cbegin() &&
		  bhcal == hcalTriggerPrimitives->cbegin() &&
		  bstdmu == standaloneMuons->cbegin() ) {
		fillDeltaEtaPhiHistograms(dttfEta,dttfPhi,
					  bcalo->hadPosition().eta(),
					  bcalo->hadPosition().phi(),
					  "dttf-calotower");
	      }
	      if (btruth == truthParticles->cbegin() &&
		  brpc == rpcTriggerPrimitives->cbegin() &&
		  bdttf == dttfTriggerPrimitives->cbegin() &&
		  bstdmu == standaloneMuons->cbegin() ) {
		fillDeltaEtaPhiHistograms(bhcal->getCMSGlobalEta(),
					  bhcal->getCMSGlobalPhi(),
					  bcalo->hadPosition().eta(),
					  bcalo->hadPosition().phi(),
					  "hcal-calotower");
	      }
	      if (btruth == truthParticles->cbegin() &&
	      	  brpc == rpcTriggerPrimitives->cbegin() &&
	      	  bdttf == dttfTriggerPrimitives->cbegin() &&
	      	  bhcal == hcalTriggerPrimitives->cbegin() &&
	      	  bstdmu == standaloneMuons->cbegin() ) {
		if (fabs(bcalo->hadPosition().eta())<1.3)
		  fillMapHistograms(bcalo->hadPosition().eta(),
				    bcalo->hadPosition().phi(),
				    bcalo->hadPosition().x(),
				    bcalo->hadPosition().y(),
				    "calotower");
	      }
	    } // end loop on calo towers
	  } // end loop on standalone muons
	} // end loop on HCAL TP
      } // end loop on DTTF
    } // end loop on RPC
  } // end loop on truth or global

  // Here we fill counters - notice that I have to avoid leaving
  // the function analyze before I get here, otherwise the counts
  // will not be correct
  if (true)       _counters->Fill(ALL);
  if (foundTruth) _counters->Fill(TRUTH);
  if (foundRpc)   _counters->Fill(RPC);
  if (foundDttf)  _counters->Fill(DTTF);
  if (foundHcal)  _counters->Fill(HCAL);
  if (foundStdMu) _counters->Fill(STDMU);

}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonCaloInvestigator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonCaloInvestigator::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonCaloInvestigator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonCaloInvestigator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonCaloInvestigator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonCaloInvestigator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonCaloInvestigator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

TriggerPrimitiveRef MuonCaloInvestigator::
getBestTriggerPrimitive(const TriggerPrimitiveList& list, 
			unsigned subsystem) const {
  TriggerPrimitiveRef result;
  unsigned bestquality = 0, qualtemp; // for CSCs / DTs / HCAL (just for fun)
  float phiavg, bestdphi, lsize; // average strip for RPCS
  auto tp = list.cbegin();
  auto tpend = list.cend();
  
  switch( subsystem ) {
  case 0: // DTs
    for( ; tp != tpend; ++tp ) {
      qualtemp = 0;
      if( (*tp)->getDTData().qualityCode != -1 ) {
	qualtemp += (~((*tp)->getDTData().qualityCode)&0x7) << 1;	
      }
      if( (*tp)->getDTData().theta_quality != -1 ) {
	qualtemp += (~((*tp)->getDTData().theta_quality)&0x1);
      }
      if( qualtemp > bestquality ) {
	bestquality = qualtemp;
	result = *tp;
      }
    }
      break;
  case 2: // CSCs
    for( ; tp != tpend; ++tp ) {
      qualtemp = (*tp)->getCSCData().quality;      
      if ( qualtemp > bestquality ) {
	bestquality = qualtemp;
	result = *tp;
      }
    }
    break;
  case 1:
  case 3: // RPCb/f
    phiavg = 0;
    lsize = list.size();
    for( ; tp != tpend; ++tp ) {
      phiavg += (*tp)->getCMSGlobalPhi();
    }
    phiavg = phiavg/lsize;    
    tp = list.cbegin();
    bestdphi = 100;
    for( ; tp != tpend; ++tp ) {      
      if( std::abs((*tp)->getCMSGlobalPhi() - phiavg) < bestdphi ) {
	result = *tp;
	bestdphi = std::abs((*tp)->getCMSGlobalPhi() - phiavg);
      }
    }
    break;
  case 4: // HCAL
    for( ; tp != tpend; ++tp ) {
      qualtemp = (*tp)->getHCALData().size;
      if ( qualtemp > bestquality ) {
	bestquality = qualtemp;
	result = *tp;
      }
    }
    break;
  default:
    break;
  }
  return result;
}

void MuonCaloInvestigator::fillDeltaEtaPhiHistograms(float eta1, float phi1,
						     float eta2, float phi2,
						     std::string key) {
  
  if(!_h1dDeltaEta.count(key))
    _h1dDeltaEta[key] = 
      _fileService->make<TH1F>(Form("deta_%s",key.c_str()),
			       Form("#Delta#eta %s",key.c_str()),
			       500,-0.5,0.5);
  _h1dDeltaEta[key]->Fill(eta1-eta2);

  if(!_h1dDeltaPhi.count(key))
    _h1dDeltaPhi[key] = 
      _fileService->make<TH1F>(Form("dphi_%s",key.c_str()),
			       Form("#Delta#phi %s",key.c_str()),
			       500,-M_PI/10.,M_PI/10.);
  _h1dDeltaPhi[key]->Fill(phi1-phi2);
  
  if(!_h1dDeltaR.count(key))
    _h1dDeltaR[key] = 
      _fileService->make<TH1F>(Form("dR_%s",key.c_str()),
			       Form("#Delta R %s",key.c_str()),
			       500,0,1.0);
  _h1dDeltaR[key]->Fill(sqrt(reco::deltaR2(eta1,phi1,eta2,phi2)));

  if(!_h2dDeltaEtaPhi.count(key))
    _h2dDeltaEtaPhi[key] = 
      _fileService->make<TH2F>(Form("detaphi_%s",key.c_str()),
			       Form("#Delta#phi vs. #Delta#eta %s",
				    key.c_str()),
			       500,-0.5,0.5,
			       500,-M_PI/10.,M_PI/10.);
  _h2dDeltaEtaPhi[key]->Fill(eta1-eta2,
			     phi1-phi2);
  
  return;
}

void MuonCaloInvestigator::fillMapHistograms(float eta, float phi,
					     float x, float y,
					     std::string key) {
  
  if(!_h2dEtaPhi.count(key))
    _h2dEtaPhi[key] = 
      _fileService->make<TH2F>(Form("etaphi_%s",key.c_str()),
			       Form("#phi vs. #eta %s",
				    key.c_str()),
			       500,-1.3,1.3,
			       500,-M_PI,M_PI);
  _h2dEtaPhi[key]->Fill(eta,phi);

  if(!_h2dXY.count(key))
    _h2dXY[key] = 
      _fileService->make<TH2F>(Form("xy_%s",key.c_str()),
			       Form("Y vs. X %s",
				    key.c_str()),
			       2050,-4100,4100,
			       2050,-4100,4100);
  _h2dXY[key]->Fill(x,y);
  
  return;
}


void MuonCaloInvestigator::fillKinematicHistograms(float eta, float phi,
						   float pt,
						   std::string key) {
  
  if(!_h1dEta.count(key))
    _h1dEta[key] = 
      _fileService->make<TH1F>(Form("eta_%s",key.c_str()),
			       Form("#eta %s",
				    key.c_str()),
			       500,-1.3,1.3);
  _h1dEta[key]->Fill(eta);
  
  if(!_h1dPhi.count(key))
    _h1dPhi[key] = 
      _fileService->make<TH1F>(Form("phi_%s",key.c_str()),
			       Form("#phi %s",
				    key.c_str()),
			       500,-M_PI,M_PI);
  _h1dPhi[key]->Fill(phi);
  
  if(!_h1dPt.count(key))
    _h1dPt[key] = 
      _fileService->make<TH1F>(Form("pt_%s",key.c_str()),
			       Form("#pt %s",
				    key.c_str()),
			       500,0,250);
  _h1dPt[key]->Fill(pt);



  
  return;
}




//define this as a plug-in
DEFINE_FWK_MODULE(MuonCaloInvestigator);
