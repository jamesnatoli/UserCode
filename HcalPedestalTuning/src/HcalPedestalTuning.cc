// -*- C++ -*-
//
// Package:    HcalPedestalTuning
// Class:      HcalPedestalTuning
// 
/**\class HcalPedestalTuning HcalPedestalTuning.cc UserCode/HcalPedestalTuning/src/HcalPedestalTuning.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri Jan 21 14:38:40 EST 2011
// $Id: HcalPedestalTuning.cc,v 1.2 2011/02/21 21:38:45 lockner Exp $
// $Id: HcalPedestalTuning.cc,v 2.0 2013/08/08 11:33:55 abelloni, cranelli Exp $
//
//

// system include files
#include <boost/filesystem.hpp>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"

#include "FWCore/Utilities/interface/tinyxml.h"

#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "DataFormats/HcalDetId/interface/HcalFrontEndId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "CondFormats/HcalObjects/interface/HcalLogicalMap.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalLogicalMapGenerator.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include <fstream>

using namespace std;

//
// class declaration
//

class HcalPedestalTuning : public edm::EDAnalyzer {
public:
  explicit HcalPedestalTuning(const edm::ParameterSet&);
  ~HcalPedestalTuning();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void readSettings();
  void tunePedestals(const double targetPedestalHPD, 
		     const double targetPedestalSiPM);
  void writeOutput();
  void writeUnstable(const string&);
  void writeTuned(const string&);
  HcalSubdetector subdet(const string&);
  string subdet(const HcalSubdetector&);

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;

  TH1D * HB_1Dpeds; 
  TH1D * HE_1Dpeds; 
  TH1D * HO_1Dpeds; 
  TH1D * HF_1Dpeds; 

  TH2D * HB_2Dpeds; 
  TH2D * HE_2Dpeds; 
  TH2D * HO_2Dpeds; 
  TH2D * HF_2Dpeds; 

  // run number
  int run;
  string run_string;
  // structure that contains channel information needed for pedestal tuning
  class ChannelData {
  public:
    bool isHOX;
    bool isCalib;
    bool isSiPM;
    uint32_t rawId;
    HcalFrontEndId frontEndId;
    int dacSetting;

    double capSum[4];
    int capEntries[4];

    double getCapPedestal(int capId) const;
    double getChannelPedestal() const;
    double getStability() const;

    static int digiSize;
    static int evtCounter;
  };

  // list of HCAL channels
  map<uint32_t, ChannelData> channelList;
  // map to go from the fake HcalFrontEndId to the rawId 
  map<HcalFrontEndId, uint32_t> rawIdMap;
  // HCAL logical map
  HcalLogicalMap *lMap;
  // needed to initialize the list of HCAL channels at the start of the job
  bool jobStart;
  // needed to determine the size of HCAL Digis
  bool digiStart;
  // InputTag to read HBHE, HO, HF and CalibHits
  string digisName_;
  // Interval Of Validity of topology map
  int mapIOV_;
  // switch on/off pedestal tuning
  bool tunePedestals_;
  // path to the old XML bricks
  string pathToOldXMLBricks_;
  // output directory for the new XML bricks
  string xmlBricksOutputDir_;
  // file with pedestal settings
  string settingsFile;
  // flag that specifies if all HCAL channels were set to a same DAC setting  
  bool uniformDACSettings_;
  // DAC setting to which all HCAL channels were set (if that was the case)
  int initialDACSetting_;
  // list of SiPM RBXes
  vector<string> SiPMRBXes_;
  // DAC setting which all SiMP channels will be set to 
  // int SiPMDACSetting_;
  // target pedestal value
  double targetPedestalHPD_;
  double targetPedestalSiPM_;
  // average change in ADC counts for the DAC change of 1
  double dac2adc_;
  // tag name
  string tagName_;
  // list of HCAL channels organized by RBX
  map<string, vector<ChannelData> > channelListByRBX;
  // list of untuned HCAL channels organized by RBX
  map<string, vector<ChannelData> > tunedChannelsByRBX;
  // list of missing HCAL channels organized by RBX
  map<string, vector<ChannelData> > unstableChannelsByRBX;
  // method to return true/false for sorting channels
  class compareRBXInfo {
  public:
    bool operator() (const ChannelData& a, const ChannelData& b) {
      if (a.frontEndId.rm() < b.frontEndId.rm()) return true;
      else if (a.frontEndId.rm() > b.frontEndId.rm()) return false;
      else { //a and b are in the same RM
	if (a.frontEndId.qieCard() < b.frontEndId.qieCard()) 
	  return true;
	else if (a.frontEndId.qieCard() > 
		 b.frontEndId.qieCard()) return false;
	else{ // a and b are in the same card
	  if (a.frontEndId.adc() < b.frontEndId.adc()) return true;
	  else return false;
	}
      }
    }
  };
  // HCAL digis
  edm::Handle<HBHEDigiCollection> hbheDigis;
  edm::Handle<HODigiCollection> hoDigis;
  edm::Handle<HFDigiCollection> hfDigis;
  edm::Handle<HcalCalibDigiCollection> calibDigis;

};

int HcalPedestalTuning::ChannelData::digiSize = 10;
int HcalPedestalTuning::ChannelData::evtCounter = 0;

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static string RBXes[132] = {
  "HBM01", "HBM02", "HBM03", "HBM04", "HBM05", "HBM06", "HBM07", "HBM08", 
  "HBM09", "HBM10", "HBM11", "HBM12", "HBM13", "HBM14", "HBM15", "HBM16", 
  "HBM17", "HBM18",
  "HBP01", "HBP02", "HBP03", "HBP04", "HBP05", "HBP06", "HBP07", "HBP08", 
  "HBP09", "HBP10", "HBP11", "HBP12", "HBP13", "HBP14", "HBP15", "HBP16", 
  "HBP17", "HBP18",
  "HEM01", "HEM02", "HEM03", "HEM04", "HEM05", "HEM06", "HEM07", "HEM08",
  "HEM09", "HEM10", "HEM11", "HEM12", "HEM13", "HEM14", "HEM15", "HEM16", 
  "HEM17", "HEM18",
  "HEP01", "HEP02", "HEP03", "HEP04", "HEP05", "HEP06", "HEP07", "HEP08",
  "HEP09", "HEP10", "HEP11", "HEP12", "HEP13", "HEP14", "HEP15", "HEP16", 
  "HEP17", "HEP18",
  "HFM01", "HFM02", "HFM03", "HFM04", "HFM05", "HFM06", 
  "HFM07", "HFM08", "HFM09", "HFM10", "HFM11", "HFM12", 
  "HFP01", "HFP02", "HFP03", "HFP04", "HFP05", "HFP06",
  "HFP07", "HFP08", "HFP09", "HFP10", "HFP11", "HFP12",
  "HO001", "HO002", "HO003", "HO004", "HO005", "HO006",
  "HO007", "HO008", "HO009", "HO010", "HO011", "HO012",
  "HO1M02", "HO1M04", "HO1M06", "HO1M08", "HO1M10", "HO1M12",
  "HO1P02", "HO1P04", "HO1P06", "HO1P08", "HO1P10", "HO1P12",
  "HO2M02", "HO2M04", "HO2M06", "HO2M08", "HO2M10", "HO2M12",
  "HO2P02", "HO2P04", "HO2P06", "HO2P08", "HO2P10", "HO2P12"
};

//
// constructors and destructor
//
HcalPedestalTuning::HcalPedestalTuning(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  HB_1Dpeds = fs->make<TH1D>("HB_1D" , "HB_1D" , 500 , 0 , 5 );
  HE_1Dpeds = fs->make<TH1D>("HE_1D" , "HE_1D" , 500 , 0 , 5 );
  HO_1Dpeds = fs->make<TH1D>("HO_1D" , "HO_1D" , 500 , 0 , 5 );
  HF_1Dpeds = fs->make<TH1D>("HF_1D" , "HF_1D" , 500 , 0 , 5 );
  
  HB_2Dpeds = fs->make<TH2D>("HB_2D" , "HB_2D" , 101 , -50 , 50, 75, 0, 75 );
  HE_2Dpeds = fs->make<TH2D>("HE_2D" , "HE_2D" , 101 , -50 , 50, 75, 0, 75  );
  HO_2Dpeds = fs->make<TH2D>("HO_2D" , "HO_2D" , 101 , -50 , 50, 75, 0, 75  );
  HF_2Dpeds = fs->make<TH2D>("HF_2D" , "HF_2D" , 101 , -50 , 50, 75, 0, 75  );
  
  digisName_          = iConfig.getParameter<string>("digisName");
  mapIOV_             = iConfig.getParameter<int>("mapIOV");
  tunePedestals_      = iConfig.getParameter<bool>("tunePedestals");
  pathToOldXMLBricks_ = 
    iConfig.getUntrackedParameter<string>("pathToOldXMLBricks");
  xmlBricksOutputDir_ = 
    iConfig.getUntrackedParameter<string>("xmlBricksOutputDir");
  uniformDACSettings_ = iConfig.getParameter<bool>("uniformDACSettings");
  initialDACSetting_  = iConfig.getParameter<int>("initialDACSetting");
  SiPMRBXes_          = iConfig.getParameter<vector <string> >("SiPMRBXes");
  //SiPMDACSetting_     = iConfig.getParameter<int>("SiPMDACSetting");
  targetPedestalHPD_  = iConfig.getParameter<double>("targetPedestalHPD");
  targetPedestalSiPM_ = iConfig.getParameter<double>("targetPedestalSiPM");
  dac2adc_            = iConfig.getParameter<double>("dac2adc");
  tagName_            = iConfig.getParameter<string>("tagName");
}


HcalPedestalTuning::~HcalPedestalTuning()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void HcalPedestalTuning::analyze(const edm::Event& iEvent,
				const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  if (jobStart) {
    // 2: HB Depth; 3: HE Depth
    HcalTopology hcalTopo(HcalTopologyMode::LHC, 2, 3);
    HcalLogicalMapGenerator lMapGen;
    lMap =  new HcalLogicalMap(lMapGen.createMap(&hcalTopo,mapIOV_));
        
    run = iEvent.id().run();
    stringstream tempstring;
    tempstring <<"Run"<<run;
    run_string = tempstring.str();
    // get the HCAL Electronics Map
    edm::ESHandle<HcalElectronicsMap> refEMap;
    iSetup.get<HcalElectronicsMapRcd>().get(refEMap);
    vector<HcalGenericDetId> listEMap = refEMap->allPrecisionId();
    // loop over the map and fill in the list of HCAL channels
    cout<<">> Filling in the list of HcalGenericDetIds...\n";
    for (vector<HcalGenericDetId>::const_iterator it = listEMap.begin(); 
	 it != listEMap.end(); it++) {
      ChannelData channelData;
      
      channelData.isHOX     = false;
      channelData.isCalib   = false;
      channelData.isSiPM    = false;
      
      if( it->isHcalCastorDetId() || 
	  it->isHcalTrigTowerDetId() || 
	  it->isHcalZDCDetId() ) continue;
      if( it->isHcalCalibDetId() ) {
	HcalCalibDetId hcalCalibDetId = HcalCalibDetId(it->rawId());
	if( hcalCalibDetId.calibFlavor() == HcalCalibDetId::HOCrosstalk ) 
	  channelData.isHOX = true;
	else 
	  channelData.isCalib = true;
      }
      
      channelData.rawId = it->rawId();
      channelData.frontEndId = 
	lMap->getHcalFrontEndId(DetId(channelData.rawId));
         
      // check if the channel is in one of the SiPM RBXes
      string theRBX = channelData.frontEndId.rbx();
      for(vector<string>::const_iterator it=SiPMRBXes_.begin(); 
	  it!=SiPMRBXes_.end(); it++) {
	if(theRBX==(*it)) channelData.isSiPM = true;
      }
      
      //if(channelData.isSiPM && !channelData.isCalib) 
      //	channelData.dacSetting = SiPMDACSetting_;
      //else 
      channelData.dacSetting = initialDACSetting_;
      for (int capId = 0; capId < 4; capId++) {
	channelData.capSum[capId] = 0;
	channelData.capEntries[capId] = 0;
      }
      
      channelList[channelData.rawId] = channelData;
      
      // fill in the map of rawIds
      HcalFrontEndId fakeFrontEndId = 
	HcalFrontEndId(channelData.frontEndId.rbx(), 
		       channelData.frontEndId.rm(), 
		       1, 1, 1, 
		       channelData.frontEndId.qieCard(), 
		       channelData.frontEndId.adc());
      rawIdMap[fakeFrontEndId] = channelData.rawId;
    }
    cout<<">> Total of "<<channelList.size()
	<<" HcalGenericDetIds filled in!\n";
    
    jobStart = false;
  }
  
  // get HCAL digis
  iEvent.getByLabel(digisName_,hbheDigis);
  iEvent.getByLabel(digisName_,hoDigis);
  iEvent.getByLabel(digisName_,hfDigis);
  iEvent.getByLabel(digisName_,calibDigis);
  // loop over HBHE digis
  for(HBHEDigiCollection::const_iterator d = hbheDigis->begin(); 
      d != hbheDigis->end(); d++) {
    uint32_t rawId = d->id().rawId();
    int size = d->size();
    if(digiStart) {
      ChannelData::digiSize = size;
      digiStart = false;
    }
    map<uint32_t, ChannelData>::iterator ch = channelList.find(rawId);
    for (int ts = 0; ts < size; ts++) {
      ch->second.capSum[d->sample(ts).capid()] += d->sample(ts).adc();
      ch->second.capEntries[d->sample(ts).capid()]++;
    }
  }
  // loop over HO digis
  for(HODigiCollection::const_iterator d = hoDigis->begin(); 
      d != hoDigis->end(); d++) {
    uint32_t rawId = d->id().rawId();
    int size = d->size();
    map<uint32_t, ChannelData>::iterator ch = channelList.find(rawId);
    for (int ts = 0; ts < size; ts++) {
      ch->second.capSum[d->sample(ts).capid()] += d->sample(ts).adc();
      ch->second.capEntries[d->sample(ts).capid()]++;
    }
  }
  // loop over HF digis
  for(HFDigiCollection::const_iterator d = hfDigis->begin(); 
      d != hfDigis->end(); d++) {
    uint32_t rawId = d->id().rawId();
    int size = d->size();
    map<uint32_t, ChannelData>::iterator ch = channelList.find(rawId);
    for (int ts = 0; ts < size; ts++) {
      ch->second.capSum[d->sample(ts).capid()] += d->sample(ts).adc();
      ch->second.capEntries[d->sample(ts).capid()]++;
    }
  }
  // loop over Calibration digis
  for(HcalCalibDigiCollection::const_iterator d = calibDigis->begin(); 
      d != calibDigis->end(); d++) {
    //if( !(d->id().calibFlavor() == HcalCalibDetId::HOCrosstalk) ) continue;
    uint32_t rawId = d->id().rawId();
    int size = d->size();
    map<uint32_t, ChannelData>::iterator ch = channelList.find(rawId);
    for (int ts = 0; ts < size; ts++) {
      ch->second.capSum[d->sample(ts).capid()] += d->sample(ts).adc();
      ch->second.capEntries[d->sample(ts).capid()]++;
    }
  }
  // increment event counter
  ChannelData::evtCounter++;

}


// ------------ method called once each job just before starting event loop  ------------
void HcalPedestalTuning::beginJob() {
  jobStart = true;
  digiStart = true;

  // Moved to analyze(): need EventSetup to create map in CMSSW>=6.1.0
  // create HCAL logical map
  //HcalLogicalMapGenerator lMapGen;
  //lMap = new HcalLogicalMap(lMapGen.createMap(3)); 
  ////IOV=3 -> May 6 2009 and onwards
  //    lMap->checkIdFunctions();
   
  for(vector<string>::const_iterator it = SiPMRBXes_.begin(); 
      it!=SiPMRBXes_.end(); it++) {
    bool RBXfound = false;
    for(int i=0; i<132; i++) {
      if( (*it)==RBXes[i] ) RBXfound = true;
    }
    if(!RBXfound) 
      cout<<">> ----> WARNING: '"
	  <<(*it)<<"' is not a valid RBX name!\n";
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void HcalPedestalTuning::endJob() {
  // read in pedestal settings
  if(!uniformDACSettings_) readSettings();
  // tune pedestals
  if(tunePedestals_) tunePedestals(targetPedestalHPD_, targetPedestalSiPM_);
  // write the output files
  writeOutput();
}

void HcalPedestalTuning::readSettings() {
  
  cout<<">> Reading in DAC settings...\n";
  
  int num_ch = 0;
  
  for(int i=0; i<132; i++) {
    
    TiXmlDocument xmlIn((pathToOldXMLBricks_ + "/" 
			 + RBXes[i] + ".xml").c_str());
    bool xmlLoaded = xmlIn.LoadFile();
    if (!xmlLoaded) {
      cout<<">> PROBLEM READING XML File:\n\t "
	  <<(pathToOldXMLBricks_ + "/" + RBXes[i] + ".xml").c_str()<<"\n";
      return;
    }
    TiXmlElement *brick, *data;
    brick = xmlIn.FirstChildElement("CFGBrick");
    if (brick) {
      data = brick->FirstChildElement("Data");
      while (data) {
	HcalFrontEndId fakeFrontEndId = 
	  HcalFrontEndId(RBXes[i],
			 atoi(data->Attribute("rm")),
			 1, 1, 1,
			 atoi(data->Attribute("card")),
			 atoi(data->Attribute("qie")));
	
	num_ch++;
	map<uint32_t, ChannelData>::iterator ch = 
	  channelList.find( (rawIdMap.find(fakeFrontEndId))->second );
	ch->second.dacSetting = atoi(data->Value()); // this was DAC
	data = data->NextSiblingElement("Data");
      }
    }
    else {
      cout<<">> PROBLEM: CFGBrick element not foundi in file:\n\t "
	  <<(pathToOldXMLBricks_ + "/" + RBXes[i] + ".xml").c_str()<<"\n";
      return;
    }

  }
  cout<<">> DAC settings for "<<num_ch<<" HCAL channels were read in!\n";

  return;
}

void HcalPedestalTuning::tunePedestals(const double targetPedestalHPD, 
				      const double targetPedestalSiPM) {
  
  cout<<">> Tuning pedestals...\n";

  double diff;
  int oldSetting;
  int newSetting;

  map<uint32_t, ChannelData>::iterator ch;
  for (ch = channelList.begin(); ch != channelList.end(); ch++) {
    if( (ch->second.getStability()==0.0)) continue; // || (ch->second.isSiPM && !ch->second.isCalib) ) continue;
   
    // check that there are entries for all 4 capId's
    int entries = 0;
    for(int i=0; i<4; i++){
      if(ch->second.capEntries[i] > 0) entries++;
    }
    // this case should never happen but it's here for completeness
    if( entries > 0 && entries < 4) 
      cout<<">> ----> WARNING: Channel problem found in RBX "
	  <<ch->second.frontEndId.rbx()<<" RM "
	  <<ch->second.frontEndId.rm()<<" QIE "
	  <<ch->second.frontEndId.qieCard()<<"!\n"; 

    double channelPedestal = 
      ch->second.getChannelPedestal(); // get the pedestal
    double targetPedestal;

    // decide whether it is a SiPM or an HPD; set target accordingly
    if (ch->second.isSiPM) targetPedestal = targetPedestalSiPM;
    if (!ch->second.isSiPM) targetPedestal = targetPedestalHPD;

    // tune the pedestal
    oldSetting = ch->second.dacSetting;
    diff = (targetPedestal - channelPedestal)/dac2adc_; 
    if (diff < 0) diff -= 0.5;
    if (diff > 0) diff +=  0.5;
    newSetting = oldSetting + (int)diff;

    if (newSetting > 7) newSetting = 7;
    if (newSetting < 0) newSetting = 0;

    if (newSetting != oldSetting) 
      tunedChannelsByRBX[ch->second.frontEndId.rbx()].push_back(ch->second);
    
    if (newSetting != oldSetting) ch->second.dacSetting = newSetting;
    
  }
  cout << ">> Pedestal tuning completed!\n";
  
  return;
}

void HcalPedestalTuning::writeOutput() {
  
  cout<< ">> Writing output XML bricks...\n";

  // set the output directory for new XML bricks
  string outputDir;
  if( xmlBricksOutputDir_ == pathToOldXMLBricks_ ) {
    outputDir = xmlBricksOutputDir_ + "/" + tagName_ + "_bricks/";
    boost::filesystem::create_directories
      (boost::filesystem::path(outputDir));

    cout<<">> ----> WARNING: 'xmlBricksOutputDir_' and "
	<<"'pathToOldXMLBricks_' point to the same location.\n"; 
    cout<<">> ----> To avoid overwriting the old XML bricks"
	<<" with the new ones, new XML bricks were written to '"
	<<outputDir.c_str()<<"'!\n";
  } else {
    outputDir = xmlBricksOutputDir_ + "/";
    boost::filesystem::create_directories
      (boost::filesystem::path(outputDir));
  }

  // get the date
  time_t rawtime;
  struct tm *timeinfo;
  char date[64];

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  strftime (date,64,"%d-%m-%y",timeinfo);
  // fill the list of HCAL channels organized by RBX
  map<uint32_t, ChannelData>::const_iterator c;
  for (c = channelList.begin(); c != channelList.end(); c++) {

    channelListByRBX[c->second.frontEndId.rbx()].push_back(c->second);
  }
  ofstream fOutAll((outputDir + run_string + "_all_channels.txt").c_str()); 
  // output file with the list of all channels
  char buffer[1024];
  // sort channels and write out the .xml bricks and the list of all channels

  map<string, vector<ChannelData> >::iterator it;
  for ( it = channelListByRBX.begin(); it != channelListByRBX.end(); it++) {
   
    sprintf(buffer, 
	    "# %8s %8s %8s %8s %8s %8s %8s %8s %12s %8s "
	    "%10s %10s %10s %10s %10s\n", 
	    "RBX", "RM", "QIE", "ADC", "subdet", "eta", "phi", "depth",
	    "stability", "DAC", "ped0", "ped1", "ped2", "ped3", "PED");
    fOutAll<<buffer;

    sort(it->second.begin(), it->second.end(), compareRBXInfo());
  
    string theRBX = it->first;
      
    TiXmlDocument xmlOut;
    xmlOut.LinkEndChild(new TiXmlDeclaration("1.0","",""));
    TiXmlElement *brick = new TiXmlElement("CFGBrick");
    TiXmlElement *rbx = new TiXmlElement("Parameter");
    rbx->SetAttribute("name","RBX");
    rbx->SetAttribute("type","string");
    rbx->LinkEndChild(new TiXmlText(theRBX));
    brick->LinkEndChild(rbx);
    TiXmlElement *info = new TiXmlElement("Parameter");
    info->SetAttribute("name","INFOTYPE");
    info->SetAttribute("type","string");
    info->LinkEndChild(new TiXmlText("PEDESTAL"));
    brick->LinkEndChild(info);
    TiXmlElement *tag = new TiXmlElement("Parameter");
    tag->SetAttribute("name","CREATIONTAG");
    tag->SetAttribute("type","string");
    tag->LinkEndChild(new TiXmlText(tagName_));
    brick->LinkEndChild(tag);
    TiXmlElement *stamp = new TiXmlElement("Parameter");
    stamp->SetAttribute("name","CREATIONSTAMP");
    stamp->SetAttribute("type","string");
    stamp->LinkEndChild(new TiXmlText(date));
    brick->LinkEndChild(tag);
    brick->LinkEndChild(stamp);
    xmlOut.LinkEndChild(brick);

    for (vector<ChannelData>::const_iterator ch = it->second.begin();
	 ch != it->second.end(); ch++) {

      //if(ch->isCalib) continue;
               
      int RM   = ch->frontEndId.rm();
      int CARD = ch->frontEndId.qieCard();
      int QIE  = ch->frontEndId.adc();
      string subDet;
      int ieta, iphi, depth;
      double stability = ch->getStability();
      int DAC     = ch->dacSetting;
      double ped0 = ch->getCapPedestal(0);
      double ped1 = ch->getCapPedestal(1);
      double ped2 = ch->getCapPedestal(2);
      double ped3 = ch->getCapPedestal(3);
      double PED  = ch->getChannelPedestal();

      if(ch->isHOX) {
	HcalCalibDetId hcalCalibDetId = HcalCalibDetId(ch->rawId);
	subDet = "HOX";
	ieta   = hcalCalibDetId.ieta();
	iphi   = hcalCalibDetId.iphi();
	depth  = 4;
      } else if(ch->isCalib) {
	HcalCalibDetId hcalCalibDetId = HcalCalibDetId(ch->rawId);
	subDet = "C_" + subdet(hcalCalibDetId.hcalSubdet());
	ieta   = hcalCalibDetId.ieta();
	iphi   = hcalCalibDetId.iphi();
	depth  = 0;
      } else {
	HcalDetId hcalDetId = HcalDetId(ch->rawId);
	subDet = subdet(hcalDetId.subdet());
	ieta   = hcalDetId.ieta();
	iphi   = hcalDetId.iphi();
	depth  = hcalDetId.depth();

	if (subDet=="HB") {
	  HB_1Dpeds->Fill(PED);
	  HB_2Dpeds->SetBinContent(ieta+51,iphi,PED);
	}
	if (subDet=="HE") {
	  HE_1Dpeds->Fill(PED);
	  HE_2Dpeds->SetBinContent(ieta+51,iphi,PED);
	}
	if (subDet=="HO") {
	  HO_1Dpeds->Fill(PED);
	  HO_2Dpeds->SetBinContent(ieta+51,iphi,PED);
	}
	if (subDet=="HF") {
	  HF_1Dpeds->Fill(PED);
	  HF_2Dpeds->SetBinContent(ieta+51,iphi,PED);
	}
      }
	 

      sprintf(buffer,
	      "  %8s %8i %8i %8i %8s %8i %8i %8i"
	      " %12.4f %8i %10.5f %10.5f %10.5f %10.5f %10.5f\n",
	      theRBX.c_str(), RM, CARD, QIE, subDet.c_str(),
	      ieta, iphi, depth, stability, DAC, ped0, ped1, ped2, ped3, PED);
      fOutAll<<buffer;

      stringstream ss;
      ss << hex << DAC;
      TiXmlElement *data = new TiXmlElement("Data");
      data->SetAttribute("elements","1");
      data->SetAttribute("encoding","hex");
      data->SetAttribute("rm",RM);
      data->SetAttribute("card",CARD);
      data->SetAttribute("qie",QIE);
      data->LinkEndChild(new TiXmlText(ss.str().c_str()));
      brick->LinkEndChild(data);
    }
    // close output file
    xmlOut.SaveFile((outputDir + theRBX + ".xml").c_str());
  }
  fOutAll<< endl;
  // close output file
  fOutAll.close();
   
  cout<< ">> XML bricks written out!\n";
   
  // write the list of missing and unstable channels
  writeUnstable(outputDir);
  // write the list of untuned channels
  if(tunePedestals_) writeTuned(outputDir);
   
  return;
}

void HcalPedestalTuning::writeUnstable(const string& fOutputDir) {
   // fill in the list of missing and unstable HCAL channels
   map<uint32_t, ChannelData>::const_iterator ch;
   for (ch = channelList.begin(); ch != channelList.end(); ch++) {
      if( ch->second.getStability() < 1. ) 
	unstableChannelsByRBX[ch->second.frontEndId.rbx()].
	  push_back(ch->second);
   }
   // sort and write out missing channels
   ofstream fOut((fOutputDir + run_string + 
		  "_missing_&_unstable_channels.txt").c_str());
   char buffer[1024];
   sprintf(buffer,
	   "# %8s %8s %8s %8s %8s %8s %8s %8s %12s %8s "
	   "%10s %10s %10s %10s %10s\n",
	   "RBX", "RM", "QIE", "ADC", "subdet", "eta", "phi", "depth",
	   "stability", "DAC", "ped0", "ped1", "ped2", "ped3", "PED");
   fOut<<buffer;
   
   map<string, vector<ChannelData> >::iterator it;
   for ( it = unstableChannelsByRBX.begin(); 
	 it != unstableChannelsByRBX.end(); it++) {
   
      sort(it->second.begin(), it->second.end(), compareRBXInfo());
  
      string theRBX = it->first;

      for (vector<ChannelData>::const_iterator ch = it->second.begin(); 
	   ch != it->second.end(); ch++) {
         
         int RM   = ch->frontEndId.rm();
         int CARD = ch->frontEndId.qieCard();
         int QIE  = ch->frontEndId.adc();
         string subDet;
         int ieta, iphi, depth;
         double stability = ch->getStability();
         int DAC     = ch->dacSetting;
         double ped0 = ch->getCapPedestal(0);
         double ped1 = ch->getCapPedestal(1);
         double ped2 = ch->getCapPedestal(2);
         double ped3 = ch->getCapPedestal(3);
         double PED  = ch->getChannelPedestal();

         if(ch->isHOX) {
            HcalCalibDetId hcalCalibDetId = HcalCalibDetId(ch->rawId);
            subDet = "HOX";
            ieta   = hcalCalibDetId.ieta();
            iphi   = hcalCalibDetId.iphi();
            depth  = 4;
         } else if(ch->isCalib) {
            HcalCalibDetId hcalCalibDetId = HcalCalibDetId(ch->rawId);
            subDet = "C_" + subdet(hcalCalibDetId.hcalSubdet());
            ieta   = hcalCalibDetId.ieta();
            iphi   = hcalCalibDetId.iphi();
            depth  = 0;
         } else {
            HcalDetId hcalDetId = HcalDetId(ch->rawId);
            subDet = subdet(hcalDetId.subdet());
            ieta   = hcalDetId.ieta();
            iphi   = hcalDetId.iphi();
            depth  = hcalDetId.depth();
         }
         
         sprintf(buffer,
		 "  %8s %8i %8i %8i %8s %8i %8i %8i %12.4f %8i"
		 " %10.5f %10.5f %10.5f %10.5f %10.5f\n",
		 theRBX.c_str(), RM, CARD, QIE, subDet.c_str(),
		 ieta, iphi, depth, stability,
		 DAC, ped0, ped1, ped2, ped3, PED);
         fOut<<buffer;
      }
   }
   fOut<<endl;
   // close output file
   fOut.close();
   
   cout<<">> For the missing HCAL channels listed in "
       <<("'" + run_string + "_missing_&_unstable_channels.txt'").c_str()
       <<", DAC settings were left unchanged!\n";
   
   return;
}

void HcalPedestalTuning::writeTuned(const string& fOutputDir) {
   // sort and write out untuned channels
   ofstream fOut((fOutputDir + run_string + "_tuned_channels.txt").c_str());
   char buffer[1024];
   sprintf(buffer,
	   "# %8s %8s %8s %8s %8s %8s %8s %8s %12s %8s "
	   "%10s %10s %10s %10s %10s\n",
	   "RBX", "RM", "QIE", "ADC", "subdet", "eta", "phi", "depth",
	   "stability", "DAC", "ped0", "ped1", "ped2", "ped3", "PED");
   fOut<<buffer;
   
   map<string, vector<ChannelData> >::iterator it;
   for ( it = tunedChannelsByRBX.begin();
	 it != tunedChannelsByRBX.end(); it++) {
   
      sort(it->second.begin(), it->second.end(), compareRBXInfo());
  
      string theRBX = it->first;

      for (vector<ChannelData>::const_iterator ch = it->second.begin();
	   ch != it->second.end(); ch++) {
    
         int RM   = ch->frontEndId.rm();
         int CARD = ch->frontEndId.qieCard();
         int QIE  = ch->frontEndId.adc();
         string subDet;
         int ieta, iphi, depth;
         double stability = ch->getStability();
         int DAC     = ch->dacSetting;
         double ped0 = ch->getCapPedestal(0);
         double ped1 = ch->getCapPedestal(1);
         double ped2 = ch->getCapPedestal(2);
         double ped3 = ch->getCapPedestal(3);
         double PED  = ch->getChannelPedestal();

         if(ch->isHOX) {
            HcalCalibDetId hcalCalibDetId = HcalCalibDetId(ch->rawId);
            subDet = "HOX";
            ieta   = hcalCalibDetId.ieta();
            iphi   = hcalCalibDetId.iphi();
            depth  = 4;
         } else if(ch->isCalib) {
	    HcalCalibDetId hcalCalibDetId = HcalCalibDetId(ch->rawId);
	    subDet = "C_" + subdet(hcalCalibDetId.hcalSubdet());
	    ieta   = hcalCalibDetId.ieta();
	    iphi   = hcalCalibDetId.iphi();
	    depth  = 0;
	 }
           else {
            HcalDetId hcalDetId = HcalDetId(ch->rawId);
            subDet = subdet(hcalDetId.subdet());
            ieta   = hcalDetId.ieta();
            iphi   = hcalDetId.iphi();
            depth  = hcalDetId.depth();
	   }

         if (ch->isCalib == true) continue;
         else sprintf(buffer, "  %8s %8i %8i %8i %8s %8i %8i %8i"
		      " %12.4f %8i %10.5f %10.5f %10.5f %10.5f %10.5f\n",
		      theRBX.c_str(), RM, CARD, QIE, subDet.c_str(),
		      ieta, iphi, depth, stability,
		      DAC, ped0, ped1, ped2, ped3, PED);
         fOut<<buffer;
      }
   }
   fOut<<endl;
   // close output file
   fOut.close();
   
   cout<<">> List of the HCAL channels that required tuning was written to "
       <<("'" + run_string + "_tuned_channels.txt'").c_str()<<"!\n";
   
   return;
}

double HcalPedestalTuning::ChannelData::getCapPedestal(int capId) const {
   if (capEntries[capId] == 0) return 0;
   else return capSum[capId]/capEntries[capId];
}

double HcalPedestalTuning::ChannelData::getChannelPedestal() const {
   int entries = capEntries[0]+capEntries[1]+capEntries[2]+capEntries[3];
   if (entries == 0) return 0;
   else return (capSum[0]+capSum[1]+capSum[2]+capSum[3])/entries;
}

double HcalPedestalTuning::ChannelData::getStability() const {
   return ((double)(capEntries[0]+capEntries[1]+
		    capEntries[2]+capEntries[3]))/(digiSize*evtCounter);
}

HcalSubdetector HcalPedestalTuning::subdet(const string& s) {
  HcalSubdetector retval = HcalEmpty;
  if      (s=="HB") retval = HcalBarrel;
  else if (s=="HE") retval = HcalEndcap;
  else if (s=="HO") retval = HcalOuter;
  else if (s=="HF") retval = HcalForward;

  return retval;
}

string HcalPedestalTuning::subdet(const HcalSubdetector& s) {
  string retval = "";
  if      (s==HcalBarrel)  retval = "HB";
  else if (s==HcalEndcap)  retval = "HE";
  else if (s==HcalOuter)   retval = "HO";
  else if (s==HcalForward) retval = "HF";

  return retval;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalPedestalTuning);
