/////////////////////////////////////////////////////////////////////////
// Class:       BeamEvent
// Plugin Type: producer (art v2_08_03)
// File:        BeamEvent_module.cc
//
// Written by Jake Calcutt (calcuttj@msu.edu)
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "IFBeam_service.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "dune/DuneObj/ProtoDUNEBeamSpill.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNETimeStamp.h"
#include <bitset>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <limits>

#include "TTree.h"
#include "TH2F.h"
#include "TVectorD.h"
#include "TPolyMarker.h"

namespace proto {
  class BeamEvent;
}

enum tofChan{
  k1A2A,
  k1B2A,
  k1A2B,
  k1B2B
};

typedef std::numeric_limits< double > dbl;

class proto::BeamEvent : public art::EDProducer {
public:
  explicit BeamEvent(fhicl::ParameterSet const & p);

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamEvent(BeamEvent const &) = delete;
  BeamEvent(BeamEvent &&) = delete;
  BeamEvent & operator = (BeamEvent const &) = delete;
  BeamEvent & operator = (BeamEvent &&) = delete;

  // Required functions.
  //void reconfigure(fhicl::ParameterSet const & p);
  void reset();
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

  uint64_t joinHighLow(double,double);

  TVector3 ConvertProfCoordinates(double x, double y, double z, double zOffset);
  void BeamMonitorBasisVectors(); 
  bool rotated = false;
  void RotateMonitorVector(TVector3 &vec); 

  void GetRawDecoderInfo(art::Event &);

  void TimeIn(art::Event &, uint64_t);
  void GetSpillInfo(art::Event &);
  void MatchBeamToTPC();
  void MatchS11ToGen();
  void SetBeamEvent(); 
  void MakeTrack(size_t);
  void MomentumSpec(size_t);
  double MomentumCosTheta(double, double, double);

  double GetPosition(std::string, int);
  TVector3 ProjectToTPC(TVector3, TVector3);
  double GetPairedPosition(std::string, size_t);
 
  void  InitXBPFInfo(beam::ProtoDUNEBeamSpill *);
  void  parseGeneralXBPF(std::string, uint64_t, size_t);
  void  parseXBPF(uint64_t);

  void  parseXTOF(uint64_t);
  void  parseXCETDB(uint64_t);


  void getS11Info(uint64_t);

  std::vector<double> FetchAndReport(long long, std::string, std::unique_ptr<ifbeam_ns::BeamFolder>& );
   
private:
  
  TTree * fOutTree;
  TH1F  * fFullMomentum;
  TH1F  * fCutMomentum;

  TTree * fGenTrigTree;
  TTree * fXTOF1ATree;
  TTree * fXTOF1BTree;
  TTree * fXTOF2ATree;
  TTree * fXTOF2BTree;


  double fGenTrigFrac;
  double fGenTrigSec;
  double fGenTrigCoarse; 
  double fXTOF1AFrac;
  double fXTOF1ACoarse; 
  double fXTOF1BFrac;
  double fXTOF1BCoarse; 
  double fXTOF2AFrac;
  double fXTOF2ACoarse; 
  double fXTOF2BFrac;
  double fXTOF2BCoarse; 

  double fXTOF1ASec;
  double fXTOF1BSec;
  double fXTOF2ASec;
  double fXTOF2BSec;
  std::vector<double> diff2A;
  std::vector<double> diff2B;
 

  long long int eventTime;
  double SpillStart;
  ULong_t SpillStart_alt;
  bool SpillStartValid;
  bool acqStampValid;
  double PrevStart=-99999.;
  double SpillEnd;
  double SpillOffset;
  double ActiveTriggerTime;
  long long RDTSTime;
  double RDTSTimeSec;
  double PrevRDTSTimeSec=-99999.;  
  double RDTSTimeNano; 

  long long valid_fetch_time; 
  long long spill_valid_fetch_time;

  double s11Sec, s11Nano;

  int RDTSTrigger;

  double acqTime;
  double acqStampMBPL;

  int C1DB;
  int C2DB;


  int eventNum;
  int runNum;
  int subRunNum;
  double CKov1Pressure;
  double CKov2Pressure;
  int CKov1Counts;
  int CKov2Counts;
  
  TVector3 fBMBasisX = TVector3(1.,0.,0.);
  TVector3 fBMBasisY = TVector3(0.,1.,0.);
  TVector3 fBMBasisZ = TVector3(0.,0.,1.);

  double  fTimeWindow;
  std::string fBundleName;
  std::string fXCETBundleName;
  std::string fOutputLabel;
  std::string fURLStr;
  double fBFEpsilon, fXCETEpsilon;
  double fXCETFetchShift;
  int fIFBeamDebug;
  uint64_t fFixedTime;

  std::string firstUpstreamName;
  std::string secondUpstreamName;
  std::string firstDownstreamName;
  std::string secondDownstreamName;

  std::string firstBPROF1;
  std::string secondBPROF1;
  std::string BPROF2;
  std::string BPROF3;
  double      fBeamBend;
  double L1, L2, L3;
  double fBProf1Shift;
  double fBProf2Shift;
  double fBProf3Shift;

  std::vector< std::string > fDevices;
  std::map<std::string, std::string > fDeviceTypes;
  std::map< std::string, double > fFiberDimension;


  // Names of the CERN Beam Devices
  // Names have the form of a prefix + device
  
  // Prefixes for different devices classes:
  // Beam Positions (BPF)
  // Time of flights (TOF)
  // and Cerenkov counters (CET)
  std::string fXBPFPrefix;
  std::string fXTOFPrefix;
  std::string fXCETPrefix;

  // Device Names for each device

  // Time of Flights
  std::string fTOF1;
  std::string fTOF1A, fTOF1B;
  std::string fTOF2;
  std::string fTOF2A, fTOF2B;

  double fTOFCalAA, fTOFCalBA, fTOFCalAB, fTOFCalBB;

  // Cerenkovs
  std::string fCKov1;
  std::string fCKov2;
  std::string fXCET1;
  std::string fXCET2;

  double fRotateMonitorXZ;
  double fRotateMonitorYZ;
  double fRotateMonitorYX;

  double fFirstTrackingProfZ;
  double fSecondTrackingProfZ;
  double fNP04FrontZ;  
  double fBeamX, fBeamY, fBeamZ;

  bool   fForceNewFetch;
  bool   fXCETDebug;
  bool   fMatchTime;
  bool   fForceRead;
  bool   fForceMatchS11;
  
  double fFillCacheUp, fFillCacheDown;

  bool   fPrintDebug;
  bool   fSaveOutTree;
  bool   fDebugTOFs;
  bool   fDebugMomentum;

  bool   gotGeneralTrigger;
  bool   gotTOFs;
  bool   gotCurrent;

  double fTimingCalibration;
  double fCalibrationTolerance;
  double fOffsetTAI;

  double fS11DiffUpper; 
  double fS11DiffLower; 
  double fRDTSToS11Upper; 
  double fRDTSToS11Lower; 

  int    fOffsetCTBtoRDTS;
  int    fToleranceCTBtoRDTS;

  double fDownstreamToGenTrig;
  double fUpstreamToDownstream;

  beam::ProtoDUNEBeamEvent * beamevt;
  beam::ProtoDUNEBeamEvent prev_beamevt;

  beam::ProtoDUNEBeamSpill * beamspill;
  beam::ProtoDUNEBeamSpill prev_beamspill;

  uint64_t cache_start = 0;
  uint64_t cache_end   = 0;

  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_xcet;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;

  art::Handle< std::vector<raw::RDTimeStamp> > RDTimeStampHandle;

//  double L1=1.980, L2=1.69472, L3=2.11666;
  double magnetLen, magnetField;
  std::vector<double> current;

  //Hardware Parameters for magnetic field stuff
  double mag_P1 = 5.82044830e-3;
  // unused double mag_P2 = 0.;
  double mag_P3 = -4.68880000e-6;
  double mag_P4 = 324.573967;

};

// Constructor
proto::BeamEvent::BeamEvent(fhicl::ParameterSet const & p)
: EDProducer(p)
{
  // Declare products this module will provide
  produces<std::vector<beam::ProtoDUNEBeamEvent>>();  

  // Implementation of optional member function here.
  fBundleName      = p.get<std::string>("BundleName");
  fXCETBundleName  = p.get<std::string>("XCETBundleName");
  fOutputLabel = p.get<std::string>("OutputLabel");
  fURLStr      = p.get<std::string>("URLStr");
  fBFEpsilon   = p.get<double>("BFEpsilon");
  fXCETEpsilon   = p.get<double>("XCETEpsilon");
  fXCETFetchShift   = p.get<double>("XCETFetchShift");
  fIFBeamDebug = p.get<int>("IFBeamDebug");
  fTimeWindow  = p.get<double>("TimeWindow");
  fFixedTime   = p.get<uint64_t>("FixedTime");

  fDevices     = p.get< std::vector< std::string > >("Devices");    

  //For Tracking/////
  firstUpstreamName    = p.get< std::string >("FirstUpstream");
  secondUpstreamName   = p.get< std::string >("SecondUpstream");
  firstDownstreamName  = p.get< std::string >("FirstDownstream");
  secondDownstreamName = p.get< std::string >("SecondDownstream");
  ///////////////////

  //For Momentum Spectrometry////
  firstBPROF1  = p.get< std::string >("FirstBPROF1");
  secondBPROF1 = p.get< std::string >("SecondBPROF1");
  BPROF2       = p.get< std::string >("BPROF2");
  BPROF3       = p.get< std::string >("BPROF3");
  fBeamBend    = p.get< double >("BeamBend");
  L1           = p.get< double >("L1"); 
  L2           = p.get< double >("L2"); 
  L3           = p.get< double >("L3"); 
  fBProf1Shift           = p.get< double >("BProf1Shift"); 
  fBProf2Shift           = p.get< double >("BProf2Shift"); 
  fBProf3Shift           = p.get< double >("BProf3Shift"); 

  magnetLen    = 1.;//(m)
  magnetField  = 1000.;//()
  ///////////////////////////////

  std::vector< std::pair<std::string, std::string> >  tempTypes = p.get<std::vector< std::pair<std::string, std::string> >>("DeviceTypes");
  fDeviceTypes  = std::map<std::string, std::string>(tempTypes.begin(), tempTypes.end() );

  //Deminsion of Fibers
  std::vector< std::pair<std::string, double> > tempFiberDims = p.get< std::vector<std::pair<std::string,double> > >("Dimension");
  fFiberDimension = std::map<std::string, double>(tempFiberDims.begin(), tempFiberDims.end());

  
  //XTOF devices 
  fTOF1 = p.get< std::string >("TOF1");
  fTOF2 = p.get< std::string >("TOF2");

  fTOF1A = fTOF1 + "A";
  fTOF1B = fTOF1 + "B";
  fTOF2A = fTOF2 + "A";
  fTOF2B = fTOF2 + "B";

  fCKov1 = p.get< std::string >("CKov1");
  fCKov2 = p.get< std::string >("CKov2");
  fXCET1 = p.get< std::string >("XCET1");
  fXCET2 = p.get< std::string >("XCET2");

  fFillCacheUp   = p.get< double >("FillCacheUp");
  fFillCacheDown = p.get< double >("FillCacheDown");

  fXBPFPrefix      = p.get<std::string>("XBPFPrefix");
  fXTOFPrefix      = p.get<std::string>("XTOFPrefix");
  fXCETPrefix      = p.get<std::string>("XCETPrefix");


  //New parameters to match Leigh's
  fRotateMonitorXZ = p.get<double>("RotateMonitorXZ");
  fRotateMonitorYZ = p.get<double>("RotateMonitorYZ");
  fRotateMonitorYX = p.get<double>("RotateMonitorYX");

  fFirstTrackingProfZ  = p.get<double>("FirstTrackingProfZ");
  fSecondTrackingProfZ = p.get<double>("SecondTrackingProfZ");
  fNP04FrontZ          = p.get<double>("NP04FrontZ");  

  fBeamX               = p.get<double>("BeamX");
  fBeamY               = p.get<double>("BeamY");
  fBeamZ               = p.get<double>("BeamZ");

  fForceNewFetch       = p.get<bool>("ForceNewFetch");
  fMatchTime           = p.get<bool>("MatchTime");
  fForceRead           = p.get<bool>("ForceRead");
  fForceMatchS11       = p.get<bool>("ForceMatchS11");
  fXCETDebug           = p.get<bool>("XCETDebug");


  fTimingCalibration      = p.get<double>("TimingCalibration");
  fCalibrationTolerance   = p.get<double>("CalibrationTolerance");
  fOffsetTAI              = p.get<double>("OffsetTAI");

  fS11DiffUpper           = p.get<double>("S11DiffUpper");
  fS11DiffLower           = p.get<double>("S11DiffLower");

  fRDTSToS11Upper         = p.get<double>("RDTSToS11Upper");
  fRDTSToS11Lower         = p.get<double>("RDTSToS11Lower");

  fOffsetCTBtoRDTS        = p.get<int>("OffsetCTBtoRDTS");
  fToleranceCTBtoRDTS     = p.get<int>("ToleranceCTBtoRDTS");

  fDownstreamToGenTrig    = p.get<double>("DownstreamToGenTrig");
  fUpstreamToDownstream   = p.get<double>("UpstreamToDownstream");

  fPrintDebug             = p.get<bool>("PrintDebug");
  fSaveOutTree            = p.get<bool>("SaveOutTree");
  fDebugMomentum          = p.get<bool>("DebugMomentum");
  fDebugTOFs              = p.get<bool>("DebugTOFs");

  fTOFCalAA               = p.get<double>("TOFCalAA");
  fTOFCalBA               = p.get<double>("TOFCalBA");
  fTOFCalAB               = p.get<double>("TOFCalAB");
  fTOFCalBB               = p.get<double>("TOFCalBB");


}
// END Constructor
////////////////////////

std::vector<double> proto::BeamEvent::FetchAndReport(long long time, std::string name, std::unique_ptr<ifbeam_ns::BeamFolder>& the_folder){

  //Note! Sometimes this won't retrieve the data from the database. We'll need
  //      catch the exception outside of this and handle it according to whichever
  //      device we're trying to retrieve from.

  std::vector<double> theResult;
  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << "Trying to grab from folder: " << name << "\n";
    MF_LOG_INFO("BeamEvent") << "At Time: " << time << "\n";    
  }

  theResult = the_folder->GetNamedVector(time, name);

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Successfully fetched " << time << "\n";

  return theResult;
}



//Gets the Timing and CTB raw-decoded info.
//Finds the triggers, and looks for a valid trigger word
//(i.e. coming from beam)
//
//Returns the timestamp of the high level trigger.
void proto::BeamEvent::GetRawDecoderInfo(art::Event & e){

  if( fPrintDebug ){ 
    MF_LOG_INFO("BeamEvent") << "\n";
    MF_LOG_INFO("BeamEvent") << "Getting Raw Decoder Info" << "\n";
  }

  e.getByLabel("timingrawdecoder","daq",RDTimeStampHandle);

  for (auto const & RDTS : *RDTimeStampHandle){

    uint64_t high = RDTS.GetTimeStamp_High();
    uint64_t low  = RDTS.GetTimeStamp_Low();
    high = high << 32; 
    uint64_t joined = (high | low);

    RDTSTime = joined;
    RDTSTrigger = RDTS.GetFlags();

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << "High: " << RDTS.GetTimeStamp_High() << "\n";
      MF_LOG_INFO("BeamEvent") << "Low: " << RDTS.GetTimeStamp_Low() << "\n"; 
      MF_LOG_INFO("BeamEvent") << "Raw Decoder Timestamp: " << joined << "\n";
      MF_LOG_INFO("BeamEvent") << "Trigger: " << RDTSTrigger << "\n"; 
    }

    //Separates seconds portion of the ticks 
    //From the nanoseconds
    long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
    RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
    long long RDTSTickNano = RDTSTime - RDTSTickSec;
  
    //Units are 20 nanoseconds ticks
    RDTSTimeSec  = 20.e-9 * RDTSTickSec;
    RDTSTimeNano = 20.    * RDTSTickNano;

  }

}

void proto::BeamEvent::GetSpillInfo(art::Event & e){

    auto const PDTStampHandle = e.getValidHandle< std::vector< dune::ProtoDUNETimeStamp > > ("timingrawdecoder:daq");  
    std::vector< dune::ProtoDUNETimeStamp > PDTStampVec(*PDTStampHandle); 
    auto PDTStamp = PDTStampVec[0];            
          
    UInt_t    ver        = PDTStamp.getVersion();   
    double RunStart   = 2.e-8*PDTStamp.getLastRunStart();            
    SpillStart = 2.e-8*PDTStamp.getLastSpillStart();          
    SpillEnd   = 2.e-8*PDTStamp.getLastSpillEnd();            
    SpillStart_alt = PDTStamp.getLastSpillStart();


    if( 0ul == ~(SpillStart_alt) ) SpillStartValid = false;
    else SpillStartValid = true;

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << PDTStamp.getLastSpillStart() << " " << SpillStart_alt << " " << SpillStartValid << "\n";
      MF_LOG_INFO("BeamEvent") << "Version:     " << ver        << "\n";      
      MF_LOG_INFO("BeamEvent") << "Run:         " << RunStart   << "\n";      
      MF_LOG_INFO("BeamEvent") << "Spill Start: " << std::setw(20) << SpillStart <<  "\n";      
      MF_LOG_INFO("BeamEvent") << "Spill End:   " << SpillEnd   <<  "\n";      
    }
          
          
}

void proto::BeamEvent::TimeIn(art::Event & e, uint64_t time){

    /////Now look at the acqStamp coming out of IFBeam
    try{
      std::vector<double> acqStamp = FetchAndReport(time, "dip/acc/NORTH/NP04/POW/MBPL022699:acqStamp[]", bfp); 

      if( acqStamp[0] < 300000000.0 ){
        if( fPrintDebug ){
          MF_LOG_INFO("BeamEvent") << "Warning: MBPL Spill Start is low " << acqStamp[0] 
                                << "\nWill need to time in with S11\n";
        }
        acqStampValid = false;                              
      }
      else{ acqStampValid = true; }

      acqStampMBPL   = 1.e-9 * joinHighLow(acqStamp[0],   acqStamp[1]); 
      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "MBPL: " << acqStampMBPL << "\n";

      //Assign the calibration offset
      SpillOffset = SpillStart - acqStampMBPL;

    }
    catch( std::exception e){
      acqStampValid = false;
      MF_LOG_WARNING("BeamEvent") << "Could not get Spill time to time in\n";
    }
    
}

void proto::BeamEvent::MatchS11ToGen(){
  
  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Matching S11 To Gen" << "\n";

  for(size_t iT = 0; iT < beamspill->GetNT0(); ++iT){
    double GenTrigSec  = beamspill->GetT0(iT).first;
    double GenTrigNano = beamspill->GetT0(iT).second;

    double diff = GenTrigSec - s11Sec;
    diff += 1.e-9*(GenTrigNano - s11Nano);

    if( fS11DiffLower < diff && diff < fS11DiffUpper ){

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "Found matching S11 and GenTrig!" << "\n";
        MF_LOG_INFO("BeamEvent") << "diff: " << diff << "\n";
      }

      beamspill->SetActiveTrigger( iT ); 
      beamevt->SetActiveTrigger( iT );
      beamevt->SetT0( beamspill->GetT0( iT ) );
      return;
    }
  }

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Could not find matching time " << "\n";

  beamspill->SetUnmatched();
}


void proto::BeamEvent::MatchBeamToTPC(){

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Matching in time between Beamline and TPC!!!" << "\n"; 
  for(size_t iT = 0; iT < beamspill->GetNT0(); ++iT){
    
    double GenTrigSec  = beamspill->GetT0(iT).first;
    double GenTrigNano = beamspill->GetT0(iT).second;

    //Separates seconds portion of the ticks 
    //From the nanoseconds
    long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
    RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
    long long RDTSTickNano = RDTSTime - RDTSTickSec;

    //Units are 20 nanoseconds ticks
    double RDTSTimeSec  = 20.e-9 * RDTSTickSec;
    double RDTSTimeNano = 20.    * RDTSTickNano;

    double diffSec = RDTSTimeSec - GenTrigSec - SpillOffset;
    double diffNano = 1.e-09*(RDTSTimeNano - GenTrigNano);

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << "RDTSTimeSec - GenTrigSec " << RDTSTimeSec - GenTrigSec << "\n";
      MF_LOG_INFO("BeamEvent") << "diff: " << diffSec << " " << diffNano << "\n";
    }

    double diff = diffSec + diffNano; 
  
    if( ( fTimingCalibration - fCalibrationTolerance < diff ) && (fTimingCalibration + fCalibrationTolerance > diff) ){

      beamspill->SetActiveTrigger( iT ); 
      beamevt->SetActiveTrigger( iT );
      beamevt->SetT0( beamspill->GetT0( iT ) );

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "FOUND MATCHING TIME!!!" << "\n";
        MF_LOG_INFO("BeamEvent") << "diff: " << diff << "\n";
        MF_LOG_INFO("BeamEvent") << "Set event T0: " << beamevt->GetT0Sec() << " " << beamevt->GetT0Nano() << "\n";
      }

      return;
    }
  }

  MF_LOG_INFO("BeamEvent") << "Could not find matching time " << "\n";
  beamspill->SetUnmatched();
}

void proto::BeamEvent::SetBeamEvent(){

  if( !beamspill->CheckIsMatched() ){
    MF_LOG_INFO("BeamEvent") << "art Event is unmatched to Beam Spill " << "\n";
    return;
  }

  size_t activeTrigger = beamspill->GetActiveTrigger();
  beamevt->SetActiveTrigger( activeTrigger );
  beamevt->SetT0( beamspill->GetT0( activeTrigger ) );
  
  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << "Setting beam event for matched event" << "\n";
    MF_LOG_INFO("BeamEvent") << "SetActiveTrigger " << beamevt->GetActiveTrigger() << "\n";
    MF_LOG_INFO("BeamEvent") << "Set T0  " << beamevt->GetT0Sec() << " " << beamevt->GetT0Nano() << "\n" << "\n"; 
  }


  for( size_t i = 0; i < fDevices.size(); ++i){
    std::string theName = fDevices[i];
    beamevt->SetFBMTrigger( theName, beamspill->GetFBM(theName, activeTrigger) );    
  }

  beamevt->SetTOFs( beamspill->GetMultipleTOFs( activeTrigger ) );
  beamevt->SetTOFChans( beamspill->GetMultipleTOFChans( activeTrigger ) );
  beamevt->SetUpstreamTriggers( beamspill->GetUpstreamTriggers( activeTrigger ) );
  beamevt->SetDownstreamTriggers( beamspill->GetDownstreamTriggers( activeTrigger ) );
  beamevt->SetCalibrations( fTOFCalAA, fTOFCalBA, fTOFCalAB, fTOFCalBB );
  beamevt->CalibrateTOFs();
  beamevt->DecodeTOF();
  beamevt->SetMagnetCurrent( beamspill->GetMagnetCurrent() );
  beamevt->SetCKov0( beamspill->GetCKov0( activeTrigger ) );
  beamevt->SetCKov1( beamspill->GetCKov1( activeTrigger ) );


  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << "beamevt has TOF " << beamevt->GetTOF() 
              << " and TOFChan "    << beamevt->GetTOFChan() << "\n";


    MF_LOG_INFO("BeamEvent") << "beamevt has Magnet Current " << beamevt->GetMagnetCurrent() << "\n";
    MF_LOG_INFO("BeamEvent") << "beamevt CKov0: " << beamevt->GetCKov0Status() << " " << beamevt->GetCKov0Pressure() << "\n";
    MF_LOG_INFO("BeamEvent") << "beamevt CKov1: " << beamevt->GetCKov1Status() << " " << beamevt->GetCKov1Pressure() << "\n";
    MF_LOG_INFO("BeamEvent") << "Finished adding info to beamevt " << "\n";

  }

  C1DB = beamspill->GetCKov0Status( activeTrigger );
  C2DB = beamspill->GetCKov1Status( activeTrigger );
}


void proto::BeamEvent::reset(){
  acqTime = 0;
  acqStampMBPL = 0;
  RDTSTrigger = -1;
  C1DB        = -1;
  C2DB        = -1;
  SpillStart  = -1;
  SpillEnd    = -1;
  SpillOffset = -1;
  
  s11Nano     = -1.;
  s11Sec      = -1.;
  
  ActiveTriggerTime = -1;
  RDTSTime   = 0;
}


////////////////////////
// Producer Method (reads in the event and derives values)
void proto::BeamEvent::produce(art::Event & e){

  reset();

  eventNum = e.event();
  runNum = e.run();
  subRunNum = e.subRun();

  if( e.time().timeHigh() == 0 ) eventTime = e.time().timeLow();
  else eventTime = e.time().timeHigh();


   
  // Create a new beam event (note the "new" here)  
  beamevt = new beam::ProtoDUNEBeamEvent();
  beamspill = new beam::ProtoDUNEBeamSpill();

  //Get the information from the timing system
  //This gets RDTSTrigger and RDTSTime
  GetRawDecoderInfo(e);
  
  //Get Spill Information
  //This stores Spill Start, Spill End, 
  //And Prev Spill Start
  GetSpillInfo(e);

  
  //Check if we have a valid beam trigger
  //If not, just place an empty beamevt
  //and move on
  //
  //Also check if we've gotten good spill info
  //
  //Or if we're forcing to read out the Beamline Info
  if( ( (RDTSTrigger == 12) ) || fForceRead ){

    //Start getting beam event info
    uint64_t fetch_time = uint64_t( RDTSTime * 2e-8 );
    uint64_t fetch_time_down = uint64_t( RDTSTime * 2e-8 );

    if( fPrintDebug )
      MF_LOG_INFO("BeamEvent") << "RDTSTime: " <<  uint64_t( RDTSTime * 2e-8 ) << "\n";

    //Check if we are still using the same spill information.
    //
    //If it's a new spill: Get new info from the database 
    //Each 'parse' command fetches new data from the database
    //
    //If it's the same spill then just pass the old BeamEvent
    //Object. Its 'active trigger' info will be updated below
    //
    //Also: Check if the SpillStart is invalid. If the RDTSTime is  
    //Greater than 5, it's in a new spill, so fetch again
    //
    //Can be overridden with a flag from the fcl
    if( ( ( PrevStart != SpillStart) || ( !SpillStartValid && ( abs(RDTSTimeSec - PrevRDTSTimeSec) > 5 ) ) )
    || fForceNewFetch){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "New spill or forced new fetch. Getting new beamspill info" << "\n";

      //Testing: printing out cache start and end 
      cache_start = bfp->GetCacheStartTime();
      cache_end   = bfp->GetCacheEndTime();

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "cache_start: " << cache_start << "\n";
        MF_LOG_INFO("BeamEvent") << "cache_end: "   << cache_end << "\n";
        MF_LOG_INFO("BeamEvent") << "fetch_time: "  << fetch_time << "\n";
      }
     
      cache_start = bfp_xcet->GetCacheStartTime();
      cache_end   = bfp_xcet->GetCacheEndTime();

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "xcet cache_start: " << cache_start << "\n";
        MF_LOG_INFO("BeamEvent") << "xcet cache_end: "   << cache_end << "\n";
        MF_LOG_INFO("BeamEvent") << "xcet fetch_time: "  << fetch_time << "\n";
      }

      //Not the first event
      if(cache_start > 0 && cache_end > 0){

        //So try filling the cache first with the 'possible' end of spill time
        //then the lower spill time.
        //
        //This is done so that the cache essentially reshuffles where it starts and ends
        //
        //All the checking is done internal to the FillCache method
        //
        //Note: I'm using a loose definition of the start and end of spills
        //      It's really just the maximum and minimum possible vales of those times
        //      since it's not possible to know for certain in any given event. 
        //      (The info does exist in the raw decoder info, but it's not always 
        //       present, so I'm just opting for this)
        try{        
          bfp->FillCache( fetch_time + fFillCacheUp );

          if( fPrintDebug ){
            cache_start = bfp->GetCacheStartTime();
            cache_end   = bfp->GetCacheEndTime();
            MF_LOG_INFO("BeamEvent") << "interim cache_start: " << cache_start << "\n";
            MF_LOG_INFO("BeamEvent") << "interim cache_end: "   << cache_end << "\n";
          }

          bfp->FillCache( fetch_time - fFillCacheDown );
          if( fPrintDebug ){
            cache_start = bfp->GetCacheStartTime();
            cache_end   = bfp->GetCacheEndTime();
            MF_LOG_INFO("BeamEvent") << "new cache_start: " << cache_start << "\n";
            MF_LOG_INFO("BeamEvent") << "new cache_end: "   << cache_end << "\n";
          }

          bfp_xcet->FillCache( fetch_time + fFillCacheUp );
          if( fPrintDebug ){
            cache_start = bfp_xcet->GetCacheStartTime();
            cache_end   = bfp_xcet->GetCacheEndTime();
            MF_LOG_INFO("BeamEvent") << "interim xcet cache_start: " << cache_start << "\n";
            MF_LOG_INFO("BeamEvent") << "interim xcet cache_end: "   << cache_end << "\n";
          }

          bfp_xcet->FillCache( fetch_time - fFillCacheDown );
          if( fPrintDebug ){
            cache_start = bfp_xcet->GetCacheStartTime();
            cache_end   = bfp_xcet->GetCacheEndTime();
            MF_LOG_INFO("BeamEvent") << "new xcet cache_start: " << cache_start << "\n";
            MF_LOG_INFO("BeamEvent") << "new xcet cache_end: "   << cache_end << "\n";
          }
        }
        catch( std::exception e ){
          MF_LOG_WARNING("BeamEvent") << "Could not fill cache\n"; 
        }
      }      
      else{
        //First event, let's get the start of spill info 

        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "First Event: Priming cache\n";

        try{        
          bfp->FillCache( fetch_time - fFillCacheDown );
        }
        catch( std::exception e ){
          MF_LOG_WARNING("BeamEvent") << "Could not fill cache\n"; 
        }
        try{
          bfp_xcet->FillCache( fetch_time - fFillCacheDown );
        }
        catch( std::exception e ){
          MF_LOG_WARNING("BeamEvent") << "Could not fill xcet cache\n"; 
        }

        if( fPrintDebug ){
          cache_start = bfp->GetCacheStartTime();
          cache_end   = bfp->GetCacheEndTime();
          MF_LOG_INFO("BeamEvent") << "new cache_start: " << cache_start << "\n";
          MF_LOG_INFO("BeamEvent") << "new cache_end: "   << cache_end << "\n";

          cache_start = bfp_xcet->GetCacheStartTime();
          cache_end   = bfp_xcet->GetCacheEndTime();
          MF_LOG_INFO("BeamEvent") << "new xcet cache_start: " << cache_start << "\n";
          MF_LOG_INFO("BeamEvent") << "new xcet cache_end: "   << cache_end << "\n";
        }

      }

      // Parse the Time of Flight Counter data for the list
      // of times that we are using
      parseXTOF(fetch_time);
      

      // Parse the Beam postion counter information for the list
      // of time that we are using
      InitXBPFInfo(beamspill);

      if( gotGeneralTrigger ){ 
        parseXBPF(fetch_time);
      }

      parseXCETDB(fetch_time);

      try{
        current = FetchAndReport(fetch_time_down, "dip/acc/NORTH/NP04/POW/MBPL022699:current", bfp);
        gotCurrent = true;
        beamspill->SetMagnetCurrent(current[0]);
      }
      catch( std::exception e){
        MF_LOG_WARNING("BeamEvent") << "Could not get magnet current\n";
        gotCurrent = false;
      }
   
      //Set PrevStart to SpillStart here
      //so that we don't skip in the case 
      //the first event in the spill did not
      //have a good beamline trigger
      PrevStart = SpillStart;
      PrevRDTSTimeSec = RDTSTimeSec;

    }
    else{
      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "Same spill. Reusing beamspill info" << "\n";

      *beamspill = prev_beamspill;
    }

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << "NGoodParticles: " << beamspill->GetNT0()            << "\n";
      MF_LOG_INFO("BeamEvent") << "NTOF0: "          << beamspill->GetNTOF0Triggers()  << "\n";
      MF_LOG_INFO("BeamEvent") << "NTOF1: "          << beamspill->GetNTOF1Triggers()  << "\n";
      MF_LOG_INFO("BeamEvent") << "acqTime: "        << acqTime                      << "\n";
      MF_LOG_INFO("BeamEvent") << "NXBPF: "          << beamspill->GetNFBMTriggers(fDevices[0]) << "\n";
    }

    if( fMatchTime ){
   
      //Now do the matching in Time:
      //Get the conversion from the ProtoDUNE Timing system
      //To the one in the SPS.
      //

      if(SpillStartValid && !fForceMatchS11){
        TimeIn(e, fetch_time_down);

        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "SpillOffset " << SpillOffset << "\n";
        
        //If not successfully timed in, 
        //Oh well. It won't cause a crash
        //But it won't be matched. 
        //
        //Additionally, check if the MBPL timestamp is valid,
        //There are some instances in the database of it being abnormally low
        //If so, then just do the S11 matching
        if( acqStampValid ){
          MatchBeamToTPC();
        }
        else{
          try{
            getS11Info(fetch_time); 
          }
          catch( std::exception e ){
            MF_LOG_WARNING("BeamEvent") << "Could not get S11 Info\n";
          }

          //Again, it won't crash, but it won't match          
          MatchS11ToGen();
        }
      }
      else{
        try{
          getS11Info(fetch_time); 
        }
        catch( std::exception e ){
          MF_LOG_WARNING("BeamEvent") << "Could not get S11 Info\n";
        }

        //Again, it won't crash, but it won't match          
        MatchS11ToGen();
      }


      if( beamspill->CheckIsMatched() ){
        std::pair<double,double> theTime = beamspill->GetT0(beamspill->GetActiveTrigger());
        ActiveTriggerTime = theTime.first + theTime.second*1.e-9;

        if( fPrintDebug ) 
          MF_LOG_INFO("BeamEvent") << "Trigger: " << beamspill->GetActiveTrigger() << " " << ActiveTriggerTime << "\n";       
      
        //Pass the information to the beamevent
        SetBeamEvent();

        MakeTrack( beamspill->GetActiveTrigger() );

        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "Added " << beamevt->GetNBeamTracks() << " tracks to the beam spill" << "\n";

        //Momentum
        if( gotCurrent ){
          MomentumSpec( beamspill->GetActiveTrigger() ); 
        }

        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "Got NRecoBeamMomenta: " << beamevt->GetNRecoBeamMomenta() << "\n";
      }
    }


    //Pass beamspill to the next event;
    //Erase the Track and Reco Momentum info
    prev_beamspill = *beamspill;
    prev_beamspill.ClearBeamTracks();
    prev_beamspill.ClearRecoBeamMomenta();
    prev_beamspill.SetUnmatched();

  }
  //Start of a new spill, but the first event was 
  //Not a beam trigger. In this case, it would not
  //have been filled with info in the block above
  //So let's make it empty so we aren't putting 
  //old spill info in the new event

  //Or. If this was not a beam trigger, and the RDTSTime
  //Comes from a new spill while the SpillStart was invalid, pass it. 
  else if( ( ( PrevStart != SpillStart ) || ( !SpillStartValid && ( abs(RDTSTimeSec - PrevRDTSTimeSec) > 5 ) ) ) 
    && RDTSTrigger != 12 ){
    prev_beamspill = *beamspill;   
  }
  
  //Can Remove BITrigger,CTBTimestamp
  beamevt->SetBITrigger(-1);
  beamevt->SetTimingTrigger(RDTSTrigger);
  beamevt->SetSpillStart(SpillStart);
  beamevt->SetSpillOffset(SpillOffset);
  beamevt->SetCTBTimestamp( -1. );
  beamevt->SetRDTimestamp( RDTSTime );

  std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent> > beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
  beamData->push_back(beam::ProtoDUNEBeamEvent(*beamevt));
  e.put(std::move(beamData));
  delete beamevt;
  delete beamspill;
 
  // Write out the to tree
  if( fSaveOutTree )fOutTree->Fill();
 
}
// END BeamEvent::produce
////////////////////////

void proto::BeamEvent::InitXBPFInfo(beam::ProtoDUNEBeamSpill * beamspill){
  // Places a dummy trigger vector for each device

  // Make a vector the names of each of the devices being readout
  std::vector<std::string> monitors;
  size_t nDev = 0; 

  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
    // Put the current device name on the list
    monitors.push_back(name);
    nDev++;
  }

  beamspill->InitFBMs(monitors);
}
// END BeamEvent::InitXBPFInfo
////////////////////////

void proto::BeamEvent::getS11Info(uint64_t time){
  MF_LOG_INFO("BeamEvent") << "Getting S11 Info " << "\n";

  std::vector<double> coarseS11         = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:coarse[]", bfp);
  std::vector<double> fracS11           = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:frac[]", bfp); 
  std::vector<double> secondsS11        = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:seconds[]", bfp); 
  std::vector<double> timestampCountS11 = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:timestampCount", bfp); 
  
  int s11Count = (int)timestampCountS11[0];

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Found " << s11Count << " returned signals" << "\n";

  //Separates seconds portion of the ticks 
  //From the nanoseconds
  long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
  RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
  long long RDTSTickNano = RDTSTime - RDTSTickSec;

  //Units are 20 nanoseconds ticks
  double RDTSTimeSec  = 20.e-9 * RDTSTickSec;
  double RDTSTimeNano = 20.    * RDTSTickNano;


  for(int i = 0; i < s11Count; ++i){
    double nano = 8.*coarseS11[i] + fracS11[i]/512.;

    double diffSec = secondsS11[2*i + 1] - RDTSTimeSec - fOffsetTAI;
    double diffNano = nano - RDTSTimeNano;

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << i << " diffSec " << diffSec << "\n"; 
      MF_LOG_INFO("BeamEvent") << i << " diffNano " << diffNano << "\n"; 
    }

    double diff = diffSec + 1.e-9*diffNano;
    if( fRDTSToS11Lower < diff && diff < fRDTSToS11Upper ){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "FOUND Match between S11 and RDTS" << "\n";  

      s11Nano = nano;
      s11Sec  = secondsS11[2*i + 1] - fOffsetTAI;
      return;
    }
  }

  MF_LOG_WARNING("BeamEvent") << "Could not match RDTS to S11\n";

}


void proto::BeamEvent::parseXTOF(uint64_t time){

  std::vector<double> coarseGeneralTrigger;         
  std::vector<double> fracGeneralTrigger;           
  std::vector<double> acqStampGeneralTrigger;       
  std::vector<double> secondsGeneralTrigger;        
  std::vector<double> timestampCountGeneralTrigger; 
  

  try{
    coarseGeneralTrigger         = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:coarse[]", bfp);
    fracGeneralTrigger           = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:frac[]", bfp); 
    acqStampGeneralTrigger       = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:acqStamp[]", bfp); 
    secondsGeneralTrigger        = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:seconds[]", bfp); 
    timestampCountGeneralTrigger = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:timestampCount", bfp); 

    if( fPrintDebug )
      MF_LOG_INFO("BeamEvent") << "timestampCounts: " << timestampCountGeneralTrigger[0] << "\n";

    gotGeneralTrigger = true;

    uint64_t low = (uint64_t)acqStampGeneralTrigger[1];
    uint64_t high = (uint64_t)acqStampGeneralTrigger[0];
    high = high << 32;
    acqTime = ( high | low ) / 1000000000.; 

  }
  catch(std::exception e){
    MF_LOG_WARNING("BeamEvent") << "Could not get GeneralTrigger information!!" << "\n";
    gotGeneralTrigger = false;
    return;
  }
  
  std::vector<double> coarseTOF1A;  
  std::vector<double> fracTOF1A;    
  std::vector<double> secondsTOF1A; 

  std::vector<double> coarseTOF1B;  
  std::vector<double> fracTOF1B;    
  std::vector<double> secondsTOF1B; 

  std::vector<double> coarseTOF2A;  
  std::vector<double> fracTOF2A;    
  std::vector<double> secondsTOF2A; 

  std::vector<double> coarseTOF2B;  
  std::vector<double> fracTOF2B;    
  std::vector<double> secondsTOF2B; 

  try{
    coarseTOF1A  = FetchAndReport(time, fXTOFPrefix + fTOF1A + ":coarse[]", bfp);
    fracTOF1A    = FetchAndReport(time, fXTOFPrefix + fTOF1A + ":frac[]", bfp);
    secondsTOF1A = FetchAndReport(time, fXTOFPrefix + fTOF1A + ":seconds[]", bfp);

    coarseTOF1B  = FetchAndReport(time, fXTOFPrefix + fTOF1B + ":coarse[]", bfp);
    fracTOF1B    = FetchAndReport(time, fXTOFPrefix + fTOF1B + ":frac[]", bfp);
    secondsTOF1B = FetchAndReport(time, fXTOFPrefix + fTOF1B + ":seconds[]", bfp);

    coarseTOF2A  = FetchAndReport(time, fXTOFPrefix + fTOF2A + ":coarse[]", bfp);
    fracTOF2A    = FetchAndReport(time, fXTOFPrefix + fTOF2A + ":frac[]", bfp);
    secondsTOF2A = FetchAndReport(time, fXTOFPrefix + fTOF2A + ":seconds[]", bfp);

    coarseTOF2B  = FetchAndReport(time, fXTOFPrefix + fTOF2B + ":coarse[]", bfp);
    fracTOF2B    = FetchAndReport(time, fXTOFPrefix + fTOF2B + ":frac[]", bfp);
    secondsTOF2B = FetchAndReport(time, fXTOFPrefix + fTOF2B + ":seconds[]", bfp);

    gotTOFs = true;
  }
  catch(std::exception e){
    MF_LOG_WARNING("BeamEvent") << "Could not get TOF information!!" << "\n";
    gotTOFs = false;
  }

  std::vector<std::pair<double,double>> unorderedGenTrigTime;
  std::vector<std::pair<double,double>> unorderedTOF1ATime;
  std::vector<std::pair<double,double>> unorderedTOF1BTime;
  std::vector<std::pair<double,double>> unorderedTOF2ATime;
  std::vector<std::pair<double,double>> unorderedTOF2BTime;


  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << i << " " << secondsGeneralTrigger[2*i + 1] << " "  << 8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512. << "\n";
      MF_LOG_INFO("BeamEvent") << "\t" << std::setw(15) << secondsGeneralTrigger[2*i + 1] + 1.e-9*(8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512.) << "\n";
    }

    //2*i + 1 because the format is weird
    fGenTrigSec    = secondsGeneralTrigger[2*i + 1];
    fGenTrigCoarse = coarseGeneralTrigger[i];
    fGenTrigFrac   = fracGeneralTrigger[i];

    unorderedGenTrigTime.push_back( std::make_pair(fGenTrigSec, (fGenTrigCoarse*8. + fGenTrigFrac/512.)) );

    if (fGenTrigCoarse == 0.0 && fGenTrigFrac == 0.0 && fGenTrigSec == 0.0) break;
    if(fDebugTOFs)fGenTrigTree->Fill();
  }

  if( gotTOFs ){

    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "TOF1A " << i << " " << secondsTOF1A[2*i+1] << " "  << 8.*coarseTOF1A[i] +  fracTOF1A[i]/512. << "\n";

      fXTOF1ACoarse = coarseTOF1A[i];
      fXTOF1AFrac   = fracTOF1A[i];
      fXTOF1ASec    = secondsTOF1A[2*i + 1];

      if(fXTOF1ASec < secondsTOF1A[1]) break; 

      if (fXTOF1ACoarse == 0.0 && fXTOF1AFrac == 0.0 && fXTOF1ASec == 0.0) break;
      unorderedTOF1ATime.push_back(std::make_pair(fXTOF1ASec, (fXTOF1ACoarse*8. + fXTOF1AFrac/512.)) );

      if(fDebugTOFs)fXTOF1ATree->Fill();
    }
      
    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "TOF1B " << i << " " << secondsTOF1B[2*i+1] << " "  << 8.*coarseTOF1B[i] + fracTOF1B[i]/512. << "\n";

      fXTOF1BCoarse = coarseTOF1B[i];
      fXTOF1BFrac   = fracTOF1B[i];
      fXTOF1BSec    = secondsTOF1B[2*i + 1];

      if(fXTOF1BSec < secondsTOF1B[1]) break; 

      if (fXTOF1BCoarse == 0.0 && fXTOF1BFrac == 0.0 && fXTOF1BSec == 0.0) break;
      unorderedTOF1BTime.push_back(std::make_pair(fXTOF1BSec, (fXTOF1BCoarse*8. + fXTOF1BFrac/512.)) );
      if(fDebugTOFs)fXTOF1BTree->Fill();
    }

    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "TOF2A " << i << " " << secondsTOF2A[2*i+1] << " "  << 8.*coarseTOF2A[i] +  fracTOF2A[i]/512. << "\n";

      fXTOF2ACoarse = coarseTOF2A[i];
      fXTOF2AFrac   = fracTOF2A[i];
      fXTOF2ASec    = secondsTOF2A[2*i + 1];

      if(fXTOF2ASec < secondsTOF2A[1]) break; 

      if (fXTOF2ACoarse == 0.0 && fXTOF2AFrac == 0.0 && fXTOF2ASec == 0.0) break;
      unorderedTOF2ATime.push_back(std::make_pair(fXTOF2ASec, (fXTOF2ACoarse*8. + fXTOF2AFrac/512.)) );

      if(diff2A.size()) diff2A.clear();
      for(size_t j = 0; j < timestampCountGeneralTrigger[0]; ++j){
        diff2A.push_back( 1.e9*(unorderedTOF2ATime.back().first - unorderedGenTrigTime[j].first) + (unorderedTOF2ATime.back().second - unorderedGenTrigTime[j].second) );
      }
      
      if(fDebugTOFs)fXTOF2ATree->Fill();

    }  

    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "TOF2B " << i << " " << secondsTOF2B[2*i+1] << " "  << 8.*coarseTOF2B[i] +  fracTOF2B[i]/512. << "\n";

      fXTOF2BCoarse = coarseTOF2B[i];
      fXTOF2BFrac   = fracTOF2B[i];
      fXTOF2BSec    = secondsTOF2B[2*i + 1];

      if(fXTOF2BSec < secondsTOF2B[1]) break; 
      if (fXTOF2BCoarse == 0.0 && fXTOF2BFrac == 0.0 && fXTOF2BSec == 0.0) break;

      unorderedTOF2BTime.push_back(std::make_pair(fXTOF2BSec, (fXTOF2BCoarse*8. + fXTOF2BFrac/512.) ));

      if(diff2B.size()) diff2B.clear();

      for(size_t j = 0; j < timestampCountGeneralTrigger[0]; ++j){
        diff2B.push_back( 1.e9*(unorderedTOF2BTime.back().first - unorderedGenTrigTime[j].first) + (unorderedTOF2BTime.back().second - unorderedGenTrigTime[j].second) );
      }

      if(fDebugTOFs)fXTOF2BTree->Fill();
    }
  }

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "NGenTrigs: " << timestampCountGeneralTrigger[0] << " NTOF2s: " << unorderedTOF2ATime.size() + unorderedTOF2BTime.size() << "\n";

  for(size_t iT = 0; iT < unorderedGenTrigTime.size(); ++iT){
    
    bool found_TOF = false;

    double the_gen_sec = unorderedGenTrigTime[iT].first;
    double the_gen_ns = unorderedGenTrigTime[iT].second;

    //1A2A = 0; 1B2A = 1, 1A2B = 2,  1B2B = 3
    //Add 1 for 1B, add 2 for 2B
    int channel = 0;

    std::vector< double > possibleTOF; 
    std::vector< int >    possibleTOFChan; 
    std::vector< size_t > UpstreamTriggers;
    std::vector< size_t > DownstreamTriggers;

    if( gotTOFs ){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "Gen: " << the_gen_sec << " " << the_gen_ns << "\n";

      //First check 2A
      for(size_t ip2A = 0; ip2A < unorderedTOF2ATime.size(); ++ip2A){
        double TOF2A_sec = unorderedTOF2ATime[ip2A].first;
        double TOF2A_ns  = unorderedTOF2ATime[ip2A].second;         
        double delta_2A = 1.e9*(the_gen_sec - TOF2A_sec) + the_gen_ns - TOF2A_ns;

        if( delta_2A < 0. ) break;
        else if( delta_2A > fDownstreamToGenTrig ) continue;

        //if here, 0. < delta < 50ns
        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "Found match 2A to Gen" << "\n";

        //So check 1A and 1B
        for(size_t ip1A = 0; ip1A < unorderedTOF1ATime.size(); ++ip1A){
          double TOF1A_sec = unorderedTOF1ATime[ip1A].first;
          double TOF1A_ns  = unorderedTOF1ATime[ip1A].second;         
          double delta = 1.e9*( TOF2A_sec - TOF1A_sec ) + TOF2A_ns - TOF1A_ns;

          if( delta < 0. ) break;
          else if( delta > 0. && delta < fUpstreamToDownstream){
            if( fPrintDebug )
              MF_LOG_INFO("BeamEvent") << "Found match 1A to 2A " << delta << "\n";

            found_TOF = true;
            channel = k1A2A;

            possibleTOF.push_back( delta ); 
            possibleTOFChan.push_back( channel );

            UpstreamTriggers.push_back( ip1A );
            DownstreamTriggers.push_back( ip2A );
          }
          else continue; 
        }
        for(size_t ip1B = 0; ip1B < unorderedTOF1BTime.size(); ++ip1B){
          double TOF1B_sec = unorderedTOF1BTime[ip1B].first;
          double TOF1B_ns  = unorderedTOF1BTime[ip1B].second;         
          double delta = 1.e9*( TOF2A_sec - TOF1B_sec ) + TOF2A_ns - TOF1B_ns;

          if( delta < 0. ) break;
          else if( delta > 0. && delta < fUpstreamToDownstream){

            if( fPrintDebug )
              MF_LOG_INFO("BeamEvent") << "Found match 1B to 2A " << delta << "\n";

            found_TOF = true;
            channel = k1B2A;

            possibleTOF.push_back( delta ); 
            possibleTOFChan.push_back( channel );

            UpstreamTriggers.push_back( ip1B );
            DownstreamTriggers.push_back( ip2A );
          }
          else continue; 
        }
      }

      //Then check 2B 
      for(size_t ip2B = 0; ip2B < unorderedTOF2BTime.size(); ++ip2B){
        double TOF2B_sec = unorderedTOF2BTime[ip2B].first;
        double TOF2B_ns  = unorderedTOF2BTime[ip2B].second;         
        double delta_2B = 1.e9*(the_gen_sec - TOF2B_sec) + the_gen_ns - TOF2B_ns;


        if( delta_2B < 0. ) break;
        else if( delta_2B > fDownstreamToGenTrig ) continue;

        //if here, 0. < delta < 50ns
        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "Found match 2B to Gen" << "\n";

        //So check 1A and 1B
        for(size_t ip1A = 0; ip1A < unorderedTOF1ATime.size(); ++ip1A){
          double TOF1A_sec = unorderedTOF1ATime[ip1A].first;
          double TOF1A_ns  = unorderedTOF1ATime[ip1A].second;         
          double delta = 1.e9*( TOF2B_sec - TOF1A_sec ) + TOF2B_ns - TOF1A_ns;


          if( delta < 0. ) break;
          else if( delta > 0. && delta < fUpstreamToDownstream){
                              
            if( fPrintDebug )
              MF_LOG_INFO("BeamEvent") << "Found match 1A to 2B " << delta << "\n";
             
            found_TOF = true;
            channel = k1A2B;

            possibleTOF.push_back( delta ); 
            possibleTOFChan.push_back( channel );

            UpstreamTriggers.push_back( ip1A );
            DownstreamTriggers.push_back( ip2B );
          }
          else continue; 
        }
        for(size_t ip1B = 0; ip1B < unorderedTOF1BTime.size(); ++ip1B){
          double TOF1B_sec = unorderedTOF1BTime[ip1B].first;
          double TOF1B_ns  = unorderedTOF1BTime[ip1B].second;         
          double delta = 1.e9*( TOF2B_sec - TOF1B_sec ) + TOF2B_ns - TOF1B_ns;

          
          if( delta < 0. ) break;
          else if( delta > 0. && delta < fUpstreamToDownstream){

            if( fPrintDebug )
              MF_LOG_INFO("BeamEvent") << "Found match 1B to 2B " << delta << "\n";

            found_TOF = true;
            channel = k1B2B;

            possibleTOF.push_back( delta ); 
            possibleTOFChan.push_back( channel );

            UpstreamTriggers.push_back( ip1B );
            DownstreamTriggers.push_back( ip2B );
          }
          else continue;
        }
      }

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "Found " << possibleTOF.size() << " matched TOFs" << "\n";
    }

    if( !found_TOF ){

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "No matching TOFs found. Placing dummy\n";

      possibleTOF.push_back( 0. );
      possibleTOFChan.push_back( -1 );
      UpstreamTriggers.push_back(0);
      DownstreamTriggers.push_back(0);
    }

    beamspill->AddT0(std::make_pair(the_gen_sec - fOffsetTAI, the_gen_ns));
    beamspill->AddMultipleTOFs( possibleTOF );
    beamspill->AddMultipleTOFChans( possibleTOFChan );
    beamspill->AddUpstreamTriggers( UpstreamTriggers );
    beamspill->AddDownstreamTriggers( DownstreamTriggers );

  }

}
// END BeamEvent::parseXTOF
////////////////////////

void proto::BeamEvent::parseXCETDB(uint64_t time){

  if(fCKov1 != ""){  

    try{
      std::vector< double > pressureCKov1 = FetchAndReport(time, fXCETPrefix + fCKov1 + ":pressure", bfp);
      std::vector< double > countsCKov1   = FetchAndReport(time, fXCETPrefix + fCKov1 + ":counts", bfp);

      CKov1Pressure = pressureCKov1[0];
      CKov1Counts   = countsCKov1[0];
    }
    catch( std::exception e){
      MF_LOG_WARNING("BeamEvent") << "Could not get Cerenkov 1 Pressure\n";
      CKov1Pressure = 0.; 
      CKov1Counts   = 0.;
    }

  }
  
  if(fCKov2 != ""){
    try{
      std::vector< double > pressureCKov2 = FetchAndReport(time, fXCETPrefix + fCKov2 + ":pressure", bfp);
      std::vector< double > countsCKov2   = FetchAndReport(time, fXCETPrefix + fCKov2 + ":counts", bfp);
      CKov2Pressure = pressureCKov2[0];
      CKov2Counts   = countsCKov2[0];
    }
    catch( std::exception e){
      MF_LOG_WARNING("BeamEvent") << "Could not get Cerenkov 2 Pressure\n";
      CKov2Pressure = 0.; 
      CKov2Counts   = 0.;
    }  
  }

  std::vector< double > XCET1_timestamps, XCET2_timestamps;
  std::vector< double > XCET1_seconds,    XCET2_seconds;   
  std::vector< double > XCET1_frac,       XCET2_frac;      
  std::vector< double > XCET1_coarse,     XCET2_coarse;    

  bool fetched_XCET1=false;
  bool fetched_XCET2=false; 
  if( fXCET1 != "" ){
    try{ 
      XCET1_seconds    = FetchAndReport( time + fXCETFetchShift, fXCET1 + ":SECONDS" , bfp_xcet);
      XCET1_frac       = FetchAndReport( time + fXCETFetchShift, fXCET1 + ":FRAC" , bfp_xcet);
      XCET1_coarse     = FetchAndReport( time + fXCETFetchShift, fXCET1 + ":COARSE" , bfp_xcet);
      fetched_XCET1 = true;

      if( fXCETDebug ){
        for( size_t i = 0; i < XCET1_seconds.size(); ++i ){
          std::cout << i << " " << XCET1_seconds[i] - fOffsetTAI << " " << (8.*XCET1_coarse[i] + XCET1_frac[i] / 512.) << std::endl;
        }
      }

    }
    catch( const std::exception &e ){
      MF_LOG_WARNING("BeamEvent") << "Could not get XCET1 info\n";
      fetched_XCET1 = false;
    }
  }
  if( fXCET2 != "" ){
    try{ 
      XCET2_seconds    = FetchAndReport( time + fXCETFetchShift, fXCET2 + ":SECONDS" , bfp_xcet);
      XCET2_frac       = FetchAndReport( time + fXCETFetchShift, fXCET2 + ":FRAC" , bfp_xcet);
      XCET2_coarse     = FetchAndReport( time + fXCETFetchShift, fXCET2 + ":COARSE" , bfp_xcet);
      fetched_XCET2 = true;

      if( fXCETDebug ){
        for( size_t i = 0; i < XCET2_seconds.size(); ++i ){
          std::cout << i << " " << XCET2_seconds[i] - fOffsetTAI << " " << (8.*XCET2_coarse[i] + XCET2_frac[i] / 512.) << std::endl;
        }
      }
    }
    catch( const std::exception &e ){
      MF_LOG_WARNING("BeamEvent") << "Could not get XCET2 info\n";

      fetched_XCET2 = false;
    }
  }



  //Go through the general triggers and try to match. If one can't be found, then just add a 0
  for( size_t i = 0; i < beamspill->GetNT0(); ++i ){

    if( fXCETDebug ) std::cout << "GenTrig: " << i << " " << beamspill->GetT0Sec(i) << " " << beamspill->GetT0Nano(i) << std::endl;
    
    beam::CKov status_1;
    status_1.pressure = CKov1Pressure;
    status_1.timeStamp = -1.;

    if( !fetched_XCET1 ) status_1.trigger = -1;
    else{

      status_1.trigger = 0;

      for( size_t ic1 = 0; ic1 < XCET1_seconds.size(); ++ic1 ){
        double delta = 1.e9 * ( beamspill->GetT0Sec(i) - (XCET1_seconds[ic1] - fOffsetTAI) );
        delta += ( beamspill->GetT0Nano(i) - (8.*XCET1_coarse[ic1] + XCET1_frac[ic1] / 512.) );

        if( fXCETDebug ) std::cout << "XCET1 delta: " << delta << std::endl;

        if( fabs(delta) < 500. ){
          if( fXCETDebug ) std::cout << "Found matching XCET1 trigger " << XCET1_seconds[ic1] - fOffsetTAI << " " << (8.*XCET1_coarse[ic1] + XCET1_frac[ic1] / 512.) << " " << delta << std::endl;
          status_1.trigger = 1;
          status_1.timeStamp = delta;
          break;
        }     
      }
    }
    beamspill->AddCKov0( status_1 );

    beam::CKov status_2;
    status_2.pressure = CKov2Pressure;
    status_2.timeStamp = -1.;

    if( !fetched_XCET2 ) status_2.trigger = -1;
    else{

      status_2.trigger = 0;

      for( size_t ic2 = 0; ic2 < XCET2_seconds.size(); ++ic2 ){

        double delta = 1.e9 * ( beamspill->GetT0Sec(i) - (XCET2_seconds[ic2] - fOffsetTAI) );
        delta += ( beamspill->GetT0Nano(i) - (8.*XCET2_coarse[ic2] + XCET2_frac[ic2] / 512.) );

        if( fXCETDebug ) std::cout << "XCET2 delta: " << delta << std::endl;

        if( fabs(delta) < 500. ){
          if( fXCETDebug ) std::cout << "Found matching XCET2 trigger " << XCET2_seconds[ic2] - fOffsetTAI << " " << (8.*XCET2_coarse[ic2] + XCET2_frac[ic2] / 512.) << " " << delta << std::endl;
          status_2.trigger = 1;
          status_2.timeStamp = delta;
          break;
        }     
      }
    }
    beamspill->AddCKov1( status_2 );

  }

  if( fXCETDebug ){
    MF_LOG_INFO("BeamEvent") << "GeneralTriggers: " << beamspill->GetNT0() << std::endl;
    MF_LOG_INFO("BeamEvent") << "XCET1: " << beamspill->GetNCKov0() << std::endl;
    MF_LOG_INFO("BeamEvent") << "XCET2: " << beamspill->GetNCKov1() << std::endl;

    int nxcet1 = 0, nxcet2 = 0;

    for( size_t i = 0; i < beamspill->GetNT0(); ++i ){
      if( beamspill->GetCKov0Status(i) == 1 ) nxcet1++ ;
      if( beamspill->GetCKov1Status(i) == 1 ) nxcet2++ ;
    }

    if( nxcet1 != CKov1Counts ) MF_LOG_WARNING("BeamEvent") << "CKov1 counts differ. In spill: " << nxcet1 << " From counts: " << CKov1Counts << "\n";
    if( nxcet2 != CKov2Counts ) MF_LOG_WARNING("BeamEvent") << "CKov2 counts differ. In spill: " << nxcet2 << " From counts: " << CKov2Counts << "\n";
    if( fetched_XCET1 ) MF_LOG_INFO("BeamEvent") << "XCET1: " << XCET1_seconds.size() << " " << nxcet1 << std::endl;
    if( fetched_XCET2 ) MF_LOG_INFO("BeamEvent") << "XCET2: " << XCET2_seconds.size() << " " << nxcet2 << std::endl;
  }

}


////////////////////////
// 
void proto::BeamEvent::parseGeneralXBPF(std::string name, uint64_t time, size_t ID){

  // Retrieve the number of counts in the BPF
  std::vector<double> counts;

  try{
    counts = FetchAndReport(time, fXBPFPrefix + name + ":countsRecords[]", bfp);
  }
  catch( std::exception e){
    MF_LOG_WARNING("BeamEvent") << "Could not fetch " << fXBPFPrefix + name + ":countsRecords[]\n";
    return;
  }
  
  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Counts: " << counts[1] << "\n";

  std::vector<double> data;
  try{
    data = FetchAndReport(time, fXBPFPrefix + name + ":eventsData[]", bfp);
  }
  catch( std::exception e){
    MF_LOG_WARNING("BeamEvent") << "Could not fetch " << fXBPFPrefix + name + ":eventsData[]\n";
    return;
  }

  if( fPrintDebug ){
    // If the number of counts is larger than the number of general triggers
    // make note
    if(counts[1] != beamspill->GetNT0()){
      MF_LOG_WARNING("BeamEvent") << "WARNING MISMATCH " << counts[1] << " " << beamspill->GetNT0() << "\n";
    }
  }
  
  
  beam::FBM fbm;
  fbm.ID = ID;
  fbm.fibers = {};
  std::uninitialized_fill( std::begin(fbm.fiberData), std::end(fbm.fiberData), 0. );
  std::uninitialized_fill( std::begin(fbm.timeData), std::end(fbm.timeData), 0. );
  fbm.timeStamp = 0.;
  fbm.decoded = false;
  fbm.active = std::vector<short>();
    
  //Use this just in case any are out of sync?
  //Shouldn't be, but just to be safe...
  //Helps cut down on time
  std::vector<size_t> leftOvers;      
  for(size_t lo = 0; lo < beamspill->GetNT0(); ++lo){
    leftOvers.push_back(lo);
  }
 
  for(size_t i = 0; i < counts[1]; ++i){      
    for(int j = 0; j < 10; ++j){

      double theData = data[20*i + (2*j + 1)];

      if(j < 4)
	fbm.timeData[j] = theData;           
      else
	fbm.fiberData[j - 4] = theData;
       
    } 

    // Check the time data for corruption
    if(fbm.timeData[1] < .0000001)
      continue;
     

    for(std::vector<size_t>::iterator ip = leftOvers.begin(); ip != leftOvers.end(); ++ip){

      // Compute the time delta between the timeStamp and the T0, see if it's less than 500ns away

      double delta = (beamspill->GetT0Sec(*ip)  - (fbm.timeData[3] - fOffsetTAI) );
      delta += 1.e-9*( beamspill->GetT0Nano(*ip) - 8.*fbm.timeData[2] ); 

      if( delta < -1.e-7 && delta > -10.e-7 ){

	beamspill->ReplaceFBMTrigger(name, fbm, *ip);
	leftOvers.erase(ip);
	break;
      } 
    } 
  } 

  if( fPrintDebug ){
    if( leftOvers.size() ){
      MF_LOG_WARNING("BeamEvent") << "Warning! Could not match to Good Particles: " << "\n";
      for( size_t ip = 0; ip < leftOvers.size(); ++ip){
        MF_LOG_WARNING("BeamEvent") << leftOvers[ip] << " ";
      }
      MF_LOG_WARNING("BeamEvent") << "\n";
    }
  }

  for(size_t i = 0; i < beamspill->GetNFBMTriggers(name); ++i){
    beamspill->DecodeFibers(name,i);
  } 

}

////////////////////////
// 
void proto::BeamEvent::parseXBPF(uint64_t time){
  for(size_t d = 0; d < fDevices.size(); ++d){
    std::string name = fDevices[d];
    parseGeneralXBPF(name, time, d);
  }  
}
// END BeamEvent::parseXBFP
////////////////////////


void proto::BeamEvent::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  
  if( fDebugMomentum ){
    fFullMomentum = tfs->make<TH1F>("FullMomentum", "", 150, 0, 15.);
    fCutMomentum = tfs->make<TH1F>("CutMomentum", "", 150, 0, 15.);
  }

  if( fSaveOutTree ){
    fOutTree = tfs->make<TTree>("tree", "lines"); 
    fOutTree->Branch("Time", &eventTime);
    fOutTree->Branch("RDTS", &RDTSTime);
    fOutTree->Branch("RDTSTrigger", &RDTSTrigger);
    fOutTree->Branch("Event", &eventNum);
    fOutTree->Branch("SpillStart", &SpillStart);
    fOutTree->Branch("SpillEnd", &SpillEnd);
    fOutTree->Branch("SpillOffset", &SpillOffset);
    fOutTree->Branch("ActiveTriggerTime", &ActiveTriggerTime);
    fOutTree->Branch("Run",   &runNum);
    fOutTree->Branch("Subrun", &subRunNum);
    fOutTree->Branch("Pressure1", &CKov1Pressure);
    fOutTree->Branch("Counts1", &CKov1Counts);
    fOutTree->Branch("Pressure2", &CKov2Pressure);
    fOutTree->Branch("Counts2", &CKov2Counts);
    fOutTree->Branch("acqTime",        &acqTime);
    fOutTree->Branch("acqStampMBPL",   &acqStampMBPL);
    fOutTree->Branch("C1DB",        &C1DB);
    fOutTree->Branch("C2DB",        &C2DB);
    fOutTree->Branch("s11Sec",      &s11Sec);
    fOutTree->Branch("s11Nano",      &s11Nano);
  }    

  if( fDebugTOFs ){
    fGenTrigTree = tfs->make<TTree>("GenTrig","");
    fGenTrigTree->Branch("coarse", &fGenTrigCoarse);
    fGenTrigTree->Branch("frac", &fGenTrigFrac);
    fGenTrigTree->Branch("sec", &fGenTrigSec);
  
    fXTOF1ATree = tfs->make<TTree>("TOF1A","");
    fXTOF1ATree->Branch("coarse", &fXTOF1ACoarse);
    fXTOF1ATree->Branch("frac", &fXTOF1AFrac);
    fXTOF1ATree->Branch("sec", &fXTOF1ASec);
  
    fXTOF1BTree = tfs->make<TTree>("TOF1B","");
    fXTOF1BTree->Branch("coarse", &fXTOF1BCoarse);
    fXTOF1BTree->Branch("frac", &fXTOF1BFrac);
    fXTOF1BTree->Branch("sec", &fXTOF1BSec);
  
    fXTOF2ATree = tfs->make<TTree>("TOF2A","");
    fXTOF2ATree->Branch("coarse", &fXTOF2ACoarse);
    fXTOF2ATree->Branch("frac", &fXTOF2AFrac);
    fXTOF2ATree->Branch("sec", &fXTOF2ASec);
    fXTOF2ATree->Branch("diff2A", &diff2A);
  
    fXTOF2BTree = tfs->make<TTree>("TOF2B","");
    fXTOF2BTree->Branch("coarse", &fXTOF2BCoarse);
    fXTOF2BTree->Branch("frac", &fXTOF2BFrac);
    fXTOF2BTree->Branch("sec", &fXTOF2BSec);
    fXTOF2BTree->Branch("diff2B", &diff2B);
  }

  //Tells IFBeam to print out debug statements
  ifbeam_ns::BeamFolder::_debug = fIFBeamDebug;

  bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
  if( fPrintDebug ){ 
    MF_LOG_INFO("BeamEvent") << "%%%%%%%%%% Got beam folder %%%%%%%%%%\n"; 
    MF_LOG_INFO("BeamEvent") << "%%%%%%%%%% Setting TimeWindow: " << fTimeWindow << " %%%%%%%%%%\n";
  }

  bfp->set_epsilon( fBFEpsilon );

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "%%%%%%%%%% Set beam epislon " << fBFEpsilon << " %%%%%%%%%%\n";

  bfp_xcet = ifb->getBeamFolder(fXCETBundleName,fURLStr,fTimeWindow);

  if( fPrintDebug ){ 
    MF_LOG_INFO("BeamEvent") << "%%%%%%%%%% Got beam folder %%%%%%%%%%\n"; 
    MF_LOG_INFO("BeamEvent") << "%%%%%%%%%% Setting TimeWindow: " << fTimeWindow << " %%%%%%%%%%\n";
  }
 
  bfp_xcet->set_epsilon( fXCETEpsilon );

  if( fPrintDebug ) 
    MF_LOG_INFO("BeamEvent") << "%%%%%%%%%% Set beam epislon " << fBFEpsilon << " %%%%%%%%%%\n";


  //Rotate the basis vectors of the FBMs
  BeamMonitorBasisVectors();
}

uint64_t proto::BeamEvent::joinHighLow(double high, double low){

  uint64_t low64 = (uint64_t)low;

  uint64_t high64 = (uint64_t)high;

  high64 = high64 << 32;
  uint64_t joined = high64 | low64;

  return joined; 
}

TVector3 proto::BeamEvent::ConvertProfCoordinates(double x, double y, double z, double zOffset){
  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*fBMBasisX.X() + y*fBMBasisY.X() + /*(z-zOffset)*fBMBasisZ.X()*/ + off*fabs(fBMBasisZ.X());
  double newY = x*fBMBasisX.Y() + y*fBMBasisY.Y() + /*(z-zOffset)*fBMBasisZ.Y()*/ + off*fabs(fBMBasisZ.Y());
  double newZ = x*fBMBasisX.Z() + y*fBMBasisY.Z() + /*(z-zOffset)              */ - off*fabs(fBMBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void proto::BeamEvent::BeamMonitorBasisVectors(){
  RotateMonitorVector(fBMBasisX);
  RotateMonitorVector(fBMBasisY);
  RotateMonitorVector(fBMBasisZ);

  rotated = true;
}

void proto::BeamEvent::RotateMonitorVector(TVector3 &vec){
  vec.RotateX( fRotateMonitorYZ * TMath::Pi()/180. );
  vec.RotateY( fRotateMonitorXZ * TMath::Pi()/180. );
}

void proto::BeamEvent::MakeTrack(size_t theTrigger){
  
  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Making Track for time: " << beamspill->GetT0Sec(theTrigger) << " " << beamspill->GetT0Nano(theTrigger) << "\n";

  //Get the active fibers from the upstream tracking XBPF
  std::vector<short> firstUpstreamFibers  = beamspill->GetActiveFibers(firstUpstreamName, theTrigger);
  std::vector<short> secondUpstreamFibers = beamspill->GetActiveFibers(secondUpstreamName, theTrigger);

  if( fPrintDebug ){ 
    MF_LOG_INFO("BeamEvent") << firstUpstreamName << " has " << firstUpstreamFibers.size() << " active fibers"<< "\n";
    for(size_t i = 0; i < firstUpstreamFibers.size(); ++i){
      MF_LOG_INFO("BeamEvent") << firstUpstreamFibers[i] << " ";
    }
    MF_LOG_INFO("BeamEvent") << "\n";

    MF_LOG_INFO("BeamEvent") << secondUpstreamName << " has " << secondUpstreamFibers.size() << " active fibers" << "\n";
    for(size_t i = 0; i < secondUpstreamFibers.size(); ++i){
      MF_LOG_INFO("BeamEvent") << secondUpstreamFibers[i] << " ";
    }
    MF_LOG_INFO("BeamEvent") << "\n";
  } 
  //////////////////////////////////////////////

  //Get the active fibers from the downstream tracking XBPF
  std::vector<short> firstDownstreamFibers = beamspill->GetActiveFibers(firstDownstreamName, theTrigger);
  std::vector<short> secondDownstreamFibers = beamspill->GetActiveFibers(secondDownstreamName, theTrigger);

  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << firstDownstreamName << " has " << firstDownstreamFibers.size() << " active fibers" << "\n";
    for(size_t i = 0; i < firstDownstreamFibers.size(); ++i){
      MF_LOG_INFO("BeamEvent") << firstDownstreamFibers[i] << " ";
    }
    MF_LOG_INFO("BeamEvent") << "\n";

    MF_LOG_INFO("BeamEvent") << secondDownstreamName << " has " << secondDownstreamFibers.size() << " active fibers" << "\n";
    for(size_t i = 0; i < secondDownstreamFibers.size(); ++i){
      MF_LOG_INFO("BeamEvent") << secondDownstreamFibers[i] << " ";
    }
    MF_LOG_INFO("BeamEvent") << "\n";
  }
  //////////////////////////////////////////////

  if( (firstUpstreamFibers.size() < 1) || (secondUpstreamFibers.size() < 1) || (firstDownstreamFibers.size() < 1) || (secondDownstreamFibers.size() < 1) ){
    if( fPrintDebug )
      MF_LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not making track" << "\n";
    return;
  }

  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< std::pair<size_t, size_t> > upstreamPairedFibers;
  std::vector< std::pair<size_t, size_t> > downstreamPairedFibers;
  std::string firstUpstreamType    = fDeviceTypes[firstUpstreamName]; 
  std::string secondUpstreamType   = fDeviceTypes[secondUpstreamName]; 
  std::string firstDownstreamType  = fDeviceTypes[firstDownstreamName]; 
  std::string secondDownstreamType = fDeviceTypes[secondDownstreamName]; 

  std::vector< TVector3 > upstreamPositions;
  std::vector< TVector3 > downstreamPositions; 
  
  //Pair the upstream fibers together
  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Upstream" << "\n";

  for(size_t iF1 = 0; iF1 < firstUpstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstUpstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondUpstreamFibers.size(); ++iF2){
      size_t secondFiber = secondUpstreamFibers[iF2];

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "Paired: " << firstFiber << " " << secondFiber << "\n"; 

      upstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondUpstreamFibers.size() - 1){
        if (secondUpstreamFibers[iF2] == (secondUpstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstUpstreamFibers.size() - 1){
      if (firstUpstreamFibers[iF1] == (firstUpstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  for(size_t iF = 0; iF < upstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = upstreamPairedFibers.at(iF);
 
    if(firstUpstreamType == "horiz" && secondUpstreamType == "vert"){
      double xPos = GetPosition(firstUpstreamName, thePair.first);
      double yPos = GetPosition(secondUpstreamName, thePair.second);
      
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      upstreamPositions.push_back( posInDet );

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
        MF_LOG_INFO("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      }

    }
    else if(firstUpstreamType == "vert" && secondUpstreamType == "horiz"){
      double yPos = GetPosition(firstUpstreamName, thePair.first);
      double xPos = GetPosition(secondUpstreamName, thePair.second);

      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      upstreamPositions.push_back( posInDet );

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
        MF_LOG_INFO("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      }
    }

  }
  
  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "Downstream" << "\n";

  //Pair the downstream fibers together
  for(size_t iF1 = 0; iF1 < firstDownstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstDownstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondDownstreamFibers.size(); ++iF2){
      size_t secondFiber = secondDownstreamFibers[iF2];

      if( fPrintDebug )
        MF_LOG_INFO("BeamEvent") << "Paired: " << firstFiber << " " << secondFiber << "\n"; 
      downstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondDownstreamFibers.size() - 1){
        if (secondDownstreamFibers[iF2] == (secondDownstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstDownstreamFibers.size() - 1){
      if (firstDownstreamFibers[iF1] == (firstDownstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }

  for(size_t iF = 0; iF < downstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = downstreamPairedFibers.at(iF);
 
    if(firstDownstreamType == "horiz" && secondDownstreamType == "vert"){
      double xPos = GetPosition(firstDownstreamName, thePair.first);
      double yPos = GetPosition(secondDownstreamName, thePair.second);

      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      downstreamPositions.push_back( posInDet );

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
        MF_LOG_INFO("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      }
    }
    else if(firstDownstreamType == "vert" && secondDownstreamType == "horiz"){
      double yPos = GetPosition(firstDownstreamName, thePair.first);
      double xPos = GetPosition(secondDownstreamName, thePair.second);

      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      downstreamPositions.push_back( posInDet );

      if( fPrintDebug ){
        MF_LOG_INFO("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
        MF_LOG_INFO("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      }

    }

  }
 

  for(size_t iU = 0; iU < upstreamPositions.size(); ++iU){
    for(size_t iD = 0; iD < downstreamPositions.size(); ++iD){
      std::vector<TVector3> thePoints;
      thePoints.push_back(upstreamPositions.at(iU));
      thePoints.push_back(downstreamPositions.at(iD));

      //Now project the last point to the TPC face
      thePoints.push_back( ProjectToTPC(thePoints[0],thePoints[1]) );    
      std::vector<TVector3> theMomenta;

      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );

      recob::Track tempTrack(recob::TrackTrajectory(recob::tracking::convertCollToPoint(thePoints),
                                                                        recob::tracking::convertCollToVector(theMomenta),
                                                                        recob::Track::Flags_t(thePoints.size()), false),
                                                 0, -1., 0, recob::tracking::SMatrixSym55(), recob::tracking::SMatrixSym55(), 1);

      beamevt->AddBeamTrack( tempTrack );
    }
  }

}

void proto::BeamEvent::MomentumSpec(size_t theTrigger){
  
  if( fPrintDebug ) 
    MF_LOG_INFO("BeamEvent") << "Doing momentum spectrometry for trigger " << beamspill->GetT0Sec(theTrigger) << " " << beamspill->GetT0Nano(theTrigger) << "\n";

  double LB = mag_P1*fabs(current[0]);
  double deltaI = fabs(current[0]) - mag_P4;
  if(deltaI>0) LB+= mag_P3*deltaI*deltaI;

  //Get the active fibers from the upstream tracking XBPF
  std::string firstBPROF1Type    = fDeviceTypes[firstBPROF1]; 
  std::string secondBPROF1Type   = fDeviceTypes[secondBPROF1]; 
  std::vector<short> BPROF1Fibers;
  std::string BPROF1Name;

  if (firstBPROF1Type == "horiz" && secondBPROF1Type == "vert"){
    BPROF1Fibers = beamspill->GetActiveFibers(firstBPROF1, theTrigger);
    BPROF1Name = firstBPROF1;

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << firstBPROF1 << " has " << BPROF1Fibers.size() << " active fibers" << "\n";
      for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
        MF_LOG_INFO("BeamEvent") << BPROF1Fibers[i] << " ";
      }
      MF_LOG_INFO("BeamEvent") << "\n";
    }

  }
  else if(secondBPROF1Type == "horiz" && firstBPROF1Type == "vert"){
    BPROF1Fibers = beamspill->GetActiveFibers(secondBPROF1, theTrigger);
    BPROF1Name = secondBPROF1;

    if( fPrintDebug ){
      MF_LOG_INFO("BeamEvent") << secondBPROF1 << " has " << BPROF1Fibers.size() << " active fibers" << "\n";
      for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
        MF_LOG_INFO("BeamEvent") << BPROF1Fibers[i] << " ";
      }
      MF_LOG_INFO("BeamEvent") << "\n";
    }
  }
  else{
    MF_LOG_WARNING("BeamEvent") << "Error: type is not correct" << "\n";
    return;
  }
  //////////////////////////////////////////////

  if( (BPROF1Fibers.size() < 1) ){
    if( fPrintDebug )
      MF_LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not checking momentum" << "\n";
    return;
  }
  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< short > strippedFibers; 
  for(size_t iF1 = 0; iF1 < BPROF1Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF1Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[iF1] == (BPROF1Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
  
  double X1, X2, X3;

  //X-direction convention is reversed for the spectrometer   
  X1 = -1.*GetPosition(BPROF1Name, strippedFibers[0])/1.E3;

  //BPROF2////
  //
  std::vector<short>  BPROF2Fibers = beamspill->GetActiveFibers(BPROF2, theTrigger);
  if( (BPROF2Fibers.size() < 1) ){
    if( fPrintDebug )
      MF_LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not checking momentum" << "\n";
    return;
  }

  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << BPROF2 << " has " << BPROF2Fibers.size() << " active fibers at time " << beamspill->GetFiberTime(BPROF2,theTrigger) << "\n";
    for(size_t i = 0; i < BPROF2Fibers.size(); ++i){
      MF_LOG_INFO("BeamEvent") << BPROF2Fibers[i] << " ";
    }
    MF_LOG_INFO("BeamEvent") << "\n";
  }

  strippedFibers.clear(); 

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "BPROF2" << "\n";

  for(size_t iF1 = 0; iF1 < BPROF2Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF2Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF2Fibers.size() - 1){
      if (BPROF2Fibers[iF1] == (BPROF2Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X2 = -1.*GetPosition(BPROF2, strippedFibers[0])/1.E3; 
  ////////////

  //BPROF3////
  //
  std::vector<short>  BPROF3Fibers = beamspill->GetActiveFibers(BPROF3, theTrigger);
  if( (BPROF3Fibers.size() < 1) ){
    if( fPrintDebug )
      MF_LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not checking momentum" << "\n";
    return;
  }

  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << BPROF3 << " has " << BPROF3Fibers.size() << " active fibers at time " << beamspill->GetFiberTime(BPROF3,theTrigger) << "\n";
    for(size_t i = 0; i < BPROF3Fibers.size(); ++i){
      MF_LOG_INFO("BeamEvent") << BPROF3Fibers[i] << " ";
    }
    MF_LOG_INFO("BeamEvent") << "\n";
  }

  strippedFibers.clear(); 

  if( fPrintDebug )
    MF_LOG_INFO("BeamEvent") << "BPROF3" << "\n";

  for(size_t iF1 = 0; iF1 < BPROF3Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF3Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF3Fibers.size() - 1){
      if (BPROF3Fibers[iF1] == (BPROF3Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X3 = -1.*GetPosition(BPROF3, strippedFibers[0])/1.E3; 
  ////////////
  
  if( (BPROF1Fibers.size() == 1) && (BPROF2Fibers.size() == 1) && (BPROF3Fibers.size() == 1) ){
    //Calibrate the positions
    //-1.*( FiberPos ) -> -1.*( FiberPos + ShiftDist )
    // = -1.*FiberPos - ShiftDist
    X1 = X1 - fBProf1Shift*1.e-3; 
    X2 = X2 - fBProf2Shift*1.e-3; 
    X3 = X3 - fBProf3Shift*1.e-3; 
    double cosTheta = MomentumCosTheta(X1,X2,X3);
    double momentum = 299792458*LB/(1.E9 * acos(cosTheta));


    if( fDebugMomentum ){
      fCutMomentum->Fill(momentum);
      MF_LOG_INFO("BeamEvent") << "Filling Cut Momentum Spectrum" << "\n";
    }

  }


  if( fPrintDebug ){
    MF_LOG_INFO("BeamEvent") << "Getting all trio-wise hits" << "\n";
    MF_LOG_INFO("BeamEvent") << "N1,N2,N3 " << BPROF1Fibers.size()
              << " "         << BPROF2Fibers.size() 
              << " "         << BPROF3Fibers.size() << "\n";
  }

  for(size_t i1 = 0; i1 < BPROF1Fibers.size(); ++i1){
    double x1,x2,x3;

    x1 = -1.*GetPosition(BPROF1Name, BPROF1Fibers[i1])/1.E3;
    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Add .5 mm
        x1 += .0005;
      }
    }

    for(size_t i2 = 0; i2 < BPROF2Fibers.size(); ++i2){
      x2 = -1.*GetPosition(BPROF2, BPROF2Fibers[i2])/1.E3;
      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Add .5 mm
          x2 += .0005;
        }
      }

      for(size_t i3 = 0; i3 < BPROF3Fibers.size(); ++i3){
        
        if( fPrintDebug )
          MF_LOG_INFO("BeamEvent") << "\t" << i1 << " " << i2 << " " << i3 << "\n";

        x3 = -1.*GetPosition(BPROF3, BPROF3Fibers[i3])/1.E3;
        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Add .5 mm
            x3 += .0005;
          }
        }


        //Calibrate the positions
        //-1.*( FiberPos ) -> -1.*( FiberPos + ShiftDist )
        // = -1.*FiberPos - ShiftDist
        x1 = x1 - fBProf1Shift*1.e-3; 
        x2 = x2 - fBProf2Shift*1.e-3; 
        x3 = x3 - fBProf3Shift*1.e-3; 

        double cosTheta_full = MomentumCosTheta(x1,x2,x3);        
        double momentum_full = 299792458*LB/(1.E9 * acos(cosTheta_full));
        if( fDebugMomentum ) fFullMomentum->Fill(momentum_full);

        beamevt->AddRecoBeamMomentum(momentum_full);

        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Skip the next
            ++i3;
          }
        }
        
      }

      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Skip the next
          ++i2;
        }
      }
    }

    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Skip the next
        ++i1;
      }
    }
  }
}

double proto::BeamEvent::MomentumCosTheta(double X1, double X2, double X3){
  double a =  (X2*L3 - X3*L2)*cos(fBeamBend)/(L3-L2);

 
  double numTerm = (a - X1)*( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ) + L1*( L3 - L2 );

  double denomTerm1, denomTerm2, denom;
  denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
  denomTerm2 = sqrt( TMath::Power( ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ),2)
                   + TMath::Power( ( (L3 - L2) ),2) );
  denom = denomTerm1 * denomTerm2;

  double cosTheta = numTerm/denom;  
  return cosTheta;
}

TVector3 proto::BeamEvent::ProjectToTPC(TVector3 firstPoint, TVector3 secondPoint){
  TVector3 dR = (secondPoint - firstPoint);
  
  double deltaZ = -1.*secondPoint.Z();
  double deltaX = deltaZ * (dR.X() / dR.Z());
  double deltaY = deltaZ * (dR.Y() / dR.Z());

  TVector3 lastPoint = secondPoint + TVector3(deltaX, deltaY, deltaZ);
  return lastPoint;
}

double proto::BeamEvent::GetPosition(std::string deviceName, int fiberIdx){
  if(fiberIdx > 192){ MF_LOG_WARNING("BeamEvent") << "Please select fiber in range [0,191]" << "\n"; return -1.;}
  double size = fFiberDimension[deviceName];
  
  //Define 0th fiber as farthest positive. Last fiber is farthest negative. Center is between 96 and 97 
  double pos = size*(96 - fiberIdx) - size/2.;
  return pos;
}


DEFINE_ART_MODULE(proto::BeamEvent)
