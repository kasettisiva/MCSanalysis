#include "ProtoDUNEBeamCuts.h"
#include "dune/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

protoana::ProtoDUNEBeamCuts::ProtoDUNEBeamCuts( const fhicl::ParameterSet & pset ){


  std::vector< fhicl::ParameterSet > DataCuts_sets = pset.get< std::vector< fhicl::ParameterSet > >
                                                             ("DataCuts");
  std::vector< fhicl::ParameterSet > MCCuts_sets   = pset.get< std::vector< fhicl::ParameterSet > >
                                                             ("MCCuts");

  //Go through the data cuts 
  for( size_t i = 0; i < DataCuts_sets.size(); ++i ){
    std::string momentum = DataCuts_sets[i].get<std::string>("Momentum");
    if( std::find( valid_momenta.begin(), valid_momenta.end(), momentum ) == valid_momenta.end() ){
      std::exception e;
      throw e;
    }

    if( DataCuts.find( momentum ) != DataCuts.end() ){
      std::cerr << "Error: Attempting to duplicate a cut parameter set" << std::endl;
      std::exception e;
      throw e;
    }
     
    DataCuts[momentum] = DataCuts_sets[i];
  }


  //Now go through the MC cuts
  for( size_t i = 0; i < MCCuts_sets.size(); ++i ){
    
    std::string momentum = DataCuts_sets[i].get<std::string>("Momentum");
    if( std::find( valid_momenta.begin(), valid_momenta.end(), momentum ) == valid_momenta.end() ){
      std::exception e;
      throw e;
    }

    if( MCCuts.find( momentum ) != MCCuts.end() ){
      std::cerr << "Error: Attempting to duplicate a cut parameter set" << std::endl;
      std::exception e;
      throw e;
    }
     
    MCCuts[momentum] = MCCuts_sets[i];
  }


}

bool protoana::ProtoDUNEBeamCuts::IsBeamlike( const recob::Track & track, art::Event const & evt, std::string momentum ){

   if( std::find( valid_momenta.begin(), valid_momenta.end(), momentum ) == valid_momenta.end() ){
     std::cerr << "Error. Momentum provided not in range" << std::endl;
     std::exception e;
     throw e;
   }

   fhicl::ParameterSet BeamPars = ( evt.isRealData() ? DataCuts[momentum] : MCCuts[momentum] );   
   BeamVals theBeamVals = ( evt.isRealData() ? GetDataBeam( evt ) : GetMCBeam( evt ) );
   if( !theBeamVals.Valid ){ 
      return false;
   }

   double startX = track.Trajectory().Start().X();
   double startY = track.Trajectory().Start().Y();
   double startZ = track.Trajectory().Start().Z();

   double endX = track.Trajectory().End().X();
   double endY = track.Trajectory().End().Y();
   double endZ = track.Trajectory().End().Z();

   auto startDir = track.StartDirection();
   auto endDir   = track.EndDirection();
   double trackDirX = 0.;
   double trackDirY = 0.;
   double trackDirZ = 0.;

   //'Flip' the track if endZ < startZ
   if( endZ < startZ ){
     startX = endX;
     startY = endY;
     startZ = endZ;

     trackDirX = -1. * endDir.X();
     trackDirY = -1. * endDir.Y();
     trackDirZ = -1. * endDir.Z();
   }
   else{
     trackDirX = startDir.X();
     trackDirY = startDir.Y();
     trackDirZ = startDir.Z();
   }

   double deltaX = startX - theBeamVals.X;
   double deltaY = startY - theBeamVals.Y;

   double costheta = theBeamVals.DirX*trackDirX + theBeamVals.DirY*trackDirY + theBeamVals.DirZ*trackDirZ;

   std::pair< double, double > startX_cut = BeamPars.get< std::pair< double, double > >("TrackStartXCut");
   std::pair< double, double > startY_cut = BeamPars.get< std::pair< double, double > >("TrackStartYCut");
   std::pair< double, double > startZ_cut = BeamPars.get< std::pair< double, double > >("TrackStartZCut");
   double costheta_cut = BeamPars.get< double >("TrackDirCut");

   static bool first = true;
   if (first){
     first = false;
     std::cout<<"Beam momentum = "<<momentum<<std::endl;
     std::cout<<"Cut on deltaX: ["<<startX_cut.first<<","<<startX_cut.second<<"]"<<std::endl;
     std::cout<<"Cut on deltaY: ["<<startY_cut.first<<","<<startY_cut.second<<"]"<<std::endl;
     std::cout<<"Cut on deltaZ: ["<<startZ_cut.first<<","<<startZ_cut.second<<"]"<<std::endl;
   }

   if( deltaX < startX_cut.first || deltaX > startX_cut.second )
     return false;
   if( deltaY < startY_cut.first || deltaY > startY_cut.second )
     return false;
   if( startZ < startZ_cut.first || startZ > startZ_cut.second )
     return false;
   if( costheta < costheta_cut )
     return false;   

   //If here, the track is in the good beam region
   return true;

}

bool protoana::ProtoDUNEBeamCuts::IsBeamlike( const recob::Shower & shower, art::Event const & evt, std::string momentum ){

   if( std::find( valid_momenta.begin(), valid_momenta.end(), momentum ) == valid_momenta.end() ){
     std::cerr << "Error. Momentum provided not in range" << std::endl;
     std::exception e;
     throw e;
   }

   fhicl::ParameterSet BeamPars = ( evt.isRealData() ? DataCuts[momentum] : MCCuts[momentum] );   
   BeamVals theBeamVals = ( evt.isRealData() ? GetDataBeam( evt ) : GetMCBeam( evt ) );
   if( !theBeamVals.Valid ){ 
      return false;
   }

   double startX = shower.ShowerStart().X();
   double startY = shower.ShowerStart().Y();
   double startZ = shower.ShowerStart().Z();

   auto startDir = shower.Direction();

   //flip the direction if pointing backward
   double flip = ( startDir.Z() < 0. ? -1. : 1. );

   double showerDirX = flip * startDir.X();
   double showerDirY = flip * startDir.Y();
   double showerDirZ = flip * startDir.Z();

   double deltaX = startX - theBeamVals.X;
   double deltaY = startY - theBeamVals.Y;

   double costheta = theBeamVals.DirX*showerDirX + theBeamVals.DirY*showerDirY + theBeamVals.DirZ*showerDirZ;
   std::pair< double, double > startX_cut = BeamPars.get< std::pair< double, double > >("TrackStartXCut");
   std::pair< double, double > startY_cut = BeamPars.get< std::pair< double, double > >("TrackStartYCut");
   std::pair< double, double > startZ_cut = BeamPars.get< std::pair< double, double > >("TrackStartZCut");
   double costheta_cut = BeamPars.get< double >("TrackDirCut");

   if( deltaX < startX_cut.first || deltaX > startX_cut.second )
     return false;
   if( deltaY < startY_cut.first || deltaY > startY_cut.second )
     return false;
   if( startZ < startZ_cut.first || startZ > startZ_cut.second )
     return false;
   if( costheta < costheta_cut )
     return false;   

   //If here, the track is in the good beam region
   return true;

}



protoana::BeamVals protoana::ProtoDUNEBeamCuts::GetMCBeam( const art::Event & evt ){

  protoana::ProtoDUNETruthUtils truthUtil;
  auto mcTruths = evt.getValidHandle< std::vector< simb::MCTruth > >("generator");
  const simb::MCParticle* true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

  BeamVals result;
  result.Valid = false;

  if( !true_beam_particle ){
    std::cout << "No true beam particle" << std::endl;       
    return result;
  }

  result.Valid = true;
  result.DirX = true_beam_particle->Px() / true_beam_particle->P();
  result.DirY = true_beam_particle->Py() / true_beam_particle->P();
  result.DirZ = true_beam_particle->Pz() / true_beam_particle->P();

  //Project the beam to Z = 0
  result.X = true_beam_particle->Position(0).X() + (-1.*true_beam_particle->Position(0).Z())*(result.DirX / result.DirZ);
  result.Y = true_beam_particle->Position(0).Y() + (-1.*true_beam_particle->Position(0).Z())*(result.DirY / result.DirZ);
  return result;
}

protoana::BeamVals protoana::ProtoDUNEBeamCuts::GetDataBeam( const art::Event & evt ){
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  auto beamHandle = evt.getValidHandle< std::vector< beam::ProtoDUNEBeamEvent > >("beamevent");
                        
  BeamVals result;
  result.Valid = false;

  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }
  //Should just have one
  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0));

  const std::vector< recob::Track > & beamTracks = beamEvent.GetBeamTracks();
  if( beamTracks.size() == 0 ){
    std::cout << "Warning: no tracks associated to beam data" << std::endl;
    return result;  
  }
  else if( beamTracks.size() > 1 ){
    std::cout << "Warning: mutiple tracks associated to beam data" << std::endl;
    return result;  
  }
  
  result.Valid = true;
  result.X = beamTracks.at(0).Trajectory().End().X();
  result.Y = beamTracks.at(0).Trajectory().End().Y();

  result.DirX = beamTracks.at(0).EndDirection().X();
  result.DirY = beamTracks.at(0).EndDirection().Y();
  result.DirZ = beamTracks.at(0).EndDirection().Z();
  return result;
}
