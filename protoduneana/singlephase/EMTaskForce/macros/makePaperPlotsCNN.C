void makePaperPlotsCNN(){

  // Protons  : {2212}
  // Pions    : {211,13} (data gives 13 for all muons / pions at 1 GeV
  // Positrons: {11}
  const std::vector<std::vector<short>> pdgCode = {{2212},{211,13},{11}};
  std::map<const short,std::string> particleNames;
  particleNames.insert(std::make_pair(2,"positron"));
  particleNames.insert(std::make_pair(1,"pion"));
  particleNames.insert(std::make_pair(0,"proton"));
  std::map<const short,std::string> particleNamesCap;
  particleNamesCap.insert(std::make_pair(2,"Positrons"));
  particleNamesCap.insert(std::make_pair(1,"Pions"));
  particleNamesCap.insert(std::make_pair(0,"Protons"));

  gStyle->SetOptStat(0);

//  TChain *c = new TChain("cnnana/ftree");
  TChain *c = new TChain("cnnana1/ftree");

  // We can divide things up based on run numbers, so just use one chain
  c->Add("checkcnn_data_5809_production4_feb2021.root"); // Prod 4 Data
  c->Add("checkcnn_data_5387_production4_feb2021.root"); // Prod 4 Data
  c->Add("checkcnn_mc_production4_feb2021.root"); // Prod 4 MC

  int fRun = 0;
  double fParticleEMScore = 0.;
  std::vector<short> *fTPC = 0x0;
  std::vector<short> *fPlane = 0x0;
  std::vector<short> *fWire = 0x0;
  std::vector<double> *fScoreEM = 0x0;
  int fBeamPDG = 0;

  TBranch *rootB=0;
  c->SetBranchAddress("run",&fRun);
  c->SetBranchAddress("average_score_em",&fParticleEMScore);
  c->SetBranchAddress("wire",&fWire,&rootB);
  c->SetBranchAddress("score_em",&fScoreEM,&rootB);
  c->SetBranchAddress("tpc",&fTPC,&rootB);
  c->SetBranchAddress("plane",&fPlane,&rootB);
  c->SetBranchAddress("beampdg",&fBeamPDG);
  

  // Plots for hit-level distributions
  const int nbins = 100;

  std::vector<TH1F*> hNoCutMC;
  hNoCutMC.push_back(new TH1F("hNoCutMC_Proton",";CNN EM Score",nbins,0,1));
  hNoCutMC.push_back(new TH1F("hNoCutMC_Pion",";CNN EM Score",nbins,0,1));
  hNoCutMC.push_back(new TH1F("hNoCutMC_Positron",";CNN EM Score",nbins,0,1));

  std::vector<TH1F*> hNoCutData;
  hNoCutData.push_back(new TH1F("hNoCutData_Proton",";CNN EM Score",nbins,0,1));
  hNoCutData.push_back(new TH1F("hNoCutData_Pion",";CNN EM Score",nbins,0,1));
  hNoCutData.push_back(new TH1F("hNoCutData_Positron",";CNN EM Score",nbins,0,1));

  std::vector<TH1F*> h90MC;
  h90MC.push_back(new TH1F("h90MC_Proton",";CNN EM Score",nbins,0,1));
  h90MC.push_back(new TH1F("h90MC_Pion",";CNN EM Score",nbins,0,1));
  h90MC.push_back(new TH1F("h90MC_Positron",";CNN EM Score",nbins,0,1));

  std::vector<TH1F*> h90Data;
  h90Data.push_back(new TH1F("h90Data_Proton",";CNN EM Score",nbins,0,1));
  h90Data.push_back(new TH1F("h90Data_Pion",";CNN EM Score",nbins,0,1));
  h90Data.push_back(new TH1F("h90Data_Positron",";CNN EM Score",nbins,0,1));

  // Plots for particle-level distributions

  const int nbinsP = 50;

  std::vector<TH1F*> hParticleScoreMC;
  hParticleScoreMC.push_back(new TH1F("hParticleScoreMC_Proton",";Average CNN EM Score",nbinsP,0,1));
  hParticleScoreMC.push_back(new TH1F("hParticleScoreMC_Pion",";Average CNN EM Score",nbinsP,0,1));
  hParticleScoreMC.push_back(new TH1F("hParticleScoreMC_Positron",";Average CNN EM Score",nbinsP,0,1));

  std::vector<TH1F*> hParticleScoreData;
  hParticleScoreData.push_back(new TH1F("hParticleScoreData_Proton",";Average CNN EM Score",nbinsP,0,1));
  hParticleScoreData.push_back(new TH1F("hParticleScoreData_Pion",";Average CNN EM Score",nbinsP,0,1));
  hParticleScoreData.push_back(new TH1F("hParticleScoreData_Positron",";Average CNN EM Score",nbinsP,0,1));

  const short wireCut = 90;

  for(unsigned int e = 0; e < c->GetEntries(); ++e){

    c->GetEntry(e);

    unsigned int particleIndex = 999;
    if(std::abs(fBeamPDG) == 2212) particleIndex = 0;
    else if(std::abs(fBeamPDG) == 211 || std::abs(fBeamPDG) == 13) particleIndex = 1;
    else if(std::abs(fBeamPDG) == 11) particleIndex = 2;
    else continue;
    
    if(fRun != 5809 && fRun != 5387) hParticleScoreMC.at(particleIndex)->Fill(fParticleEMScore); 
    else hParticleScoreData.at(particleIndex)->Fill(fParticleEMScore);

    for(unsigned int h = 0; h < fWire->size(); ++h){

      if(fTPC->at(h) != 1) continue;
      if(fPlane->at(h) != 2) continue;

      // MC
      if(fRun != 5809 && fRun != 5387){
        hNoCutMC.at(particleIndex)->Fill(fScoreEM->at(h));
        if(fWire->at(h) >= wireCut) h90MC.at(particleIndex)->Fill(fScoreEM->at(h));
      }
      // Data
      else{
        hNoCutData.at(particleIndex)->Fill(fScoreEM->at(h));
        if(fWire->at(h) >= wireCut) h90Data.at(particleIndex)->Fill(fScoreEM->at(h));
      }
    }
  }

  std::cout << "MC particle totals  : " << hParticleScoreMC.at(0)->GetEntries() << ", " << hParticleScoreMC.at(1)->GetEntries() << ", " << hParticleScoreMC.at(2)->GetEntries() << std::endl;
  std::cout << "Data particle totals: " << hParticleScoreData.at(0)->GetEntries() << ", " << hParticleScoreData.at(1)->GetEntries() << ", " << hParticleScoreData.at(2)->GetEntries() << std::endl;

  std::cout << "MC hit totals  : " << h90MC.at(0)->GetEntries() << ", " << h90MC.at(1)->GetEntries() << ", " << h90MC.at(2)->GetEntries() << std::endl;
  std::cout << "Data hit totals: " << h90Data.at(0)->GetEntries() << ", " << h90Data.at(1)->GetEntries() << ", " << h90Data.at(2)->GetEntries() << std::endl;

  // Totals
  std::vector<unsigned int> nNoCutMCEvents = {0,0,0};
  std::vector<unsigned int> n90MCEvents = {0,0,0};
  std::vector<unsigned int> nNoCutDataEvents = {0,0,0};
  std::vector<unsigned int> n90DataEvents = {0,0,0};

  std::vector<unsigned int> nNoCutMCCorrect = {0,0,0};
  std::vector<unsigned int> n90MCCorrect = {0,0,0};
  std::vector<unsigned int> nNoCutDataCorrect = {0,0,0};
  std::vector<unsigned int> n90DataCorrect = {0,0,0};

  const unsigned int cutValue = 72;
  for(unsigned int h = 0; h < hNoCutMC.size(); ++h){
    // Above a cut value of 0.72
    unsigned int lowerBin = 0;
    unsigned int upperBin = cutValue - 1;
    if(h == 2){
      lowerBin = cutValue;
      upperBin = 101;
    }

    // Get the total number of events
    nNoCutMCEvents.at(h) = hNoCutMC.at(h)->Integral();
    n90MCEvents.at(h)    = h90MC.at(h)->Integral();
    nNoCutDataEvents.at(h) = hNoCutData.at(h)->Integral();
    n90DataEvents.at(h)    = h90Data.at(h)->Integral();
    // Get the number of correctly identified events
    nNoCutMCCorrect.at(h)   = hNoCutMC.at(h)->Integral(lowerBin,upperBin);
    n90MCCorrect.at(h)      = h90MC.at(h)->Integral(lowerBin,upperBin);
    nNoCutDataCorrect.at(h) = hNoCutData.at(h)->Integral(lowerBin,upperBin);
    n90DataCorrect.at(h)    = h90Data.at(h)->Integral(lowerBin,upperBin);

    // Normalise the histograms now
    hNoCutMC.at(h)->Sumw2();
    hNoCutMC.at(h)->Scale(1.0 / static_cast<float>(hNoCutMC.at(h)->Integral()));
    h90MC.at(h)->Sumw2();
    h90MC.at(h)->Scale(1.0 / static_cast<float>(h90MC.at(h)->Integral()));
    hNoCutData.at(h)->Sumw2();
    hNoCutData.at(h)->Scale(1.0 / static_cast<float>(hNoCutData.at(h)->Integral()));
    h90Data.at(h)->Sumw2();
    h90Data.at(h)->Scale(1.0 / static_cast<float>(h90Data.at(h)->Integral()));

    // Set the line colours and scales etc
    hNoCutMC.at(h)->SetMinimum(1e-5);
    hNoCutMC.at(h)->SetMaximum(1e0);
    hNoCutMC.at(h)->SetLineColor(kRed);
    hNoCutMC.at(h)->SetLineWidth(2);
    h90MC.at(h)->SetMinimum(1e-5);
    h90MC.at(h)->SetMaximum(1e0);
    h90MC.at(h)->SetLineColor(kRed);
    h90MC.at(h)->SetLineWidth(2);
    hNoCutData.at(h)->SetLineColor(kBlack);
    hNoCutData.at(h)->SetLineWidth(2);
    h90Data.at(h)->SetLineColor(kBlack);
    h90Data.at(h)->SetLineWidth(2);

    // Same for particle level plots

    hParticleScoreMC.at(h)->Sumw2();
    hParticleScoreMC.at(h)->Scale(1.0 / static_cast<float>(hParticleScoreMC.at(h)->Integral()));
    hParticleScoreData.at(h)->Sumw2();
    hParticleScoreData.at(h)->Scale(1.0 / static_cast<float>(hParticleScoreData.at(h)->Integral()));

    hParticleScoreMC.at(h)->SetMinimum(1e-5);
    hParticleScoreMC.at(h)->SetMaximum(1.0);
    hParticleScoreMC.at(h)->SetLineColor(kRed);
    hParticleScoreMC.at(h)->SetLineWidth(2);
    hParticleScoreData.at(h)->SetLineColor(kBlack);
    hParticleScoreData.at(h)->SetLineWidth(2);

  }

  TLegend *leg = 0x0; 
  TCanvas *can = new TCanvas("can","",0,0,800,600);
  TLatex tt;
  tt.SetNDC();
  tt.SetTextSize(0.04);
  for(unsigned int pIndex = 0; pIndex < 3; ++pIndex){
    leg = new TLegend(0.3,0.65,0.5,0.85,("1.0 GeV "+particleNamesCap.at(pIndex)).c_str());
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hNoCutData.at(pIndex),"Data","lp");
    leg->AddEntry(hNoCutMC.at(pIndex),"Simulation","l"); 

    can->SetLogy(1);
    hNoCutMC.at(pIndex)->Draw("hist");
    hNoCutData.at(pIndex)->Draw("same");
    leg->Draw("same");
    tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP} Preliminary");
//    const std::string plotBaseName = "finalPlots/default/emscore_"+particleNames.at(pIndex);
    const std::string plotBaseName = "finalPlots/emscore_"+particleNames.at(pIndex);
    can->Print((plotBaseName+"_nocut_log.C").c_str());
    can->Print((plotBaseName+"_nocut_log.png").c_str());
    can->Print((plotBaseName+"_nocut_log.pdf").c_str());

    h90MC.at(pIndex)->Draw("hist");
    h90Data.at(pIndex)->Draw("same");
    leg->Draw("same");
    tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP} Preliminary");
    can->Print((plotBaseName+"_90_log.C").c_str());
    can->Print((plotBaseName+"_90_log.png").c_str());
    can->Print((plotBaseName+"_90_log.pdf").c_str());

    const float effNoCutMC = nNoCutMCCorrect.at(pIndex) / static_cast<float>(nNoCutMCEvents.at(pIndex));
    const float errNoCutMC = std::sqrt((1-effNoCutMC) / nNoCutMCEvents.at(pIndex));
    const float eff90MC = n90MCCorrect.at(pIndex) / static_cast<float>(n90MCEvents.at(pIndex));
    const float err90MC = std::sqrt((1-eff90MC) / n90MCEvents.at(pIndex));
    const float effNoCutData = nNoCutDataCorrect.at(pIndex) / static_cast<float>(nNoCutDataEvents.at(pIndex));
    const float errNoCutData = std::sqrt((1-effNoCutData) / nNoCutDataEvents.at(pIndex));
    const float eff90Data = n90DataCorrect.at(pIndex) / static_cast<float>(n90DataEvents.at(pIndex));
    const float err90Data = std::sqrt((1-eff90Data) / n90DataEvents.at(pIndex));
  
    std::cout << "Summary for " << particleNames.at(pIndex) << " candidates:" << std::endl;
    std::cout << " - Fraction correct MC   = " << effNoCutMC << " +/- " << errNoCutMC << " :: " << eff90MC << " +/- " << err90MC << std::endl;
    std::cout << " - Fraction correct Data = " << effNoCutData << " +/- " << errNoCutData << " :: " << eff90Data << " +/- " << err90Data << std::endl;

    leg->SetX1NDC(0.4);
    leg->SetX2NDC(0.6);

    // Now for drawing the particle level plots
    hParticleScoreMC.at(pIndex)->Draw("hist");
    hParticleScoreData.at(pIndex)->Draw("same");
    leg->Draw("same");
    tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP} Preliminary");
    can->Print((plotBaseName+"_per_particle_log.C").c_str());
    can->Print((plotBaseName+"_per_particle_log.png").c_str());
    can->Print((plotBaseName+"_per_particle_log.pdf").c_str());
    
    can->SetLogy(0);
    hParticleScoreMC.at(pIndex)->SetMinimum(0.0);    
    hParticleScoreMC.at(pIndex)->SetMaximum(0.3);    
    hParticleScoreMC.at(pIndex)->Draw("hist");
    hParticleScoreData.at(pIndex)->Draw("same");
    leg->Draw("same");
    tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP} Preliminary");
    can->Print((plotBaseName+"_per_particle.C").c_str());
    can->Print((plotBaseName+"_per_particle.png").c_str());
    can->Print((plotBaseName+"_per_particle.pdf").c_str()); 

    delete leg;
    leg = 0x0;
  }

  

}


