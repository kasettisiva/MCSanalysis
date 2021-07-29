/*
  Functions used in angle distribution fits, sigmaHL, sigmaRMS, etc.
 */

// paras[0] = kappa_a, paras[1] = kappa_c, paras[2] = sigmaRES, paras[3] = gamma_3D
Double_t Sigma(Double_t* x, Double_t* paras) {
  // Fit Parameters
  Double_t kappa_a = paras[0];
  Double_t kappa_c = paras[1];
  Double_t sigmaRES = paras[2];
  Double_t gamma_3D = paras[3];

  // x-value
  Double_t momentum = x[0];
  Double_t momentum2 = momentum*momentum;

  // Permanenetly Fixed Parameters
  Double_t mass = 0.106;
  Double_t mass2 = mass*mass;

  // Calculate sigmaHL
  Double_t kappa = kappa_a / momentum2 + kappa_c;
  Double_t pbetac = momentum2 / pow(mass2 + momentum2, 0.5); // p*beta*c
  Double_t sigmaHL = kappa/pbetac;

  // Calculate sigmaRMS
  Double_t sigmaRMS = std::pow(sigmaHL*sigmaHL + sigmaRES*sigmaRES, 0.5);

  // Scale sigmaRMS (typically for 3D angles)
  Double_t scaled_sigmaRMS = gamma_3D*sigmaRMS;

  return scaled_sigmaRMS;
}

// Get the most general sigma function as a TF1, with set parameter names and initial values.
TF1* GetSigma_function(const  char* name, Double_t fitMin, Double_t fitMax) {

  // Define sigma_function
  TF1* sigma_function = new TF1(name, Sigma, fitMin, fitMax, 4);
  sigma_function->SetLineColor(kGreen);

  // Set Parameter Names
  sigma_function->SetParName(0, "kappa_a");
  sigma_function->SetParName(1, "kappa_c");
  sigma_function->SetParName(2, "sigma_RES");
  sigma_function->SetParName(3, "gamma_3D");

  // Set Initial Parameters, using kappa_a & kappa_c from MicroBooNE, sigmaRES = 3 (MicroBooNE), & gamma_3D = 0.8525 (typical conversion for Half-Gaussian distribution).
  // These correspond to estimates of parameters that would show up for reco3D angles distributions.
  sigma_function->SetParameters(0.105, 11.004, 3, 0.8525);

  // Return sigma_function
  return sigma_function;
}

// TODO: Documentation

// Fit sigmaRMS using kappa_a, kappa_c, and sigmaRES as fit parameters. (raw sigmaRMS fit)
TF1* GetSigmaRMSraw_function(const char* name, Double_t fitMin, Double_t fitMax) {
  TF1* sigmaRMSraw_function = GetSigma_function(name, fitMin, fitMax);

  sigmaRMSraw_function->FixParameter(3, 1); // Fix gamma_3D to 1.

  return sigmaRMSraw_function;
}

TF1* GetSigmaHL_function(const char* name, Double_t fitMin, Double_t fitMax) {
  TF1* sigmaHL_function = GetSigma_function(name, fitMin, fitMax);

  sigmaHL_function->FixParameter(2, 0); // Fix sigmaRES to 0.
  sigmaHL_function->FixParameter(3, 1); // Fix gamma_3D to 1.

  return sigmaHL_function;
}

TF1* GetSigmaHL_3D_function(const char* name, Double_t fitMin, Double_t fitMax, Double_t kappa_a, Double_t kappa_c) {
  TF1* sigmaHL_3D_function = GetSigma_function(name, fitMin, fitMax);

  sigmaHL_3D_function->FixParameter(0, kappa_a); // Fix kappa_a
  sigmaHL_3D_function->FixParameter(1, kappa_c); // Fix kappa_c
  sigmaHL_3D_function->FixParameter(2, 0); // Fix sigmaRES to 0.

  return sigmaHL_3D_function;
}

TF1* GetSigmaRMS_function(const char* name, Double_t fitMin, Double_t fitMax, Double_t kappa_a, Double_t kappa_c) {
  TF1* sigmaRMS_function = GetSigma_function(name, fitMin, fitMax);
  sigmaRMS_function->SetLineColor(kRed);

  sigmaRMS_function->FixParameter(0, kappa_a); // Fix kappa_a
  sigmaRMS_function->FixParameter(1, kappa_c); // Fix kappa_c
  sigmaRMS_function->FixParameter(3, 1.0); // Fix gamma_3D to 1.

  return sigmaRMS_function;
}

TF1* GetSigmaRMS_3D_function(const char* name, Double_t fitMin, Double_t fitMax, Double_t kappa_a, Double_t kappa_c) {
  TF1* sigmaRMS_3D_function = GetSigma_function(name, fitMin, fitMax);
  sigmaRMS_3D_function->SetLineColor(kRed);

  sigmaRMS_3D_function->FixParameter(0, kappa_a); // Fix kappa_a
  sigmaRMS_3D_function->FixParameter(1, kappa_c); // Fix kappa_c

  return sigmaRMS_3D_function;
}

TF1* GetSigmaRMS_3D_function_fixSigmaRES(const char* name, Double_t fitMin, Double_t fitMax, Double_t kappa_a, Double_t kappa_c, Double_t sigmaRES) {
  TF1* sigmaRMS_3D_function = GetSigmaRMS_3D_function(name, fitMin, fitMax, kappa_a, kappa_c);
  
  sigmaRMS_3D_function->FixParameter(2, sigmaRES); // Fix sigmaRES

  return sigmaRMS_3D_function;
}

TF1* GetSigmaRMS_3D_function_fixGamma3D(const char* name, Double_t fitMin, Double_t fitMax, Double_t kappa_a, Double_t kappa_c, Double_t gamma3D) {
  TF1* sigmaRMS_3D_function = GetSigmaRMS_3D_function(name, fitMin, fitMax, kappa_a, kappa_c);

  sigmaRMS_3D_function->FixParameter(3, gamma3D); // Fix gamma3D

  return sigmaRMS_3D_function;
}
