find_ups_product( dunetpc )
art_make(BASENAME_ONLY
  MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${ART_PERSISTENCY_COMMON}
  canvas
  art_Persistency_Provenance
  ${ART_UTILITIES}
  ${MF_MESSAGELOGGER}
  cetlib cetlib_except
  lardataobj_RecoBase
  lardata_Utilities
  lardata_ArtDataHelper
  larcorealg_Geometry
  nusimdata_SimulationBase
  larsim_MCCheater_PhotonBackTrackerService_service           
  larsim_MCCheater_BackTrackerService_service           
  larsim_MCCheater_ParticleInventoryService_service
  larreco_RecoAlg
  ${ROOT_BASIC_LIB_LIST}
  ProtoDUNEDataUtils
  )
install_headers()
install_fhicl()
install_source()

