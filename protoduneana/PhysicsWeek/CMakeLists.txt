simple_plugin(MichelHitCounter "module"
              larcorealg_Geometry
              larreco_Calorimetry
              larreco_RecoAlg_PMAlg
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              lardataobj_RecoBase
              lardata_ArtDataHelper
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ${CETLIB}
              ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(VtxScore "module"
              larcorealg_Geometry
              larreco_Calorimetry
              larreco_RecoAlg_PMAlg
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              lardataobj_RecoBase
              lardata_ArtDataHelper
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ${CETLIB}
              ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)


simple_plugin(MichelEventSelection "module"
              larcorealg_Geometry
              larreco_Calorimetry
              larreco_RecoAlg_PMAlg
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              lardataobj_RecoBase
              lardata_ArtDataHelper
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ${CETLIB}
              ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
