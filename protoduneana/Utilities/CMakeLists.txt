
art_make( LIBRARY_NAME ProtoDUNEUtilities
	  LIB_LIBRARIES
                        dune_ArtSupport
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larreco_RecoAlg
                        larreco_Calorimetry
                        larreco_RecoAlg_PMAlg
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_ArtDataHelper
                        lardata_Utilities
                        larsim_Simulation 
                        lardataobj_Simulation
                        #ProtoDUNEDataUtils
                        dune_DuneObj
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        nug4_ParticleNavigation 
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES
          ProtoDUNEUtilities
)

install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(FitUtils)
