art_make( MODULE_LIBRARIES
                        lardataobj_RawData
			larcorealg_Geometry
                        larcorealg_CoreUtils
			larcore_Geometry_Geometry_service
                        lardataobj_Simulation
                        nusimdata_SimulationBase
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larsim_MCCheater_PhotonBackTrackerService_service
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GENVECTOR}
                        ${ROOT_GEOM}
)

install_headers()
install_fhicl()
install_source()
install_scripts()
