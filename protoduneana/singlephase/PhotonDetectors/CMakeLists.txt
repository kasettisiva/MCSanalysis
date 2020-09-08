simple_plugin(PDSPmatch "module"
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune_Geometry
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)

simple_plugin(PDSPmatchmc "module"
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune_Geometry
                        nusimdata_SimulationBase
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)


install_fhicl()
install_source()
install_scripts()
