simple_plugin(		HitNormCheck "module"
			lardataobj_RecoBase
                        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

simple_plugin(		CosmicsdQdx "module"
			lardataobj_RecoBase
                        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			ProtoDUNEDataUtils
			ProtoDUNEUtilities
                        BASENAME_ONLY
)

install_fhicl()
install_source()
