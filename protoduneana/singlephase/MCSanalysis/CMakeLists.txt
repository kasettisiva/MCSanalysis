art_make(
LIBRARY_NAME
	MCSanalysis
LIB_LIBRARIES
	# Root
	${ROOT_BASIC_LIB_LIST}

	# Framework (not sure if all are needed, but this works and don't want to touch it)
  	${ART_FRAMEWORK_CORE}
  	${ART_FRAMEWORK_PRINCIPAL}
  	${ART_FRAMEWORK_SERVICES_REGISTRY}
  	${ART_FRAMEWORK_SERVICES_BASIC}
  	${ART_FRAMEWORK_SERVICES_OPTIONAL}
  	${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  	${ART_PERSISTENCY_COMMON}
  	${ART_UTILITIES}
	art_root_io
	${ART_ROOT_IO_TFILE_SUPPORT}
	${ART_ROOT_IO_TFILESERVICE_SERVICE}
  	art_Persistency_Provenance
	${MF_MESSAGELOGGER}
  	canvas
  	cetlib
  	cetlib_except

	# LArSoft reco
	lardataobj_RecoBase
	larreco_RecoAlg
	larcorealg_Geometry

	# LArSoft Truth
	nusimdata_SimulationBase
	larsim_MCCheater_BackTrackerService_service

MODULE_LIBRARIES
	MCSanalysis
)

install_headers()
install_fhicl()
install_source()
