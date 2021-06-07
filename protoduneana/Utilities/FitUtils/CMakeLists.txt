include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

FILE (GLOB src_files *.cxx)

art_make_library( LIBRARY_NAME protodune_fitutils
		  SOURCE ${src_files}
		  LIBRARIES ${ART_FRAMEWORK_CORE}
	       		    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_ROOT_IO_TFILE_SUPPORT} 
			    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			    ${ART_PERSISTENCY_COMMON}
			    canvas
			    art_Persistency_Provenance
			    cetlib cetlib_except
			    ${ART_UTILITIES}
			    ${MF_MESSAGELOGGER}
			    ${FHICLCPP}
			    ${CETLIB}
			    ${ROOT_BASIC_LIB_LIST}
                            
			    ${ROOT_ROOFITCORE}
			    ${ROOT_ROOFIT}
			    ${ROOT_ROOSTATS}
			    ${ROOT_HISTFACTORY}
			    ${ROOT_FOAM}
			    ${ROOT_GEOM}
			    ${ROOT_MINUIT}
			    ${ROOT_MINUIT2}
			    ${ROOT_CORE}
			    ${ROOT_MATHCORE}
                            ProtoDUNEUtilities
	    )
            #message(WARNING ${ROOT_ROOTDATAFRAME})

art_make_exec( NAME runProtoDUNEFit
	       SOURCE runProtoDUNEFit.cc
	       LIBRARIES protodune_fitutils
	    )

art_make_exec( NAME runEventSelection
	       SOURCE runEventSelection.cc
	       LIBRARIES protodune_fitutils
                         ${ROOT_ROOTDATAFRAME}
	    )

art_make_exec( NAME Geant4ReweightSysts 
	       SOURCE Geant4ReweightSysts.cc
	       LIBRARIES protodune_fitutils
	    )

art_make_exec( NAME runPDSPThinSliceFit
	       SOURCE runPDSPThinSliceFit.c++
	       LIBRARIES protodune_fitutils
	     )

install_headers()
install_fhicl()
install_source()
