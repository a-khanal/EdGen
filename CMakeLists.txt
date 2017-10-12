#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(EdGen)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Debug)
set(useROOT true)        # use true or false
include(FindROOT.cmake)

message(STATUS "CLAS6LIB " $ENV{CLAS6LIB})
set(WITHBOS $ENV{CLAS6LIB})
if(WITHBOS)
	include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} $ENV{CLAS6INC})
	set(LINK_DIRECTORIES
	  ${ROOT_LIBRARY_DIR} $ENV{CLAS6LIB} $ENV{CERNLIB}		    
	  )
	file(GLOB sources 
	  ${PROJECT_SOURCE_DIR}/EdGen.cxx
	  ${PROJECT_SOURCE_DIR}/EdInput.cxx
	  ${PROJECT_SOURCE_DIR}/EdModel.cxx
	  ${PROJECT_SOURCE_DIR}/EdHisto.cxx
	  ${PROJECT_SOURCE_DIR}/EdOutputBos.cxx
	  ${PROJECT_SOURCE_DIR}/EdPhysics.cxx
	  ${PROJECT_SOURCE_DIR}/EdGenPhaseSpace.cxx
	  ${PROJECT_SOURCE_DIR}/EdProcess.cxx
	  ${PROJECT_SOURCE_DIR}/EdCrossSection.cxx
	  )
	file(GLOB headers 
	  ${PROJECT_SOURCE_DIR}/EdGen.h
	  ${PROJECT_SOURCE_DIR}/EdInput.h
	  ${PROJECT_SOURCE_DIR}/EdModel.h
	  ${PROJECT_SOURCE_DIR}/EdHisto.h
	  ${PROJECT_SOURCE_DIR}/EdOutputBos.h
	  ${PROJECT_SOURCE_DIR}/EdPhysics.h
	  ${PROJECT_SOURCE_DIR}/EdGenPhaseSpace.h
	  ${PROJECT_SOURCE_DIR}/EdProcess.h
	  ${PROJECT_SOURCE_DIR}/EdCrossSection.h
	  )
	message("Setting also BOS output")
else()
	include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )
	set(LINK_DIRECTORIES
	  ${ROOT_LIBRARY_DIR}	
	  )
	file(GLOB sources 
	  ${PROJECT_SOURCE_DIR}/EdGen.cxx
	  ${PROJECT_SOURCE_DIR}/EdInput.cxx
	  ${PROJECT_SOURCE_DIR}/EdModel.cxx
	  ${PROJECT_SOURCE_DIR}/EdHisto.cxx
	  ${PROJECT_SOURCE_DIR}/EdOutput.cxx
	  ${PROJECT_SOURCE_DIR}/EdPhysics.cxx
	  ${PROJECT_SOURCE_DIR}/EdGenPhaseSpace.cxx
	  ${PROJECT_SOURCE_DIR}/EdProcess.cxx
	  ${PROJECT_SOURCE_DIR}/EdCrossSection.cxx
	  )
	file(GLOB headers 
	  ${PROJECT_SOURCE_DIR}/EdGen.h
	  ${PROJECT_SOURCE_DIR}/EdInput.h
	  ${PROJECT_SOURCE_DIR}/EdModel.h
	  ${PROJECT_SOURCE_DIR}/EdHisto.h
	  ${PROJECT_SOURCE_DIR}/EdOutput.h
	  ${PROJECT_SOURCE_DIR}/EdPhysics.h	
	  ${PROJECT_SOURCE_DIR}/EdGenPhaseSpace.h	
	  ${PROJECT_SOURCE_DIR}/EdProcess.h
	  ${PROJECT_SOURCE_DIR}/EdCrossSection.h
	  )
	message("No CLAS6LIB found")
endif()

link_directories( ${LINK_DIRECTORIES})


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(EdGen ${sources} ${headers})
if(WITHBOS)
	set(LINK_LIBRARIES
  ${ROOT_LIBRARIES} -lThread -lHtml -lGed -lEG -lGeom -lc_bos_io -lbankdefs -lbosio -lpacklib
)
else()
	set(LINK_LIBRARIES
  ${ROOT_LIBRARIES} -lThread -lHtml -lGed -lEG -lGeom 
)
endif()
target_link_libraries(EdGen ${LINK_LIBRARIES})


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/output) 
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS EdGen DESTINATION bin)


