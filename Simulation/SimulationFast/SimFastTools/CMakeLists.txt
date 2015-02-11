################################################################################
# Package: SimFastTools
################################################################################
gaudi_subdir(SimFastTools v1r0)

gaudi_depends_on_subdirs(GaudiAlg
                         GaudiKernel
			 GaudiSvc
			 Simulation/SimulationInterfaces
			 RecoGeometry)

find_package(ROOT COMPONENTS MathCore GenVector Geom)

find_package(DD4hep)

gaudi_install_headers(SimFastTools)

gaudi_add_module(SimFastTools
                 src/*.cxx
                 INCLUDE_DIRS ROOT GaudiKernel DD4hep SimulationInterfaces RecoGeometry 
                 LINK_LIBRARIES DataObjects ROOT GaudiKernel DD4hep RecoGeometry) 

