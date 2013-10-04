
################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS OptVesselVisModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/vesselvisraycaster.cpp
	${MOD_DIR}/processors/vesselvissliceviewer.cpp
    ${MOD_DIR}/properties/spatialtransfuncproperty.cpp
	${MOD_DIR}/properties/link/linkevaluatorspatialtransfuncid.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/vesselvisraycaster.h
	${MOD_DIR}/processors/vesselvissliceviewer.h
    ${MOD_DIR}/properties/spatialtransfuncproperty.h
	${MOD_DIR}/properties/link/linkevaluatorspatialtransfuncid.h
)
  
# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS OptVesselVisModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/spatialtransfunc1dkeyseditor.cpp
	${MOD_DIR}/qt/spatialtransfunceditor.cpp
	${MOD_DIR}/qt/spatialtransfuncmappingcanvas.cpp
	${MOD_DIR}/qt/spatialtransfuncplugin.cpp
	${MOD_DIR}/qt/spatialtransfuncpropertywidget.cpp
	${MOD_DIR}/qt/spatialtransfuncwidgetfactory.cpp
)  
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/spatialtransfunc1dkeyseditor.h
	${MOD_DIR}/qt/spatialtransfunceditor.h
	${MOD_DIR}/qt/spatialtransfuncmappingcanvas.h
	${MOD_DIR}/qt/spatialtransfuncplugin.h
	${MOD_DIR}/qt/spatialtransfuncpropertywidget.h
	${MOD_DIR}/qt/spatialtransfuncwidgetfactory.h
)

SET(MOD_QT_HEADERS_NONMOC
	${MOD_DIR}/qt/spatialtransfuncwidgetfactory.h
)
