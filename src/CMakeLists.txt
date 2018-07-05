set(coco_header_files
  ${coco_header_files}
	${CMAKE_CURRENT_SOURCE_DIR}/ToolInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/options.h
  ${CMAKE_CURRENT_SOURCE_DIR}/util.h
  ${CMAKE_CURRENT_SOURCE_DIR}/types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Lookuptable.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mathsupport.h
	PARENT_SCOPE
	)

set(coco_source_files
  ${coco_source_files}
	${CMAKE_CURRENT_SOURCE_DIR}/pcoverage.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pcreads.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/options.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Lookuptable.cpp
  PARENT_SCOPE
        )

