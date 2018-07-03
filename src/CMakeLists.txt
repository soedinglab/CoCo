set(coco_header_files
  ${coco_header_files}
	${CMAKE_CURRENT_SOURCE_DIR}/ToolInfo.h
	PARENT_SCOPE
	)

set(coco_source_files
  ${coco_source_files}
	${CMAKE_CURRENT_SOURCE_DIR}/pcoverage.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pcreads.cpp
  PARENT_SCOPE
        )

