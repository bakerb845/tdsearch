CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/tdsearch_config.h)
SET(HEADERS tdsearch_commands.h tdsearch_data.h tdsearch_event.h
            tdsearch_greens.h tdsearch_gridsearch.h tdsearch.h
            tdsearch_hudson.h tdsearch_struct.h tdsearch_config.h)
FILE(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
INSTALL(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)
