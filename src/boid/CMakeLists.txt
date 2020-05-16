add_custom_target(link_boid ALL
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_LIST_DIR}/boid.py ${CMAKE_SOURCE_DIR}/bin/boid.py
	SOURCES boid.py
)

add_dependencies(link_boid pypsim)