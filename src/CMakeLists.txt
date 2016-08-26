add_executable (hydro_analysis.e
  main.cpp 
  Stopwatch.h
  FluidcellStatistic.cpp
  Hydroinfo_h5.cpp
  Hydroinfo_MUSIC.cpp
  ParameterReader.cpp
  arsenal.cpp
  cornelius.cpp
  SurfaceFinder.cpp
  )
target_link_libraries (hydro_analysis.e ${LIBS})

install(TARGETS hydro_analysis.e DESTINATION ${CMAKE_HOME_DIRECTORY})
