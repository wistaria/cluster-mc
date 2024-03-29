enable_testing()
include(CTest)

macro(add_iotest)
  if(${ARGC} EQUAL 1)
    set(name ${ARGV0})
    set(input ${ARGV0})
    set(output ${ARGV0})
  else(${ARGC} EQUAL 1)
    if(${ARGC} EQUAL 2)
      set(name ${ARGV0})
      set(input ${ARGV1})
      set(output ${ARGV1})
    else(${ARGC} EQUAL 2)
      set(name ${ARGV0})
      set(input ${ARGV1})
      set(output ${ARGV2})
    endif(${ARGC} EQUAL 2)
  endif(${ARGC} EQUAL 1)
  enable_testing()
  add_test(${name}
    ${CMAKE_COMMAND}
      -Dcmd=${name}
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=${input}
      -Doutput=${output}
      -P ${PROJECT_SOURCE_DIR}/cmake/run_iotest.cmake
    )
endmacro(add_iotest)
