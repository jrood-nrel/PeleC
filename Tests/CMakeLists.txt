#=============================================================================
# PeleC Testing
#=============================================================================

if(PELEC_ENABLE_FCOMPARE_FOR_TESTS)
  if("${PELEC_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To reference gold files, PELEC_REFERENCE_GOLDS_DIRECTORY must be set and exist")
  else()
    if(EXISTS ${PELEC_REFERENCE_GOLDS_DIRECTORY})
      set(GOLD_FILES_DIRECTORY ${PELEC_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for reference gold files does not exist: ${PELEC_REFERENCE_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

if(PELEC_SAVE_GOLDS)
  if("${PELEC_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To save gold files, PELEC_SAVED_GOLDS_DIRECTORY must be set and exist")
  else()
    if(EXISTS ${PELEC_SAVED_GOLDS_DIRECTORY})
      set(SAVED_GOLDS_DIRECTORY ${PELEC_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${PELEC_SAVED_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

macro(setup_test)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/RegTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/PeleC-${TEST_EXE_DIR})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_EXE_DIR}/tests/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Find fcompare
    if(PELEC_ENABLE_FCOMPARE_FOR_TESTS)
      set(FCOMPARE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/fcompare)
    endif()
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
    # Copy files to test working directory
    file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.inp DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests
    set(RUNTIME_OPTIONS "amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amrex.the_arena_is_managed=0")
    if(PELEC_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(PELEC_ENABLE_MPI)
      if(PELEC_ENABLE_CUDA)
        set(PELEC_NP 2) # 1 rank per GPU on Eagle
      else()
        set(PELEC_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELEC_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(PELEC_NP 1)
      unset(MPI_COMMANDS)
    endif()
    # Use fcompare to test diffs in plots against gold files
    if(PELEC_ENABLE_FCOMPARE_FOR_TESTS AND (NOT "${TEST_NAME}" MATCHES "hdf5$"))
      if(PELEC_ENABLE_CUDA)
        set(FCOMPARE_TOLERANCE "-r 1e-12 --abs_tol 1.0e-12")
      endif()
      set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE} ${FCOMPARE_TOLERANCE} ${PLOT_TEST} ${PLOT_GOLD}")
    endif()
    if(PELEC_SAVE_GOLDS)
      file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME})
      set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME}/")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUNTIME_OPTIONS "max_step=10 ${RUNTIME_OPTIONS}")
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELEC_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME TEST_EXE_DIR)
    add_test_r(${TEST_NAME} ${TEST_EXE_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no-ci")
endfunction(add_test_re)

#=============================================================================
# Regression tests
#=============================================================================

# Run in CI
add_test_r(sod-1 Sod)
add_test_r(sod-2 Sod)

# Not run in CI
add_test_re(sedov-1 Sedov)
add_test_re(shu-osher-1 Shu-Osher)
if(PELEC_DIM GREATER 2)
  add_test_re(eb-bluffbody-1 EB-BluffBody)
endif()
