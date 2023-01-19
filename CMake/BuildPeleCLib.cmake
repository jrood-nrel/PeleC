function(build_pelec_lib pelec_lib_name)
  if (NOT (TARGET ${pelec_lib_name}))
    add_library(${pelec_lib_name} OBJECT)

    set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
    set(PELE_PHYSICS_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport")
    set(PELE_PHYSICS_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos")
    set(PELE_PHYSICS_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Mechanism/Models/${PELEC_CHEMISTRY_MODEL}")

    if(CLANG_TIDY_EXE)
      set_target_properties(${pelec_lib_name} PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
    endif()

    include(SetPeleCCompileFlags)
    
    target_include_directories(${pelec_lib_name} PUBLIC "${PELE_PHYSICS_SRC_DIR}/Source")

    target_sources(${pelec_lib_name} PRIVATE
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.cpp
                   ${PELE_PHYSICS_TRANSPORT_DIR}/TransportParams.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/TransportTypes.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Constant.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Simple.H)
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_TRANSPORT_DIR})
    if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Constant")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_CONSTANT_TRANSPORT)
    endif()
    if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_SIMPLE_TRANSPORT)
    endif()

    target_sources(${pelec_lib_name} PRIVATE
                   ${PELE_PHYSICS_EOS_DIR}/EOS.cpp
                   ${PELE_PHYSICS_EOS_DIR}/EOS.H
                   ${PELE_PHYSICS_EOS_DIR}/GammaLaw.H
                   ${PELE_PHYSICS_EOS_DIR}/Fuego.H)
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_EOS_DIR})
    if("${PELEC_EOS_MODEL}" STREQUAL "GammaLaw")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_GAMMALAW_EOS)
    endif()
    if("${PELEC_EOS_MODEL}" STREQUAL "Fuego")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_FUEGO_EOS)
    endif()

    target_sources(${pelec_lib_name} PRIVATE
                   ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.cpp
                   ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.H)
    target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${PELE_PHYSICS_MECHANISM_DIR})

    target_sources(${pelec_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorTypes.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.cpp
    )
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Reactions)

    include(AMReXBuildInfo)
    generate_buildinfo(${pelec_lib_name} ${CMAKE_SOURCE_DIR})
    target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)
    
    if(PELEC_ENABLE_MPI)
      target_link_libraries(${pelec_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
    endif()
    
    target_link_libraries(${pelec_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api AMReX::amrex)

  endif()
endfunction()
