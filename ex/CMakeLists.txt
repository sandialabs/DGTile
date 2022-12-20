set(HYDRO_SOURCES
  amr.cpp
  ics.cpp
  io.cpp
  hydro.cpp
  limit.cpp
  physics.cpp
)

set_source_files_properties(${HYDRO_SOURCES} PROPERTIES LANGUAGE ${p3a_LANGUAGE})
add_executable(hydro ${HYDRO_SOURCES})
set_target_properties(hydro PROPERTIES ${p3a_LANGUAGE}_ARCHITECTURES "${p3a_ARCHITECTURES}")
target_link_libraries(hydro PRIVATE dgtile)
install(TARGETS hydro)

function(config_advect dim)
  foreach(n e IN ZIP_LISTS Ns Es)
    set(type "unif")
    if(smr)
      set(type "smr")
    endif()
    set(name "advect_${dim}D_${type}_p${p}_n${n}")
    configure_file(advect_${dim}D.in ${name} @ONLY)
  endforeach()
endfunction()

function(config_isentropic_vortex)
  foreach(n e IN ZIP_LISTS Ns Es)
    set(name "isentropic_vortex_p${p}_n${n}")
    configure_file(isentropic_vortex.in ${name} @ONLY)
  endforeach()
endfunction()

# -----------------------
# -- 2D advect uniform --
# -----------------------

set(smr "")
set(nb 1)

set(p 0)
set(Ns 4 8 16 32 64)
set(Es 1.7025879373213904e-01 1.3634291271013158e-01 9.1113573964130959e-02 5.3492442263265577e-02 2.9124864231421629e-02)
config_advect(2)

set(p 1)
set(Ns 2 4 8 16 32)
set(Es 1.6433757420765480e-01 5.3199121238565260e-02 1.1349211798098259e-02 2.5519840382122400e-03 6.1380706352362512e-04) 
config_advect(2)

set(p 2)
set(Ns 2 4 8 16 32)
set(Es 3.9424495484157947e-02 4.6976389453536353e-03 5.9235914723061549e-04 7.3916107659276610e-05 9.2375200856581065e-06)
config_advect(2)

# -----------------------
# -- 3D advect uniform --
# -----------------------

set(smr "")
set(nb 1)

set(p 0)
set(Ns 4 8 16 32 64)
set(Es 1.7229122207284575e-01 1.4172757581753875e-01 9.6557751447742729e-02 5.7309264579106817e-02 3.1430030518590252e-02)
config_advect(3)

set(p 1)
set(Ns 2 4 8 16 32)
set(Es 1.6243116280891215e-01 5.4866865890234651e-02 1.2708718332308662e-02 3.0015833717076455e-03 7.3432632162546746e-04)
config_advect(3)

set(p 2)
set(Ns 2 4 8 16 32)
set(Es 4.3937922550760400e-02 5.8754690475273843e-03 7.2221925969691736e-04 9.0483435889345134e-05 1.1310704530546172e-05)
config_advect(3)

# -------------------
# -- 2D advect smr --
# -------------------

set(smr "first_block")
set(nb 2)

set(p 0)
set(Ns 8 16 32 64)
set(Es 9.9681826032634493e-02 6.0238198239686611e-02 3.3389085369575383e-02 1.7625438862315355e-02)
config_advect(2)

set(p 1)
set(Ns 4 8 16)
set(Es 9.6376001289121566e-03 2.1767365023989704e-03 5.2299734065147728e-04)
config_advect(2)

set(p 2)
set(Ns 2 4 8)
set(Es 4.2816587063393367e-03 5.2425590235762346e-04 6.5156118404552531e-05)
config_advect(2)

# ------------------
# -- 3D advet smr --
# ------------------

set(smr "first_block")
set(nb 2)

set(p 0)
set(Ns 8 16 32 64)
set(Es 1.0826460310762859e-01 6.6617795524005846e-02 3.7256017615535764e-02 1.9749500590837927e-02)
config_advect(3)

set(p 1)
set(Ns 4 8 16)
set(Es 1.1773097976716489e-02 2.7878105898248256e-03 6.8188429681937822e-04)
config_advect(3)

set(p 2)
set(Ns 2 4 8)
set(Es 5.5964279214205086e-03 6.8741058484760508e-04 8.5697832571347317e-05)
config_advect(3)

# -------------------
# -- 2D advect amr --
# -------------------

set(p 0)
set(e 1.4647887464513659e-01)
configure_file(advect_2D_amr.in advect_2D_amr_p${p} @ONLY)

set(p 1)
set(e 9.8918756395325342e-03)
configure_file(advect_2D_amr.in advect_2D_amr_p${p} @ONLY)

set(p 2)
set(e 5.8992188720198211e-04)
configure_file(advect_2D_amr.in advect_2D_amr_p${p} @ONLY)

# ---------------
# -- 2D RT amr --
# ---------------

set(p 1)
configure_file(rt.in rt @ONLY)

# ----------------
# -- 2D sod amr --
# ----------------

set(p 1)
configure_file(sod.in sod @ONLY)

# --
# -- Regression tests
# --

set(base ${CMAKE_CURRENT_BINARY_DIR})

set(SHORT_TESTS
  advect_2D_unif_p0_n4
  advect_2D_unif_p0_n8
  advect_2D_unif_p0_n16
  advect_2D_unif_p1_n2
  advect_2D_unif_p1_n4
  advect_2D_unif_p1_n8
  advect_2D_unif_p2_n2
  advect_2D_unif_p2_n4
  advect_2D_smr_p0_n8
  advect_2D_smr_p1_n4
  advect_2D_smr_p2_n2
  advect_2D_amr_p0
  advect_2D_amr_p1
  advect_2D_amr_p2
  advect_3D_unif_p0_n4
  advect_3D_unif_p0_n8
  advect_3D_unif_p0_n16
  advect_3D_unif_p1_n2
  advect_3D_unif_p1_n4
  advect_3D_unif_p1_n8
  advect_3D_unif_p2_n2
  advect_3D_unif_p2_n4
  advect_3D_smr_p0_n8
  advect_3D_smr_p1_n4
  advect_3D_smr_p2_n2)

foreach(name IN LISTS SHORT_TESTS)
  add_test(NAME ${name} COMMAND "hydro" "${base}/${name}")
  set_tests_properties(${name} PROPERTIES LABELS "short")
endforeach()

set(SHORT_MPI_TESTS
  advect_2D_smr_p0_n8
  advect_2D_smr_p1_n4
  advect_2D_smr_p2_n2
  advect_2D_amr_p0
  advect_2D_amr_p1
  advect_2D_amr_p2)

foreach(name IN LISTS SHORT_MPI_TESTS)
  add_test(NAME ${name}_mpi COMMAND "mpirun" "-n" "4" "hydro" "${base}/${name}")
  set_tests_properties(${name}_mpi PROPERTIES LABELS "short;mpi")
endforeach()
