######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += set_masses.f90 nuclear_verlet.f90
INCLUDES += write_geom.f90 write_energy.f90

$(obj_path)/liosubs.o : $(INCLUDES) liosubs.mk
######################################################################