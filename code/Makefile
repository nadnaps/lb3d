# Makefile created by mkmf.pl $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include /data/ernie4/schmie/src/lb3d.refac/lbe/code/Makefile.template


.DEFAULT:
	-touch $@
all: getrev lb3d
lb3d.o: ./lb3d.F90 lb3d.h lb3d_initialise.o lb3d_time_loop.o lb3d_finalise.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d.F90
lb3d_advection.o: ./lb3d_advection.F90 lb3d.h lb3d_lattice.o lb3d_bc.o lb3d_timer.o lb3d_analysis.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_advection.F90
lb3d_analysis.o: ./lb3d_analysis.F90 lb3d.h lb3d_config.o lb3d_lattice.o lb3d_mpi.o lb3d_global.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_analysis.F90
lb3d_bc.o: ./lb3d_bc.F90 lb3d.h lb3d_config.o lb3d_timer.o lb3d_log.o lb3d_lattice.o lb3d_analysis.o lb3d_invasion.o lb3d_force.o lb3d_bc_periodic.o lb3d_bc_leesedwards.o lb3d_bc_bounce.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_bc.F90
lb3d_bc_bounce.o: ./lb3d_bc_bounce.F90 lb3d.h lb3d_config.o lb3d_lattice.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_bc_bounce.F90
lb3d_bc_leesedwards.o: ./lb3d_bc_leesedwards.F90 lb3d.h lb3d_bdist.o lb3d_mpi_parameters.o lb3d_config.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_bc_leesedwards.F90
lb3d_bc_periodic.o: ./lb3d_bc_periodic.F90 lb3d.h lb3d_lattice.o lb3d_mpi_parameters.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_bc_periodic.F90
lb3d_bdist.o: ./lb3d_bdist.F90 lb3d.h lb3d_mpi.o lb3d_log.o lb3d_config.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_bdist.F90
lb3d_collision.o: ./lb3d_collision.F90 lb3d.h lb3d_global.o lb3d_config.o lb3d_bc_leesedwards.o lb3d_bc_periodic.o lb3d_timer.o lb3d_force_shanchen.o lb3d_lattice.o lb3d_bdist.o lb3d_analysis.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_collision.F90
lb3d_compatibility.o: ./lb3d_compatibility.F90 lb3d_config.o lb3d_mpi.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./lb3d_compatibility.F90
lb3d_config.o: ./lb3d_config.F90 lb3d.h lb3d_revision.h lb3d_log.o lb3d_global.o lb3d_mpi_parameters.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_config.F90
lb3d_finalise.o: ./lb3d_finalise.F90 lb3d.h lb3d_config.o lb3d_bc_leesedwards.o lb3d_timer.o lb3d_io_hdf5.o lb3d_mpi.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_finalise.F90
lb3d_force.o: ./lb3d_force.F90 lb3d.h lb3d_global.o lb3d_config.o lb3d_log.o lb3d_lattice.o lb3d_mpi.o lb3d_force_constant.o lb3d_force_kolmogorov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_force.F90
lb3d_force_constant.o: ./lb3d_force_constant.F90 lb3d.h lb3d_global.o lb3d_io_helper.o lb3d_analysis.o lb3d_config.o lb3d_lattice.o lb3d_mpi_parameters.o lb3d_mpi.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_force_constant.F90
lb3d_force_kolmogorov.o: ./lb3d_force_kolmogorov.F90 lb3d.h lb3d_global.o lb3d_mpi.o lb3d_log.o lb3d_config.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_force_kolmogorov.F90
lb3d_force_shanchen.o: ./lb3d_force_shanchen.F90 lb3d.h lb3d_global.o lb3d_lattice.o lb3d_config.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_force_shanchen.F90
lb3d_global.o: ./lb3d_global.F90 lb3d.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_global.F90
lb3d_helper.o: ./lb3d_helper.F90 lb3d.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_helper.F90
lb3d_ic.o: ./lb3d_ic.F90 lb3d.h lb3d_config.o lb3d_global.o lb3d_lattice.o lb3d_bdist.o lb3d_log.o lb3d_ic_rock.o lb3d_io_checkpoint.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_ic.F90
lb3d_ic_rock.o: ./lb3d_ic_rock.F90 lb3d.h lb3d_lattice.o lb3d_io_hdf5.o lb3d_io_xdrf.o lb3d_mpi.o lb3d_config.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_ic_rock.F90
lb3d_initialise.o: ./lb3d_initialise.F90 lb3d.h lb3d_io.o lb3d_config.o lb3d_mpi.o lb3d_random.o lb3d_lattice.o lb3d_bc.o lb3d_ic.o lb3d_force.o lb3d_log.o lb3d_timer.o lb3d_compatibility.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_initialise.F90
lb3d_invasion.o: ./lb3d_invasion.F90 lb3d.h lb3d_global.o lb3d_mpi.o lb3d_config.o lb3d_bdist.o lb3d_force_shanchen.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_invasion.F90
lb3d_io.o: ./lb3d_io.F90 lb3d.h lb3d_global.o lb3d_lattice.o lb3d_io_helper.o lb3d_io_arrstats.o lb3d_io_stress.o lb3d_io_hdf5.o lb3d_io_checkpoint.o lb3d_io_dump_data.o lb3d_config.o lb3d_mpi.o lb3d_mpi_parameters.o lb3d_timer.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io.F90
lb3d_io_arrstats.o: ./lb3d_io_arrstats.F90 lb3d.h lb3d_global.o lb3d_lattice.o lb3d_config.o lb3d_log.o lb3d_mpi.o lb3d_analysis.o lb3d_io_helper.o lb3d_io_dump_data.o lb3d_io_hdf5.o lb3d_io_xdrf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_arrstats.F90
lb3d_io_checkpoint.o: ./lb3d_io_checkpoint.F90 lb3d.h lb3d_global.o lb3d_io_arrstats.o lb3d_config.o lb3d_helper.o lb3d_io_helper.o lb3d_bc_leesedwards.o lb3d_io_xdrf.o lb3d_io_hdf5.o lb3d_lattice.o lb3d_mpi_parameters.o lb3d_mpi.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_checkpoint.F90
lb3d_io_dump_data.o: ./lb3d_io_dump_data.F90 lb3d.h lb3d_io_hdf5.o lb3d_io_xdrf.o lb3d_global.o lb3d_config.o lb3d_io_helper.o lb3d_force_shanchen.o lb3d_analysis.o lb3d_lattice.o lb3d_mpi.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_dump_data.F90
lb3d_io_hdf5.o: ./lb3d_io_hdf5.F90 lb3d.h lb3d_mpi.o lb3d_config.o lb3d_io_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_hdf5.F90
lb3d_io_helper.o: ./lb3d_io_helper.F90 lb3d.h lb3d_mpi.o lb3d_config.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_helper.F90
lb3d_io_stress.o: ./lb3d_io_stress.F90 lb3d.h lb3d_config.o lb3d_global.o lb3d_mpi.o lb3d_helper.o lb3d_io_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_stress.F90
lb3d_io_xdrf.o: ./lb3d_io_xdrf.F90 lb3d.h lb3d_mpi.o lb3d_log.o lb3d_lattice.o lb3d_config.o lb3d_io_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_io_xdrf.F90
lb3d_lattice.o: ./lb3d_lattice.F90 lb3d.h lb3d_global.o lb3d_config.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_lattice.F90
lb3d_log.o: ./lb3d_log.F90 lb3d.h lb3d_mpi_parameters.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_log.F90
lb3d_mpi.o: ./lb3d_mpi.F90 lb3d.h lb3d_config.o lb3d_global.o lb3d_lattice.o lb3d_mpi_parameters.o lb3d_log.o lb3d_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_mpi.F90
lb3d_mpi_parameters.o: ./lb3d_mpi_parameters.F90 lb3d.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_mpi_parameters.F90
lb3d_random.o: ./lb3d_random.F90 lb3d.h lb3d_global.o lb3d_log.o lb3d_config.o lb3d_mpi_parameters.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_random.F90
lb3d_sanity.o: ./lb3d_sanity.F90 lb3d.h lb3d_lattice.o lb3d_config.o lb3d_log.o lb3d_mpi.o lb3d_analysis.o lb3d_io.o lb3d_io_helper.o lb3d_global.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_sanity.F90
lb3d_time_loop.o: ./lb3d_time_loop.F90 lb3d.h lb3d_global.o lb3d_io.o lb3d_timer.o lb3d_force.o lb3d_bc.o lb3d_advection.o lb3d_collision.o lb3d_helper.o lb3d_sanity.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lb3d_time_loop.F90
lb3d_timer.o: ./lb3d_timer.F90 lb3d_global.o lb3d_mpi.o lb3d_mpi_parameters.o lb3d_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./lb3d_timer.F90
SRC = ./lb3d_config.F90 ./lb3d_io_arrstats.F90 ./lb3d_advection.F90 ./lb3d_compatibility.F90 ./lb3d_bc.F90 ./lb3d_ic_rock.F90 ./lb3d_sanity.F90 ./lb3d_time_loop.F90 ./lb3d.F90 ./lb3d_collision.F90 ./lb3d_initialise.F90 ./lb3d_finalise.F90 ./lb3d_force_constant.F90 ./lb3d_mpi.F90 ./lb3d_ic.F90 ./lb3d_timer.F90 ./lb3d_io_checkpoint.F90 ./lb3d_io_xdrf.F90 ./lb3d_io_dump_data.F90 ./lb3d_bdist.F90 ./lb3d_io_hdf5.F90 ./lb3d_bc_leesedwards.F90 ./lb3d_lattice.F90 ./lb3d_mpi_parameters.F90 ./lb3d_io_helper.F90 ./lb3d_analysis.F90 ./lb3d_force.F90 ./lb3d_global.F90 ./lb3d_bc_periodic.F90 ./lb3d_io.F90 ./lb3d_io_stress.F90 ./lb3d_force_shanchen.F90 ./lb3d_log.F90 ./lb3d_random.F90 ./lb3d_helper.F90 ./lb3d_bc_bounce.F90 ./lb3d_force_kolmogorov.F90 ./lb3d_invasion.F90 lb3d.h lb3d_revision.h
OBJ = lb3d_config.o lb3d_io_arrstats.o lb3d_advection.o lb3d_compatibility.o lb3d_bc.o lb3d_ic_rock.o lb3d_sanity.o lb3d_time_loop.o lb3d.o lb3d_collision.o lb3d_initialise.o lb3d_finalise.o lb3d_force_constant.o lb3d_mpi.o lb3d_ic.o lb3d_timer.o lb3d_io_checkpoint.o lb3d_io_xdrf.o lb3d_io_dump_data.o lb3d_bdist.o lb3d_io_hdf5.o lb3d_bc_leesedwards.o lb3d_lattice.o lb3d_mpi_parameters.o lb3d_io_helper.o lb3d_analysis.o lb3d_force.o lb3d_global.o lb3d_bc_periodic.o lb3d_io.o lb3d_io_stress.o lb3d_force_shanchen.o lb3d_log.o lb3d_random.o lb3d_helper.o lb3d_bc_bounce.o lb3d_force_kolmogorov.o lb3d_invasion.o
clean: neat
	-rm -f .cppdefs $(OBJ) lb3d
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lb3d: $(OBJ) 
	$(LD) $(OBJ) -o lb3d  $(LDFLAGS)
getrev:
	echo "#define LBEREVISION \"unknown\"" | sed "s/unknown/`svnversion`/" | sed "s/\"\"/\"unknown\"/" > lb3d_revision.new
	diff lb3d_revision.new lb3d_revision.h > /dev/null 2>&1 ; if [ $$? == 0 ] ; then rm lb3d_revision.new ; else echo "Creating new lb3d_revision.h"; mv lb3d_revision.new lb3d_revision.h ; fi
	echo "#define LBEFLAGS \"$(DEFFLAGS) $(MAKEFFLAGS)\"" >> lb3d_revision.h
	cat lb3d_revision.h
