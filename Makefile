EXEC  = L-Galaxies

# Default object files (others may be added with -D options)
OBJS  = /Users/Bam/ClionProjects/Rings/main.o \
	/Users/Bam/ClionProjects/Rings/io_tree.o \
	/Users/Bam/ClionProjects/Rings/init.o \
	/Users/Bam/ClionProjects/Rings/cool_func.o \
	/Users/Bam/ClionProjects/Rings/save.o \
	/Users/Bam/ClionProjects/Rings/save_galtree.o \
	/Users/Bam/ClionProjects/Rings/mymalloc.o \
	/Users/Bam/ClionProjects/Rings/read_parameters.o \
	/Users/Bam/ClionProjects/Rings/peano.o \
	/Users/Bam/ClionProjects/Rings/allvars.o \
	/Users/Bam/ClionProjects/Rings/age.o \
	/Users/Bam/ClionProjects/Rings/update_type_two.o \
	/Users/Bam/ClionProjects/Rings/metals.o \
	/Users/Bam/ClionProjects/Rings/model_infall.o \
	/Users/Bam/ClionProjects/Rings/model_h2fraction.o \
	/Users/Bam/ClionProjects/Rings/model_cooling.o \
	/Users/Bam/ClionProjects/Rings/model_starformation_and_feedback.o \
	/Users/Bam/ClionProjects/Rings/model_reincorporation.o \
	/Users/Bam/ClionProjects/Rings/model_mergers.o \
	/Users/Bam/ClionProjects/Rings/model_dust.o \
	/Users/Bam/ClionProjects/Rings/model_misc.o \
	/Users/Bam/ClionProjects/Rings/model_disrupt.o \
	/Users/Bam/ClionProjects/Rings/model_stripping.o \
	/Users/Bam/ClionProjects/Rings/scale_cosmology.o

# The following is used only to set dependencies
INCL  = Makefile \
    /Users/Bam/ClionProjects/Rings/allvars.h \
	/Users/Bam/ClionProjects/Rings/h_funcs.h \
	/Users/Bam/ClionProjects/Rings/h_params.h \
	/Users/Bam/ClionProjects/Rings/h_metals.h \
	/Users/Bam/ClionProjects/Rings/h_galaxy_output.h \
	/Users/Bam/ClionProjects/Rings/h_galaxy_tree_data.h \
	/Users/Bam/ClionProjects/Rings/h_galaxy.h \
	/Users/Bam/ClionProjects/Rings/h_halo_data.h \
	/Users/Bam/ClionProjects/Rings/h_halo_ids_data.h \
	/Users/Bam/ClionProjects/Rings/h_halo_aux_data.h \
	/Users/Bam/ClionProjects/Rings/h_lightcone.h \
	/Users/Bam/ClionProjects/Rings/h_variables.h \
	/Users/Bam/ClionProjects/Rings/proto.h
ifeq (ALL_SKY_LIGHTCONE,$(findstring ALL_SKY_LIGHTCONE,$(OPT)))
INCL  += /Users/Bam/ClionProjects/Rings/lightcone.h
endif

# Either include the default set of Makefile options, or define your own
# include Makefile_options
include My_Makefile_options/My_Makefile_options_rings

# Choose your system type (copy an entry from Makefile_compilers)
# SYSTYPE = "ETH"
# include Makefile_compilers
include Makefile_compilers

#LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas $(HDF5_LIBS) -lhdf5_serial -lhdf5_serial_hl
LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas $(HDF5_LIBS) -lhdf5 -lhdf5_hl

CFLAGS =   -g $(OPTIONS) $(OPT) -DCOMPILETIMESETTINGS=\""$(OPT)"\" $(OPTIMIZE) $(GSL_INCL) $(HDF5_INCL)

all: metadata $(EXEC)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) Makefile Makefile_compilers My_Makefile_options/My_Makefile_options_rings

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# uses -CC compiler option to save comments, needed for HDF5 output
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC_MD} ${OPT} ${CFLAGS} -E -CC /Users/Bam/ClionProjects/Rings/h_galaxy_output.h -o /Users/Bam/ClionProjects/Rings/h_galaxy_output.i
	${CC_MD} ${OPT} ${CFLAGS} -E -CC /Users/Bam/ClionProjects/Rings/h_metals.h -o /Users/Bam/ClionProjects/Rings/h_metals.i
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_TypeString.awk > ./AuxCode/awk/output/L-Galaxies_Types.txt
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_DDL.awk > ./AuxCode/awk/output/L-Galaxies_DDL.sql
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_BIN_2_DDL.awk >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
else
	awk -f ./AuxCode/awk/extract_SFH_Time.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_Time_2_DDL.awk >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
endif	
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_struct.awk >  ./AuxCode/awk/output/idl/LGalaxy.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_hists.awk > ./AuxCode/awk/output/idl/LGalaxy_plot.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_testfloats.awk > ./AuxCode/awk/output/idl/LGalaxy_testfloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_zerofloats.awk > ./AuxCode/awk/output/idl/LGalaxy_zerofloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_LGalaxy.awk > ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_FileFormat.awk > ./AuxCode/awk/output/L-Galaxies_FileFormat.csv
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/MOMAF_INPUT_2_MoMaFGalaxy.awk >> ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_python_struct.awk >  ./AuxCode/awk/output/python/LGalaxy.py
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_3_HDF5.awk > /Users/Bam/ClionProjects/Rings/io_hdf5.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT_props.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_prop_2_HDF5_proptable.awk > ./input/hdf5_field_props.txt

	awk -f ./AuxCode/awk/extract_struct_metals.awk /Users/Bam/ClionProjects/Rings/h_metals.i > ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_elements.awk /Users/Bam/ClionProjects/Rings/h_metals.i >> ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/Rings/h_galaxy_output.i >> ./AuxCode/awk/output/structs.dat
