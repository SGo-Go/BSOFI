# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Script with general rules for makefiles.
# ------------------------------------------------------------------

OBJ     ?= $(patsubst %.f90, %.o, \
	$(patsubst %.f, %.o, \
	$(patsubst %.c, %.o, $(SOURCES:.cpp=.o))))
LIBS    += $(if $(filter %.f, $(SOURCES)), $(FLIBS),)
LINKER  ?= $(CXX)
.PHONY  : clean install all #$(TARGET) 
ifneq   "$(TARGET)" ""
all: $(TARGET)
else
all:;@echo 'Nothing to make!'
endif

$(TARGET): $(OBJ)
	$(LINKER) $(^) $(LIBS) -o $@

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp $(HEADERSXX)
	$(CXX) $(CXXFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

#$(OBJ:.o=$(THIS_DIR)/.o)
# $(addprefix $(THIS_DIR)/, $(OBJ))

clean:
	$(RM) $(OBJ) 
	$(RM) $(addprefix $(THIS_DIR)/, *~ $(TARGET))

install: $(TARGET)
	$(MKDIR) $(INSTALL_DIR)
	$(CP) $(TARGET) $(INSTALL_DIR)/$(TARGET)
