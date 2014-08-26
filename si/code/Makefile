SHELL   = /bin/sh

# define detail root directory
EXEDIR = ..
BASEDIR = $(EXEDIR)/code
MAKEDIR = $(BASEDIR)
CODEDIR = $(BASEDIR)
EXEFILE = si_lineform

PREOBJCS= 

FF      = ifort
FFLAGS  = -vms -assume nounderscore -align dcommons -cpp -c  -g0 -O3 -shared-intel # -mcmodel=medium #nounderflow

CC      = cc
CFLAGS  = -c

LD      = ifort -o

F_FILES = $(wildcard *.f)
C_FILES = $(wildcard *.c)
FOBJCS  = $(patsubst %.f,%.o,$(F_FILES)) 
COBJCS  = $(patsubst %.c,%.o,$(C_FILES)) 

$(EXEDIR)/$(EXEFILE): $(FOBJCS) $(COBJCS) $(PREOBJCS)
	$(LD) $@ $(COBJCS) $(FOBJCS) \
	      $(PREOBJCS) 

$(FOBJCS): %.o: %.f
	$(FF) $(FFLAGS) $(CODEDIR)/$<

$(COBJCS): %.o: %.c
	$(CC)  $(CFLAGS) $(CODEDIR)/$<


lf_decl0.inc 	: lf_param.inc atm_def.inc
mol_data.inc    : lf_decl0.inc lf_param.inc 
elem_data.inc   : lf_decl0.inc lf_param.inc
iondat.inc 	: lf_param.inc 
lineform.o 	: lineform.f physcnst.inc lf_param.inc iondat.inc atm_def.inc spline.f spline_integral.f ew.f \
                  partpress.f \
	          lf_decl0.inc lf_decl1.inc elemdata.inc mol_data.inc \
                  lf_info.inc spol.f hyd_esw.f
kappa.o		: kappa.f lf_param.inc atm_def.inc mol_data.inc lf_info.inc
atmosin1.o      : atmosin1.f physcnst.inc atm_def.inc atmgrid.inc lf_param.inc \
                  lf_decl0.inc
partpress.o     : partpress.f physcnst.inc lf_decl0.inc iondat.inc ppcommon.inc
read_lines.o    : read_lines.f physcnst.inc lf_param.inc atm_def.inc \
		  lf_decl0.inc lf_decl1.inc lf_decl2.inc 
hydline.o       : hydline.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hyd_esw.o       : hyd_esw.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hydl_ag.o       : hydl_ag.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hydl_ag_con.o   : hydl_ag_con.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hydl_ag_con_fs.o : hydl_ag_con_fs.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hydl_bpo_con.o  : hydl_bpo_con.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hydl_bpo_con_fs.o : hydl_bpo_con_fs.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
hydl_ag.o       : hydl_ag.f  physcnst.inc lf_param.inc atm_def.inc lf_decl0.inc
convol.o	: convol.f   physcnst.inc
iterate_source_fct.o: iterate_source_fct.f  lf_decl3.inc deriv.o  deriv1.f
deriv.o         :         deriv1.f
read_vcs_hyd.o	: read_vcs_hyd.f  physcnst.inc lf_param.inc atm_def.inc \
		  lf_decl0.inc
read_hydprof.o	: read_hydprof.f  physcnst.inc lf_param.inc atm_def.inc \
		  lf_decl0.inc
resample.o	: resample.f physcnst.inc
resample_ref.o	: resample_ref.f physcnst.inc lf_ref.inc
reson_radia.o   : reson_radia.f  physcnst.inc lf_param.inc atm_def.inc \
                  lf_decl0.inc
check_stepwidth_crit.o : check_stepwidth_crit.f 
exec_lineform.o : exec_lineform.f

clean:
	rm -f *.o *~

cleanall:
	$(MAKE) clean
	rm -f $(TARGET)

.PHONY: all clean cleanall
