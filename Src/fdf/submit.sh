PROCS=32
INITDIR=$(PWD)
WORKDIR=$(INITDIR)/work$(PROCS)
CUES.SH=$(WORKDIR)/cues.sh
FILES=sample sample.fdf coords.fdf XY.fdf Coords.dat Otherfile

run:
	@if [ -d $(WORKDIR) ]; then   \
          rm -rf $(WORKDIR);          \
        fi
	@mkdir $(WORKDIR)
	@cp $(FILES) $(WORKDIR)
	@echo "#! /bin/bash"                       > $(CUES.SH)
	@echo "# @ class = bsc_cs"                >> $(CUES.SH)
	@echo "# @ initialdir = $(WORKDIR)"       >> $(CUES.SH)
	@echo "# @ output = ./cues.$(PROCS).out"  >> $(CUES.SH)
	@echo "# @ error  = ./cues.$(PROCS).err"  >> $(CUES.SH)
	@echo "# @ total_tasks = $(PROCS)"        >> $(CUES.SH)
	@echo "# @ wall_clock_limit = 01:00:00"   >> $(CUES.SH)
	@echo "time srun ./sample"                >> $(CUES.SH)
	@echo " "
	@echo "########################################"
	@echo "Envian a cues l'script cues.sh"
	@echo "########################################"
	@cat $(CUES.SH)
	@mnsubmit $(CUES.SH)

clean:
	@rm -rf work*
