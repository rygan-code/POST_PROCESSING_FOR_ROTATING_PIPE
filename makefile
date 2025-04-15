f77=mpif90

DIR = OBJ
opt= -O3 -J $(DIR)/ -I $(DIR)/ 
EXEC = post.out
srs= parameters.F90 allocate.F90 helmholtz_decomposition.F90 viscous.F90 integral_transform.F90 sgs_transport.F90 postproc3d-1.1_mpi.F90	 
OBJS=$(addprefix $(DIR)/,$(patsubst %.F90,%.o,$(srs)))

default: $(EXEC)

$(EXEC): $(OBJS)
	$(f77) -O3  -o $(EXEC) $(OBJS)

$(DIR)/%.o:%.F90
	@if [ ! -d $(DIR) ]; then mkdir -p $(DIR); fi
	$(f77) $(opt) -c $< -o $@

clean:
	rm -rf $(DIR) *.out
