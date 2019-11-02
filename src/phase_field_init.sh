#!/bin/bash
# init_phase_field_model.sh

deform=$1    # Deformability (d = epsilon/alpha)
peclet=$2    # Peclet number (Pe = v/(R*Dr))
run=$3       # Trial number
run_dir=$4   # Run directory

if [ "$#" != 4 ]; then
    echo "usage: deform peclet run run_dir"
    exit 1
fi

# Format input params
deform=$(python -c "print '{:.3f}'.format($deform)")
peclet=$(python -c "print '{:.3f}'.format($peclet)")

# A function for generating random numbers
max_seed=1000000
function get_rand(){
    # Generate a 4-byte random integer using urandom
    rand=$(od -vAn -N4 -tu4 < /dev/urandom)
    echo $rand
}

# Set the model parameters
ncells=100            # Number of cells in the system
ncell_x=10            # Number of cells in each row
ncell_y=10            # Number of cells in each column
confine_radius=8.0    # Length between lattice point on triangular lattice
init_radius=7.0       # Radius of the initial droplet
ideal_radius=12.0     # The radius in the volume term of the free energy
cell_Lx=41            # Width of the local cell field
cell_Ly=41            # Height of the local fcell field
phi0=2.0              # phi_0 in the free energy expression
mu=6000.0             # mu in the volume term of the free energy
epsilon=0.1           # epsilon in the repulsion term of the free energy
rotate_diff=0.0001    # Rotational diffusion coefficient
relax_rate=0.1        # Relaxation rate
nsteps=1000000        # Number of timesteps in the main simulation run
nequil=10000          # Number of timesteps for equilibration
delta_t=0.5           # Time interval between each timestep 

# Set the output frequency of each observable dump
dump_cm_freq=1000
dump_bulk_cm_freq=1000
dump_gyr_freq=1000
dump_field_freq=10000
dump_cell_field_freq=10000
dump_index_field_freq=10000
dump_shape_freq=1000
dump_neighbour_freq=1000
dump_energy_freq=1000
dump_overlap_freq=1000
dump_overlap_field_freq=10000
seed=$(get_rand)

# Set alpha based on deformability
alpha=$(python -c "print '{:f}'.format($epsilon/$deform)")
kappa=$(python -c "print '{:f}'.format($alpha*2.0)")

# Set motility based on Peclet number, rotatioal diff, and ideal radius
v=$(python -c "print '{:f}'.format($peclet*$rotate_diff*$ideal_radius)")
v_all=1.0

# Set a triangular lattice
tmp_cm_file="cm_$seed.tmp"
tmp_shape_file="shape_$seed.tmp"
size=$(python cell_position.py $ncell_x $ncell_y $confine_radius $v $seed $tmp_cm_file)
Lx=$(echo $size | awk '{print $3}')
Ly=$(echo $size | awk '{print $6}')
python cell_shape.py $cell_Lx $cell_Ly $phi0 $tmp_shape_file circle $init_radius

# Create run directory and set file names
sim_name="cell_N_${ncells}_d_${deform}_Pe_${peclet}_run_${run}"
run_dir="${run_dir}/${sim_name}/"

if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
fi

# Name of simulation files
cm_file="cm_${sim_name}.in"
shape_file="shape_${sim_name}.in"
params_file="params_${sim_name}.txt"

# Name of output files
dump_cm_file="pos_${sim_name}.dat"
dump_gyr_file="gyr_${sim_name}.dat"
dump_field_file="field_${sim_name}.dat"
dump_cell_field_file="cell-field_${sim_name}"
dump_index_field_file="index-field_${sim_name}.dat"
dump_bulk_cm_file="pos-bulk_${sim_name}.dat"
dump_shape_file="shape_${sim_name}.dat"
dump_neighbour_file="neigh_${sim_name}.dat"
dump_energy_file="energy_${sim_name}.dat"
dump_overlap_file="olap_${sim_name}.dat"
dump_overlap_field_file="olap-field_${sim_name}.dat"

# Copy the template file
params_file=${run_dir}/$params_file
cp params_template.txt $params_file
mv $tmp_cm_file ${run_dir}/$cm_file
mv $tmp_shape_file ${run_dir}/$shape_file

# Replace macros in template with input values
perl -pi -e "s/PHI0/${phi0}/g" $params_file
perl -pi -e "s/ALPHA/${alpha}/g" $params_file
perl -pi -e "s/KAPPA/${kappa}/g" $params_file
perl -pi -e "s/MU/${mu}/g" $params_file
perl -pi -e "s/EPSILON/${epsilon}/g" $params_file
perl -pi -e "s/ROTATE_DIFF/${rotate_diff}/g" $params_file
perl -pi -e "s/MOTILITY/${v_all}/g" $params_file
perl -pi -e "s/RELAX_RATE/${relax_rate}/g" $params_file
perl -pi -e "s/IDEAL_RADIUS/${ideal_radius}/g" $params_file

perl -pi -e "s/CELL_LX/${cell_Lx}/g" $params_file
perl -pi -e "s/CELL_LY/${cell_Ly}/g" $params_file
perl -pi -e "s/LX/${Lx}/g" $params_file
perl -pi -e "s/LY/${Ly}/g" $params_file

perl -pi -e "s/NCELLS/${ncells}/g" $params_file
perl -pi -e "s/NSTEPS/${nsteps}/g" $params_file
perl -pi -e "s/NEQUIL/${nequil}/g" $params_file
perl -pi -e "s/DELTA_T/${delta_t}/g" $params_file
perl -pi -e "s/SEED/${seed}/g" $params_file

perl -pi -e "s/CM_FILE/${cm_file}/g" $params_file
perl -pi -e "s/SHAPE_FILE/${shape_file}/g" $params_file

# Set dumps
function add_dump() {
    params=$1; file=$2;
    if [ "$file" ]; then
	echo "$params $file" >> $params_file 
    fi
}

# Output the centre of mass of each cell
add_dump "dump_cm $dump_cm_freq 0 main" $dump_cm_file

# Output the centre of mass of the entire system (bulk cm)
add_dump "dump_bulk_cm $dump_bulk_cm_freq main" $dump_bulk_cm_file

# Output the gyration tensor compnents of each cell
add_dump "dump_gyr $dump_gyr_freq 0 main" $dump_gyr_file

# Output the total field of the system
add_dump "dump_field $dump_field_freq 0 main" $dump_field_file

# Output the index field of the system
add_dump "dump_index_field $dump_index_field_freq 0 main" $dump_index_field_file

# Output the nearest neighbours of each cell
add_dump "dump_neighbour $dump_neighbour_freq 0 main" $dump_neighbour_file

# Output the energy of the system
#add_dump "dump_energy $dump_energy_freq 0 main" $dump_energy_file

# Output the average local overlap of each cell
#add_dump "dump_overlap $dump_overlap_freq 0 main" $dump_overlap_file

# Output the geometric properties (perimeter and area) of each cell
#add_dump "dump_shape 4 25 4.0 5 31 $dump_shape_freq 0 main" $dump_shape_file

# Output the local field for each cell
#for (( i=0; $i<$ncells; i++ ))
#do
#    add_dump "dump_cell_field $i $dump_cell_field_freq 0 main" "${dump_cell_field_file}_cell_${i}.dat" 
#done 

# Output the overlap field for each cell
#for (( i=0; $i<$ncells; i++ ))
#do
#    add_dump "dump_overlap_field $i $dump_overlap_field_freq 0 main" "${dump_overlap_field_file}_cell_${i}.dat" 
#done 
