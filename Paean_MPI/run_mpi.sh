#!/bin/bash

echo "#!/bin/bash" > Paean.slurm
echo >> Paean.slurm
echo "#SBATCH -J Paean" >> Paean.slurm
echo "#SBATCH -N "$1 >> Paean.slurm
echo "#SBATCH -n "`expr $1 \* 4` >> Paean.slurm
echo "#SBATCH -p huge" >> Paean.slurm
echo "#SBATCH -c 8" >> Paean.slurm
echo "#SBATCH --gres=dcu:4" >> Paean.slurm
echo "#SBATCH -o log/%J.out" >> Paean.slurm
echo "#SBATCH -e log/%J.err" >> Paean.slurm
echo >> Paean.slurm
echo "for i in \`scontrol show hostnames \$SLURM_JOB_NODELIST\`" >> Paean.slurm
echo "do" >> Paean.slurm
echo "  echo "\$i slots=4"" >> Paean.slurm
echo "done > hostfile" >> Paean.slurm
echo >> Paean.slurm
echo "mpirun -np "`expr $1 \* 4`" --hostfile \`pwd\`/hostfile --mca plm_rsh_no_tree_spawn 1 --mca plm_rsh_num_concurrent "$1" -mca routed_radix "$1" -mca pml ucx -x LD_LIBRARY_PATH -mca coll_hcoll_enable 0 --bind-to none \`pwd\`/mpi_bind.sh "$2 >> Paean.slurm

sbatch -x k06r2n17 Paean.slurm
