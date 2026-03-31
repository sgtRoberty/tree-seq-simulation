#!/bin/bash

NSIMS=10

# Argument parsing
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --nsims)
      NSIMS="$2"
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [--nsims N]"
      exit 0
      ;;
    *)
      echo "Unknown parameter passed: $1"
      exit 1
      ;;
  esac
  shift
done

echo "Number of simulations: $NSIMS"

# Loop through replicates
for ((i=1; i<=NSIMS; i++)); do
  rep_dir="rep$i"

  if [ ! -d "$rep_dir" ]; then
    echo "Directory $rep_dir not found, skipping..."
    continue
  fi

  cd "$rep_dir"

  job_script="job_rep${i}.sh"

  # 1. Create SLURM job script
  cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=job-rep${i}
#SBATCH --output=job-rep${i}-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:05:00

beast -overwrite -D index=${i} -df var.json -DFout run.out.xml ../XML/run.xml
EOF

  chmod +x "$job_script"

  echo "Submitting $job_script ..."
  sbatch "$job_script"

  cd ..
done