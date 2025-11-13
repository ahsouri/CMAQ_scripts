#!/bin/bash

# CMAQ Daily Loop Script with SLURM dependencies
# This script submits CMAQ jobs for each day, with dependencies between days

# Set your date range
START_DATE="2024-05-26"
END_DATE="2024-06-01"
CMAQ_RUN_SCRIPT="./run_cctm_TEMPO.csh"  # Path to your CMAQ run script

# Create directories for organization
SCRIPT_DIR="./cmaq_daily_scripts"
LOG_DIR="./cmaq_daily_logs"
mkdir -p $SCRIPT_DIR $LOG_DIR

# Create a backup of the original script
BACKUP_SCRIPT="${CMAQ_RUN_SCRIPT}.backup_$(date +%Y%m%d_%H%M%S)"
cp $CMAQ_RUN_SCRIPT $BACKUP_SCRIPT
echo "Created backup: $BACKUP_SCRIPT"

# Initialize variables
current_date=$START_DATE
is_first_day=false
previous_job_id=""

echo "Setting up CMAQ daily jobs from $START_DATE to $END_DATE"

# Loop through each day
while [[ "$current_date" < $(date -I -d "$END_DATE + 1 day") ]]; do
    
    echo "Preparing job for date: $current_date"
    
    # Calculate next day for END_DATE
    next_date=$(date -I -d "$current_date + 6 day")
    # Create batch script for this day
    batch_script="${SCRIPT_DIR}/run_cmaq_${current_date}.sh"
    
    # Create the SLURM batch script
    cat > $batch_script << EOF
#!/bin/bash

# CMAQ batch script for $current_date

#SBATCH -J CMAQ_${current_date}
##SBATCH --qos=debug
##SBATCH --constraint="mil"
#SBATCH --account=s3223
#SBATCH --ntasks=480
##SBATCH --cpus-per-task=120
##SBATCH --mem=24G
#SBATCH -t 6:30:00
#SBATCH -o ${LOG_DIR}/cmaq_${current_date}-%j.out
#SBATCH -e ${LOG_DIR}/cmaq_${current_date}-%j.err

# Load software
module load comp/intel/2021.4.0
module load wrf-deps/1.0
module load mpi/impi/2021.4.0

# Run CMAQ script for $current_date
./run_cctm_TEMPO_${current_date}.csh
EOF
    
    # Create a working copy of the CMAQ script from the original backup
    working_script="./run_cctm_TEMPO_${current_date}.csh"
    cp $BACKUP_SCRIPT $working_script
    
    # Set NEW_START based on whether this is the first day
    if [ "$is_first_day" = true ]; then
        # First day - set NEW_START to TRUE
        sed -i '/setenv NEW_START/c\setenv NEW_START TRUE             #> Set to FALSE for model restart' $working_script
        echo "  -> First day: NEW_START set to TRUE"
        is_first_day=false
    else
        # Subsequent days - set NEW_START to FALSE  
        sed -i '/setenv NEW_START/c\setenv NEW_START FALSE            #> Set to FALSE for model restart' $working_script
        echo "  -> Continuation day: NEW_START set to FALSE"
    fi
    
    # Update the dates - START_DATE is current day, END_DATE is next day
    sed -i '/set START_DATE =/c\set START_DATE = "'$current_date'"     #> beginning date' $working_script
    sed -i '/set END_DATE   =/c\set END_DATE   = "'$next_date'"     #> ending date' $working_script
    
    # Make both scripts executable
    chmod +x $batch_script
    chmod +x $working_script
    
    # Create SLURM submission command
    if [ -z "$previous_job_id" ]; then
        # First job - no dependency
        echo "  -> Submitting first day job with no dependencies"
        submit_cmd="sbatch ${batch_script}"
    else
        # Subsequent jobs - depend on previous day
        echo "  -> Submitting job with dependency on previous job (${previous_job_id})"
        submit_cmd="sbatch --dependency=afterok:${previous_job_id} ${batch_script}"
    fi
    
    # Submit the job and capture job ID
    echo "  -> Running: $submit_cmd"
    submit_output=$(eval $submit_cmd)
    job_id=$(echo $submit_output | awk '{print $4}')
    previous_job_id=$job_id
    
    echo "  -> Submitted job for $current_date (runs until $next_date) with Job ID: $job_id"
    
    # Move to next day
    current_date=$(date -I -d "$current_date + 7 day")
    echo ""
    
done

echo "All CMAQ jobs have been submitted with appropriate dependencies!"
echo "Batch scripts located in: $SCRIPT_DIR"
echo "CMAQ run scripts in current directory: run_cctm_TEMPO_YYYY-MM-DD.csh"
echo "Job logs will be in: $LOG_DIR"
echo "Original script backup preserved at: $BACKUP_SCRIPT"
