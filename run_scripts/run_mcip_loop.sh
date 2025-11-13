#!/bin/bash

# User settings
start_date="2024-06-01"
end_date="2024-08-01"
InMetDir="/discover/nobackup/asouri/MODELS/CMAQv5.5/data/wrf/"
mcip_script="./run_mcip.csh"

current_date="$start_date"
end_date_plus1=$(date -I -d "$end_date + 1 day")

while [[ "$current_date" != "$end_date_plus1" ]]; do
    echo $current_date
    day_before=$(date -I -d "$current_date - 1 day")
    day_after=$(date -I -d "$current_date + 1 day")
    yyyymmdd=$(date -d "$current_date" +%Y%m%d)
    # Format the dates properly for MCIP
    mcip_start="${current_date}-00:00:00.0000"
    mcip_end="${day_after}-00:00:00.0000"
    # Use awk to replace the block
    awk -v mcip_start="$mcip_start" -v mcip_end="$mcip_end" -v  current_date="$current_date" -v day_before="$day_before" -v day_after="$day_after" -v InMetDir="$InMetDir" '
    /# BEGIN: AUTO-MODIFIED/ {
        print "# BEGIN: AUTO-MODIFIED"
        print "set MCIP_START = " mcip_start
        print "set MCIP_END = " mcip_end
        print ""
        print "set InMetFiles = ( " InMetDir "/wrfout_d01_" day_before "_00:00:00 \\"
        print "                   " InMetDir "/wrfout_d01_" current_date "_00:00:00 \\"
        print "                   " InMetDir "/wrfout_d01_" day_after "_00:00:00 )"
        print "# END: AUTO-MODIFIED"
        # Skip until we find the END marker
        while (getline > 0 && !/# END: AUTO-MODIFIED/) continue
        next
    }
    {print}
    ' "$mcip_script" > "${mcip_script}.tmp" && mv "${mcip_script}.tmp" "$mcip_script"
    
    # Replace APPL with DOY string
    sed -i "s/^set APPL = .*/set APPL = CONUS_8km_${yyyymmdd}/" "$mcip_script"
    
    echo ">>> Running MCIP for ${current_date}..."
    chmod a+x ./run_mcip.csh
    ./run_mcip.csh

    current_date=$(date -I -d "$current_date + 1 day")
done
