#!/bin/bash

# Copyright 2022 Xanadu Quantum Technologies Inc.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#SBATCH --time=1:00:00
#SBATCH --job-name fp_sims2cloud_job
#SBATCH --mail-type=FAIL
#SBATCH --output=fp_sims2cloud_%j.out
#SBATCH  --cpus-per-task=1

# Initializing some default variables
cpu_count=2
AWS_PROFILE="default"
pyfile="default.py"
prefix="$(date +%s)"
orig_dir=`pwd`
python="python"

read -r -d '' help_message <<- EOM
	
	Copyright 2022 Xanadu Quantum Technologies Inc.
       	
	A simple BASH script to run already-parallelized FlamingPy's 
	Python scripts and automate writing/reading of simulations' 
	outputs and metadata from/to AWS S3 buckets. Note that this 
	script ASSUMES the simulation outputs are meant to be temporarily 
	and locally managed on 
	"<FTSTACKDIR>/flamingpy/sims_data/<prefix>_*.*". Please make 
	sure the intended file structure of FlamingPy is already checked 
	out by git or manually at <FTSTACKDIR> location.
	
	USAGE:
	
	  source sims2cloud.sh [OPTIONS as below]
	  
	EXAMPLE RUN (local systems and login nodes of the cloud):
        
	  source sims2cloud.sh --action "sims" --pyfile flamingpy/simulations.py --cpu 4 --ftstackdir .
	  source sims2cloud.sh --action "write" --aws-profile <xanadu_user> --prefix 1646469373 --ftstackdir .
	  source sims2cloud.sh --action "read" --aws-profile <xanadu_user> --prefix 1645633121 --ftstackdir .
	  
	EXAMPLE RUN (submitting a slurm job, e.g., valid for ComputeCanda's Niagara cluster with 32 tasks and 2 nodes totalling 64 PROCESSORs):
	
	  sbatch -n 32 -N 2 sims2cloud.sh --action "sims" --pyfile flamingpy/simulations.py --ftstackdir <FTSTACKDIR> --python <python>    
	
  OPTIONS:
	
	  --action <string>
	      Choose from "sims", "write", or "read". Note you may want 
	      to provide your own python interpreter with FlamingPy 
	      installed for "sims" purposes. Furthermore, you will 
	      need awscli package and a pre-configured aws-profile for 
	      "write" and "read" actions [REQUIRED].
	      
	  --python <path>
	      A python interpreter with FlamingPy installed 
	      [DEFAULT python].
	      
	  --pyfile <path>
	      Path to an accessible and compatible FlamingPy Python  
	      script for running simulations [DEFAULT default.py].
        
	  --cpu <integer>
	      Number of PROCESSORS to utilize for parallelization [DEFAULT 2].
	        
	  --prefix <string>
	      Name prefix to be used for storing outputs and metadata
	      [DEFAULT NAMING STYLE: [timestamp]_[pyfiletrimmed]_*.*].
	      
	  --aws-profile <string>
	      An existing aws profile in config file to be used 
	      [DEFAULT default].
	      
	  --ftstackdir <string>
	      An accessible location where FlamingPy is checked out 
	      by git or manually structured. [REQUIRED]
	      
EOM

if [ $# -eq 0 ]; then 
	echo "$help_message"
	echo
        return 1
fi    

while test $# -gt 0; do
    case "$1" in
        -h|--help)
            echo "$help_message"
	    echo
            return 1
            ;;
	--action)
            shift
            if test $# -gt 0; then
                export action=$1
            else
                echo "ERROR: no action was specified. Choose from \"sims\" \"write\", or \"read\""
                return -1
            fi
            shift
            ;;  
	--python)
            shift
            if test $# -gt 0; then
                export python=$1
            fi
            shift
            ;;    
        --pyfile)
            shift
            if test $# -gt 0; then
                export pyfile=$1
            fi
            shift
            ;;
        --cpu)
            shift
            if test $# -gt 0; then
                export cpu_count=$1
            fi
            shift
            ;;
        --prefix)
            shift
            if test $# -gt 0; then
                export prefix=$1
            fi
            shift
            ;; 
	--aws-profile)
            shift
            if test $# -gt 0; then
                export AWS_PROFILE=$1
            fi
            shift
            ;;
	--ftstackdir)
            shift
            if test $# -gt 0; then
                export FTSTACKDIR=$1
	    else
                echo "ERROR: no FTSTACKDIR was specified."
                return -1
            fi
            shift
            ;;
        *)
            break
            ;;
        esac
done	


cd ${FTSTACKDIR}


if [[ "$action" == "sims" ]]; then

	echo
	echo "\"sims\" action was selected. This scripts will now prepare for running FlamingPy" 
	echo "simulations."

        # Checking for FlamingPy
	echo
        echo "Checking if some version of FlamingPy is already installed ..."

	fp_package=`$python -m pip list | grep -F flamingpy`
	if [ ! -z "$fp_package" ]; then
		echo
		echo "	The following existing package will be utilized:" 
		echo "	${fp_package}"
	else
        	echo "ERROR: no flamingpy package as found."
		return -1
	fi
	
	if [[ "$HOSTNAME" = nia*.scinet.local ]] && [[ "$HOSTNAME" != nia-login*.scinet.local ]]; then
	  	echo
		echo "	ComputeCanda Niagara submission host was detected. Setting cpu_count ..." 
		module load gcc/9.4.0 openmpi/4.1.1
		mkdir -p mpl-config
		export MPLCONFIGDIR=./mpl-config
		cpu_count=${SLURM_ARRAY_TASK_COUNT}   
	fi

	echo
	echo "Running $pyfile simulations on $cpu_count CPU processors ..."
	#export MKL_NUM_THREADS=${cpu_count}
	#export OMP_NUM_THREADS=${cpu_count}
	mpirun --map-by node $python $pyfile
	
	# Finding filename from the full address
	pyfile_new=""
	pyfile_trimmed=$pyfile
	while [ "$pyfile_new" != "$pyfile_trimmed" ]; do
        	pyfile_new=${pyfile_trimmed}
		pyfile_trimmed=${pyfile_trimmed#*/}
	done	
	pyfile_trimmed=${pyfile_trimmed%.py*}
	prefix="${prefix}_${pyfile_trimmed}"
	
	## Rename simulations output based on prefix:
	cd flamingpy/sims_data/
	for f in ${pyfile_trimmed}_*.*; do mv "$f" $(echo "$f" | sed 's/^'"$pyfile_trimmed"'/'"$prefix"'/'); done

elif [[ "$action" == "write" ]]; then

	echo
	echo "\"write\" action was selected. Storing simulation outputs and metadata found in" 
	echo "\"src/flamingpy/sims_data\" to AWS S3 buckets using name prefix=$prefix."

        cd flamingpy/sims_data/

	awscli_package=`$python -m pip list | grep -F awscli`
	if [ -z "$awscli_package" ]; then $python -m pip install awscli; fi
	aws s3 cp ${prefix}_*.* s3://flamingpy-sims-bucket2/ --profile $AWS_PROFILE
	
elif [[ "$action" == "read" ]]; then	

	echo
	echo "\"read\" action was selected. Downloading simulation outputs and metadata from" 
	echo "AWS S3 buckets using name prefix=$prefix ..."

	cd flamingpy/sims_data/

	awscli_package=`$python -m pip list | grep -F awscli`
	if [ -z "$awscli_package" ]; then $python -m pip install awscli; fi

	aws s3 sync s3://flamingpy-sims-bucket2 . --exclude "*" --include "${prefix}_*.*" --profile $AWS_PROFILE
	
fi


cd ${orig_dir} 


echo
#return 0
