To build the cython code, use the following command in the terminal:
CC=gcc-version python BlueCrystal_setup.py build_ext -fi

Two different setup files are provided: mac_setup.py and BlueCrystal_setup.py


In the run_script.py file, if save_matrix = True, then the position matrices will be 
Stored locally. This matrix is used to generate the animations using animations.py

The run_script.py can be executed by using the following command in the terminal  
python ./run_script.py {number of threads}


To run run_script.py on BlueCrystal, use the cython.sh make file provided in the folder and submit the job.  


