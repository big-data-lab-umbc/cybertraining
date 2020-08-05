Machine learning models for generating cloud properties

SERIAL VERSION (cnn.py): works, run with
    python cnn.py
Currently...
    - Generates COT for given cloud profiles in DATA/Out1
    - Important parameters:
        - l2r determines how many spatial "slices" to run, currently works with l2r between 1 and 500
        - num is number of cloud profiles to use
        - tpr is the percentage of profiles to use as training data, rest are testing data
    - Generates plots for comparison for last 3 testing data points

PARALLEL VERSION (par_cnn.py): works, runs with
    sbatch run.slurm
Currently...
    - Same as CNN, but splits data over spatial slices
    - In theory should work for all node, proc, size combos but it sometimes hangs (memory issues)
    - Don't change parameters in run.slurm for best results, should run in ~1min
