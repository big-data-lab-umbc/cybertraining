Owing to the volume of dataset, the analysis for this project has conducted on Jupyter notebooks in taki. To connect Jupyter notebook to your taki environment, please following the below-mentioned steps:

Step1: Connect to taki through terminal using your credentials

Step 2: In your taki user-space, type the following command:

  jupyter notebook --no-browser --port=8889

Running this, will prompt a localhost url similar to one below:

  http://localhost:8889/?token=cfccfc7e16af59f0d23be71bff7f729a91fcf5172672d3f3

Step3: In another terminal on your local machine, type the following with your taki username inplace of 'sali9':

  ssh -N -f -L localhost:8890:localhost:8889 sali9@taki.rs.umbc.edu 

Step4: In your browser, paste the link you get from Step 2 and replace the port with '8890', like this:

  http://localhost:8890/?token=cfccfc7e16af59f0d23be71bff7f729a91fcf5172672d3f3 (dummy URL)

You should be able to view Jupyter Lab interface in the browser if your connection is successful. Here you can move to different repositories and even upload/download notebooks directly to your taki space.
