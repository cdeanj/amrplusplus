Installation
------------

This section will help you get started with running the AmrPlusPlus pipeline with Nextflow and Docker. This tutorial assumes you will be running the pipeline from a POSIX compatible system such as Linux, Solaris, or OS X.

Setup
-----

We will go over a typical pipeline setup scenario in which you connect to a remote server, install Nextflow, and download the pipeline source code.

```bash
# username and host address (only necessary if you will be running the pipeline from a remove server)
$ ssh [USER]@[HOST]

# install Nextflow
$ curl -s https://get.nextflow.io | bash

# set executable permissions for the nextflow executable
$ chmod u+x nextflow

# move nextflow executable to a folder in your PATH environment variable
$ mv nextflow $HOME/bin

# create a test directory and change into it
$ mkdir amr_test && cd amr_test

# download pipeline source code
$ git clone https://github.com/cdeanj/amrplusplus .
```

Run a Simple Test
-----------------

We will run a small dataset that comes with the pipeline source code. As such, we will not be specifying any input paths as they have already been included. The first step, though, is to pull each Docker container from DockerHub. This will download each of the pipeline tool dependencies. Note, the download may halt for some of the larger containers. In that case, simply rerun the Docker pull command. As there are many tool dependencies, this could take some time depending on your connection speed.

```bash
# pull docker containers from Dockerhub (if this command times out, try running it once more)
$ docker pull colostatemeg/amrplusplus -a

# command to run the amrplusplus pipeline
$ nextflow run main.nf -profile docker --threads 4 --output test

# change directories to view pipeline outputs
$ cd test/
```
