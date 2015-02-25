Requirements:
    - UCSC tools installed on local machine, not too old (make sure that your
      liftUp command supports the bed8 format for -type), at least the tools:
      faSize, pslToBed, liftUp, 
    - UCSC source code on local machine 
    - variable KENTSRC which points to the base directory of the UCSC srccode
      e.g. if you are using bourne-like shells with the command
        export KENTSRC=/usr/local/src/kent/src
    - an SGE compute cluster 
    - python (at least version 2.5), gawk, cut

How to run the pipeline (e.g. for Drosophila pseudoobscura):

- edit config.mk and edit username, name of your compute cluster and the target
  directory on it.
- copy your ssh public key to the cluster (if it's not there yet)
    make clusterPublicKey
  this avoids passwords. Run ssh-keygen and retry if you do not have a public
  key yet.
- put the file dm3.dp3.synteny.chain from the supplemental file 2 into the
  directory engstromPipeline/dm3dp3Example/
- download the genomes from the UCSC download server:
    cd engstromPipeline
    make getGenomes
- split the genome into pieces according to the chain file:
    cd ../fastaFragments
    make
        (this can take several minutes)
    You can inspect the output of this step with
        less dm3dp3Example/fragments.bed
        less dm3dp3Example/chr2L-25*.fa

- run contemplate on the pieces:
    cd ../contemplate
        (edit the makefile to point to the target SGE compute cluster
        directory, variables CLUSTERHOST and CLUSTERDIR)
    make toCluster
        ( this copies all files to the cluster)
    log into your cluster, go to the target directory defined in 
    config.mk and run:
        make submitJobs
    When the jobs are finished, copy them back with:
        make convert
- Inspect the result with
        less dm3dp3Example/constrainedElements.filtered.minBlocks3.bed

Description of files and directories:

config.mk:
    central config file which specifies the name of the current dataset to
    use and its parameters. The name of the dataset will be used to create 
    subdirectories in all of the following directories. All contemplate-
    related data will be stored in these subdirectories. This scheme allows
    us to switch between contemplate runs by modifying only config.mk

engstromPipeline:
    modified sourcecode of Per Engstrom's synteny pipeline. 
    cd into this directory, run make in all
    subdirectories of src/, then run "make" in this directory to generate
    a list of syntenic regions between dm3 and droVir3.

    I have modified the chainNetSynteny program to output
    the native chain format instead of the processed 
    bed file.

fastaFragments:
    a C program based on the UCSC toolchain which gets FASTA fragments for
    chains and cuts them into suitable fragments. Ouput is fasta, bed (to
    inspect where it cut) and a liftUp file (to lift fragment annotations to
    the chromosome level).
    
    cd into this, run "make compile", then
    "make run"

contemplate:
    The contemplate software (biojava) written by David Huen with input from
    Thomas Down and Max Haeussler.

    Adapt the makefile with the name of your cluster and the directory where you
    want to store the files there.

    Change into the "contemplate" directory, run "make toCluster", ssh onto your
    cluster, run "make submitJobs", go back to your machine, run "make
    fromCluster" and "make liftUp"

    If you have a local UCSC mirror, you can load the results into it with
    "make loadTracks" (assuming that your ~/.hg.conf is working)

    check the log directory on your cluster for "Exception" messages.

