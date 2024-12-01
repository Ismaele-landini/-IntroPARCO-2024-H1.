# -IntroPARCO-2024-H1.
 Exploring Implicit and Explicit Parallelism with OpenMP.

# -Requirements
GNU Compiler: Ensure that versions 5.4, 7.5, or 9.1 are available on the cluster.
SSH Client: Use PuTTY or MobaXterm for Windows. For macOS and Linux, the built-in SSH client is sufficient.
VPN: To access the University network from an external network, establish a secure connection using the Virtual Private Network (VPN).

# -Instructions for Compilation and Execution

1) Download Files: Obtain the instructions.pbs and parallel_matrix_transposition.c files. These will be used for job submission and program execution.
2) Access the HPC Cluster:
   a) Use a VPN to establish a secure connection to the Trento University network.
   b) Open your SSH client and connect to the cluster with the following command:
     ssh username@hpc.unitn.it
Enter your university credentials when prompted.
3) Upload Files: Navigate to your desired directory on the cluster or create a new one.
Upload instructions.pbs and parallel_matrix_transposition.c to the chosen directory using the file transfer feature of your SSH client (e.g., drag-and-drop in MobaXterm or scp for Linux/macOS).
4) Reserve a Node and Enter an Interactive Session:
   a) Move to the directory containing your files:
   cd path_to_your_directory
   b) Request an interactive session on a node with the following specifications: 64 cores, 64
    OpenMP threads and 1 MB of memory. Submit the request using:
    [username@hpc-head-n1 ~]$ qsub -I instructions.pbs
5) Modify the Program Input (Optional): once inside the node, you can modify the instructions.pbs file to change the matrix size. Open the file in a text editor and locate the section marked, instead of the default number (128):
    ./parallel_matrix_transposition 128 >> $OUTPUT_FILE 
6) Submit the job to the queue for compile and execution:
   [username@hpc-head-n1 ~]$ qsub instructions.pbs
7) View the Results: once the program completes, its output will be saved in a file named results.txt in the same directory.

# -View the paper
Download the Exploring_Implicit_and_Explicit_Parallelism_with_OpenMP.pdf
