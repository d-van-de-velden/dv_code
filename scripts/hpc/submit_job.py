"""Functions to grab a function, prepare it, and send it to HPC cluster for execution."""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#          Alexander Enge          (XX@YY.de)
#          Anne-Sophie Kieslinger  (XX@YY.de)
#
# License: BSD (3-clause)
# 
# 
#

from simple_slurm import Slurm
from dv_code.scripts.misc import check_make_dir

def submit_job(args_list, cpus=2, mem=32000, time='24:00:00', log_dir='logs/',
            dependency_jobs=[], dependency_type='afterok', job_name='job'):
    
    """Submits a single batch job via SLURM, which can depend on other jobs.

    Parameters
    ----------
    args_list : list
        A list of shell commands and arguments. The first element will usually
        be the path of a shell script and the following elements the input
        arguments to this script.
    cpus : int, default=8
        The number of CPUs that the batch job should use.
    mem : int, default=320000
        The amount of memory (in MB) that the batch job should use.
    time : str, default='24:00:00'
        The maximum run time (in format 'HH:MM:SS') that the batch job can use.
        Must not exceed 24 hours.
    log_dir : str or Path, default='logs/'
        Directory to which the standard error and output messages of the batch
        job should be written.
    dependency_jobs : int or list, default=[]
        Other SLURM batch job IDs on which the current job depends. Can be used
        to create a pipeline of jobs that are executed after one another.
    dependency_type : str, default='afterok
        How to handle the 'dependency_jobs'. Must be one of ['after',
        'afterany', 'afternotok', 'afterok', 'singleton']. See [1] for further
        information. 
    job_name : str, default='job'
        Name of the slurm job that will submitted. 

    Returns
    -------
    job_id : int
        The job ID of the submitted SLURM batch job.

    Notes
    -----
    [1] https://hpc.nih.gov/docs/job_dependencies.html
    """


    # Workflow:
    # 1. Check out function and dependencies
    # 2. Check if container for said function is already available
    # 3. Create container if not ready
    # 4. Create script that executes command in container on HPC cluster
    # 5. Send container to cluster
    # 6. Execute 

    # https://apptainer.org/docs/user/latest/build_a_container.html
    # https://apptainer.org/docs/user/latest/definition_files.html
    
    
    # Join arguments to a single bash command
    cmd = ' '.join(str(arg) for arg in args_list)

    # Create directory for output logs
    log_dir = check_make_dir(log_dir)
    
    error = f'{log_dir}/slurm-%j-{job_name}.out'
    output = f'{log_dir}/slurm-%j-{job_name}.out'

    # Prepare job scheduler
    slurm = Slurm(cpus_per_task=cpus, error=error, mem=mem, nodes=1, ntasks=1,
                    output=output, time=time, job_name=job_name)

    # Make the current job depend on previous jobs
    if dependency_jobs != []:
        if isinstance(dependency_jobs, int):
            dependency_jobs = [dependency_jobs]
        dependency_str = ':'.join([str(job_id) for job_id in dependency_jobs])
        dependency = {dependency_type: dependency_str}
        slurm.set_dependency(dependency)

    # Submit
    print('Submitting', cmd)
    job_id = slurm.sbatch(cmd)

    return job_id