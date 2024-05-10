import numpy as np
from memory_usage import *
from time import time as runtimer
from read_functions import *
from mpi4py import MPI


# MPI data read/transfer test

comm = MPI.COMM_WORLD
# nproc = comm.Get_size()

datdir="/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/haiyan/memb_01/ctl/post/d02/"
t0=0
t1=30

nt, nz, nx1, nx2, pres = get_file_dims(datdir)
nt=t1
# Account for dropped edges
buffer = 80
nx1-=buffer*2
nx2-=buffer*2


# start = runtimer()

if comm.rank != 0:
    variable = np.ndarray((nt,nz,nx1,nx2), dtype=np.float64)

varname='T'
if comm.rank == 0:
    start = runtimer()
    variable = np.ascontiguousarray(var_read_3d(datdir,varname,t0,t1,mask=True,drop=True), dtype=np.float64)
    end = runtimer()
    time_elapsed = end - start
    print("Time elapsed part1: ", time_elapsed)

start = runtimer()

comm.Bcast(variable, root=0)

end = runtimer()
time_elapsed = end - start
print("Time elapsed part2: ", time_elapsed)


# Memory test

# print()
# print("Before function:")
# memory_usage()

# def init_vars():
#     dims = 500000000
#     var1 = np.arange(dims, dtype=np.float64)
#     var2 = np.arange(dims, dtype=np.float64)
#     var3 = np.arange(dims, dtype=np.float64)
#     print()
#     print("Inside function:")
#     memory_usage()
#     return

# init_vars()

# print()
# print("After function:")
# memory_usage()
