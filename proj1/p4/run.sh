#!/usr/bin/env bash

NUM_OF_PROCESS=4
REDIRECTION="" #"> /dev/null"

if [ ! -z "$1" ]
then
	NUM_OF_PROCESS=$1
fi

echo ">> Run with $NUM_OF_PROCESS processes."
echo

echo "------- Using MPI_Scan -------"
eval "mpicc mpi_scan.c -o mpi_scan.out && mpirun -np \"$NUM_OF_PROCESS\" -hostfile host mpi_scan.out $REDIRECTION"
echo
sleep 3

echo "------- Using My Solution -------"
eval "mpicc my_prefix.c -o my_prefix.out && mpirun -np \"$NUM_OF_PROCESS\" -hostfile host my_prefix.out $REDIRECTION"
echo
