mpic++ -lm main_mpi.c lab4_mpi.cpp lab4_io.c -o ppm
mpic++ -lm main_mpi.c lab4_mpi_bf.cpp lab4_io.c -o ppmbf

# cd testcase
# python gen_testcase.py
# cd ..

time mpirun -np 8 ppm $1 > out.txt
time mpirun -np 4 ppmbf $1 > outbf.txt

diff out.txt outbf.txt