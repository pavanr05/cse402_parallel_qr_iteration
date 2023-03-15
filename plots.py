from cProfile import label
import matplotlib.pyplot as plt
import numpy as np

#Speedups
size = np.array([16, 32, 64, 128, 256, 512])
serial_time = np.array([0.00398183, 0.041914, 0.437767, 6.54397, 104.739, 1629.93])
omp_time = np.array([0.108697, 0.225223, 0.472626, 1.54925, 11.5782, 170.722])
mpi_time = np.array([0.0151141, 0.092443, 0.352974, 2.79907, 22.825, 218.502])

omp_speedup = np.divide(serial_time, omp_time)
mpi_speedup = np.divide(serial_time, mpi_time)

plt.figure(1)
plt.plot(size, omp_speedup, label="OpenMP (20 threads)")
plt.plot(size, mpi_speedup, label="MPI (16 processes)")
plt.xlabel("Matrix Size")
plt.ylabel("Speedup")
plt.title("Speedup in OpenMP and MPI")
plt.legend()
plt.savefig("speedup.png")
plt.show()


#OpenMP scaling
omp_ss_time = np.array([1527.86, 688.164, 387.963, 202.423, 125.593, 172.783, 150.778, 153.028])
nthreads = np.array([1, 2, 4, 8, 16, 20, 32, 64])
ideal_ss = np.array(1527.86/nthreads)

plt.figure(2)
plt.loglog(nthreads, omp_ss_time, label="Observed scaling")
plt.loglog(nthreads, ideal_ss, '--', label = "Ideal Scaling")
plt.xlabel("Number of threads")
plt.ylabel("Time (s)")
plt.title("OpenMP scaling using (512 x 512) matrix")
plt.legend()
plt.savefig("omp_ss.png")
plt.show()

#MPI scaling
mpi_ss_time = np.array([1279.92, 665.061, 392.332, 254.447, 218.652])
num_procs = np.array([1, 2, 4, 8, 16])
ideal_ss_mpi = np.array(1279.92/num_procs)
plt.figure(3)
plt.loglog(num_procs, mpi_ss_time, label="Observed scaling")
plt.loglog(num_procs, ideal_ss_mpi, '--', label = "Ideal Scaling")
plt.xlabel("Number of processes")
plt.ylabel("Time (s)")
plt.title("MPI scaling using (512 x 512) matrix")
plt.legend()
plt.savefig("mpi_ss.png")
plt.show()
