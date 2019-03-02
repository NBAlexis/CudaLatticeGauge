# CudaLatticeGauge

lattice gauge simulation with CUDA

Testing on GTX-1060 (6G)

Using GMRES with 10 orthrognal basis, maximum supported (32 x 32 x 32 x 16)
For the release built
One trajectory with 50 Omelyan sub-steps cost about 470 secs
(The speed is related with kappa and beta, 470secs is tested using kappa=0.1355 and beta=2.5, 
GMRES with 10 orthrognal basis will converge with about 3 restarts)

See detailed.pdf

The default Code generation is "compute_61,sm_61";. To test on a GTX970m, please add "compute_52,sm_52;" . (NOTE, to support double float, must use a sm60+.)

put
CUDACXX=/usr/local/cuda-10.0/bin/nvcc
into /etc/environment so that it is applied when using sudo.
