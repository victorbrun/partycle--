#!/bin/bash

# Changing dir to build/ to get all the file access stuff correct
cd ./build/

echo "Benching 209 particles"
hyperfine './bin/partycle-- -d [0,10]x[0,10]x[0,10] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_209p.hyperfine

echo "Benching 1676 particles"
hyperfine './bin/partycle-- -d [0,20]x[0,20]x[0,20] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_1676p.hyperfine

echo "Benching 5658 particles"
hyperfine './bin/partycle-- -d [0,30]x[0,30]x[0,30] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_5658p.hyperfine

echo "Benching 13413 particles"
hyperfine './bin/partycle-- -d [0,40]x[0,40]x[0,40] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_13413p.hyperfine

#echo "Benching 26198 particles"
#hyperfine './bin/partycle-- -d [0,50]x[0,50]x[0,50] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_26198p.hyperfine

echo "Benching 45272 particles"
hyperfine './bin/partycle-- -d [0,60]x[0,60]x[0,60] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_45272p.hyperfine

echo "Benching 71890 particles"
hyperfine './bin/partycle-- -d [0,70]x[0,70]x[0,70] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_71890p.hyperfine

echo "Benching 107311 particles"
hyperfine './bin/partycle-- -d [0,80]x[0,80]x[0,80] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_107311p.hyperfine

echo "Benching 152794 particles"
hyperfine './bin/partycle-- -d [0,90]x[0,90]x[0,90] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_152794p.hyperfine

echo "Benching 209594 particles"
hyperfine './bin/partycle-- -d [0,100]x[0,100]x[0,100] -cf ../tests/example_particle_distributions.csv -ct 1e-2' > bench_209594p.hyperfine
