threads=64;
final_core=$((threads-1));
# Section 6.3: Mode: Custom Gates
#Table 5

taskset -c 0-"$final_core" ./bin/main 63/stackexchange/1.txt -o 63/stackexchange/1 -gs stackexchange/1 --threads "$threads" --save --circuits 100 --time 43200 -cr 1.999 --depth-gates U 
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/2.txt -o 63/stackexchange/2 --threads "$threads" --expand --save --circuits 100 --time 43200 -tc 3
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/3.txt -o 63/stackexchange/3 --threads "$threads" --expand --save --circuits 100 -tc 7 --time 43200
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/4.txt -o 63/stackexchange/4 -gs stackexchange/4 --threads "$threads" --save --circuits 100 --time 43200 -cr 2.999 --depth-gates scz
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/5.txt -o 63/stackexchange/5 --threads "$threads" --expand --save --time 43200 --circuits 100 -tc 4

# Section 6.4: Mode: Clifford+T Gate Set

taskset -c 0-"$final_core" ./bin/main 64/comparison/ccx.txt -o 64/comparison/tcount/ccx --threads "$threads" --circuits 100 -tc 7 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/ccx.txt -o 64/comparison/tdepth/ccx --threads "$threads" --circuits 100 -td 3 --time 43200 --save

# Section 6.5: Mode: Approximate Circuit Synthesis

taskset -c 0-"$final_core" ./bin/main 65/rz_4.txt -o 65/005/rz_4 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/rz_5.txt -o 65/005/rz_5 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05

taskset -c 0-"$final_core" ./bin/main 65/crz_2.txt -o 65/005/crz_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/crz_2.txt -o 65/10_2/crz_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01

taskset -c 0-"$final_core" ./bin/main 65/crk_2.txt -o 65/005/crk_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/crk_2.txt -o 65/10_2/crk_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01

# Section 6.6: Mode: Relative Phase Operators

taskset -c 0-"$final_core" ./bin/main 66/carry/carry_relative.txt -o 66/carry_relative --threads "$threads" --expand --save --circuits 100 -tc 4 -td 2 --time 1800 -cg composite_rccx
taskset -c 0-"$final_core" ./bin/main 66/carry/carry_relative.txt -o 66/carry_relative_cnot --threads "$threads" -gs CliffordTCNOT --expand --save --circuits 1000 --time 1800 --depth-gates cx