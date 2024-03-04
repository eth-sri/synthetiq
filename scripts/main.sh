declare -A counts=( [three_superpose]=3 [M_valued]=7 [GHZ_from_100]=6 [GHZ_by_iSWAP]=8 [GHZ_by_QFT]=7 [GHZ_Game]=8 [W_orthog]=8 [W_phased]=7 [W_four]=7 [cluster]=13 [bit_measure]=7 [flip]=6 [teleportation]=4 [indexed_bell]=4 [toffoli_by_sqrt_X]=5 [QFT]=7 [draper]=5 );
declare -A optimal_counts=( [three_superpose]=3 [M_valued]=6 [GHZ_from_100]=6 [GHZ_by_iSWAP]=8 [GHZ_by_QFT]=5 [GHZ_Game]=8 [W_orthog]=6 [W_phased]=7 [W_four]=6 [cluster]=9 [bit_measure]=7 [flip]=6 [teleportation]=4 [indexed_bell]=4 [toffoli_by_sqrt_X]=5 [QFT]=7 [draper]=5 );
i=0;
threads=64;
final_core=$((threads-1));
# Section 6.1: Hyperparameter Optimization
python scripts/optimizing.py --cores "$threads"

# Section 6.2: Better Operator Decompositions
taskset -c 0-"$final_core" ./bin/main 62/rcccx.txt -o 62/rcccx --threads "$threads" --circuits 10000 --time 28800 --save -cg composite_rccx --expand
taskset -c 0-"$final_core" ./bin/main 62/cct.txt -o 62/cct --threads "$threads" --circuits 10000 --time 7200 --save -cg composite_rccx -a 1 --expand
taskset -c 0-"$final_core" ./bin/main 62/cciswap.txt -o 62/cciswap --threads "$threads" --circuits 10000 --time 14400 --save -cg composite_ccx --expand
taskset -c 0-"$final_core" ./bin/main 62/csqrtswap.txt -o 62/csqrtswap --threads "$threads" --circuits 10000 --time 14400 --save -cg composite_rccx -a 1 --expand
taskset -c 0-"$final_core" ./bin/main 62/csqrtiswap.txt -o 62/csqrtiswap --threads "$threads" --circuits 10000 --time 14400 --save -cg composite_rccx --expand

# Section 6.3: Mode: Custom Gates

#Figure 3
for file in ./data/input/63/comparison/*   
do  
    filename=$(basename "$file");
    name=${filename%.txt};
    echo "$name";
    taskset -c "$i" ./bin/main 63/comparison/"$name".txt -o 63/theirs/"$name" --threads 1 --circuits 100 --time 10800 -gs composite_gatesets/"$name" --gates 4 --save -gc "${counts[$name]}"  &
    i=$((i+1));
    if [ "$i" -eq "$threads" ]
    then
        wait
        i=0;
    fi
done

for file in ./data/input/63/comparison/*   
do  
    filename=$(basename "$file");
    name=${filename%.txt};
    taskset -c "$i" ./bin/main 63/comparison/"$name".txt -o 63/optimal/"$name" --threads 1 --circuits 100 --time 10800 -gs composite_gatesets/"$name" --gates 4 --save -gc "${optimal_counts[$name]}"  &
    i=$((i+1));
    if [ "$i" -eq "$threads" ]
    then
        wait
        i=0;
    fi
done

wait

taskset -c 0-"$final_core" ./bin/main 63/stackexchange/1.txt -o 63/theirs/1 -gs stackexchange/1 --threads "$threads" --save --circuits 100 --time 43200 -gc 9
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/2.txt -o 63/theirs/2 --threads "$threads" --expand --save --circuits 100 --time 43200 -gc 9
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/3.txt -o 63/theirs/3 --threads "$threads" --expand --save --circuits 100 -gc  19 --time 43200
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/4.txt -o 63/theirs/4 -gs stackexchange/4 --threads "$threads" --save --circuits 100 --time 43200 -gc 14
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/5.txt -o 63/theirs/5 --threads "$threads" --expand --save --time 43200 --circuits 100 -gc 13

#Table 5

taskset -c 0-"$final_core" ./bin/main 63/stackexchange/1.txt -o 63/stackexchange/1 -gs stackexchange/1 --threads "$threads" --save --circuits 100 --time 43200 -cr 1.999 --depth-gates U 
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/2.txt -o 63/stackexchange/2 --threads "$threads" --expand --save --circuits 100 --time 43200 -tc 3
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/3.txt -o 63/stackexchange/3 --threads "$threads" --expand --save --circuits 100 -tc 7 --time 43200
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/4.txt -o 63/stackexchange/4 -gs stackexchange/4 --threads "$threads" --save --circuits 100 --time 43200 -cr 2.999 --depth-gates scz
taskset -c 0-"$final_core" ./bin/main 63/stackexchange/5.txt -o 63/stackexchange/5 --threads "$threads" --expand --save --time 43200 --circuits 100 -tc 4

# Section 6.4: Mode: Clifford+T Gate Set

taskset -c 0-"$final_core" ./bin/main 64/comparison/ccx.txt -o 64/comparison/tcount/ccx --threads "$threads" --circuits 100 -tc 7 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/cch.txt -o 64/comparison/tcount/cch --threads "$threads" --circuits 100 -tc 9 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/adder.txt -o 64/comparison/tcount/adder --threads "$threads" --circuits 100 -tc 7 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/U1_var.txt -o 64/comparison/tcount/U1_var --threads "$threads" --circuits 100 -tc 11 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/U1.txt -o 64/comparison/tcount/U1 --threads "$threads" --circuits 100 -tc 11 --time 43200 --save  
taskset -c 0-"$final_core" ./bin/main 64/comparison/U2.txt -o 64/comparison/tcount/U2 --threads "$threads" --circuits 100 -tc 7 --time 43200 --save  


taskset -c 0-"$final_core" ./bin/main 64/comparison/ccx.txt -o 64/comparison/tdepth/ccx --threads "$threads" --circuits 100 -td 3 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/cch.txt -o 64/comparison/tdepth/cch --threads "$threads" --circuits 100 -td 4 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/adder.txt -o 64/comparison/tdepth/adder --threads "$threads" --circuits 100 -td 2 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/U1_var.txt -o 64/comparison/tdepth/U1_var --threads "$threads" --circuits 100 -td 5 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/U1.txt -o 64/comparison/tdepth/U1 --threads "$threads" --circuits 100 -td 5 --time 43200 --save
taskset -c 0-"$final_core" ./bin/main 64/comparison/U2.txt -o 64/comparison/tdepth/U2 --threads "$threads" --circuits 100 -td 3 --time 43200 --save

for i in {0,}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tdepth/"$i" --threads "$threads" --circuits 100 --time 7200 --save -td 0 --expand
done

for i in {1,2,5,6,13,22,29}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tdepth/"$i" --threads "$threads" --circuits 100 --time 7200 --save -td 3 --expand
done

for i in {3,4,7,8,11,12,15,16,19,20,23,24,27,28}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tdepth/"$i" --threads "$threads" --circuits 100 --time 10800 --save -td 3 --expand
done

for i in {9,10,14,17,18,21,25,26}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tdepth/"$i" --threads "$threads" --circuits 100 --time 7200 --save -td 7 -cg composite_rccx --expand
done

for i in {0,}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tcount/"$i" --threads "$threads" --circuits 100 --time 7200 --save -tc 0 --expand  
done

for i in {1,2,5,6,13,22,29}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tcount/"$i" --threads "$threads" --circuits 100 --time 7200 --save --save -tc 7 --expand  
done

for i in {3,4,7,8,11,12,15,16,19,20,23,24,27,28}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tcount/"$i" --threads "$threads" --circuits 100 --time 7200 --save -tc 8 --expand  
done

for i in {9,10,14,17,18,21,25,26}
do    
    echo "$i";
    taskset -c 0-"$final_core" ./bin/main 64/permutations/"$i".txt -o 64/permutations/tcount/"$i" --threads "$threads" --circuits 100 --time 7200 --save -tc 15 -cg composite_rccx --expand  
done

# Section 6.5: Mode: Approximate Circuit Synthesis

taskset -c 0-"$final_core" ./bin/main 65/rz_4.txt -o 65/005/rz_4 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/rz_5.txt -o 65/005/rz_5 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/rz_4.txt -o 65/10_2/rz_4 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01
taskset -c 0-"$final_core" ./bin/main 65/rz_5.txt -o 65/10_2/rz_5 --threads "$threads" --circuits 100 --time 21600 --save -eps 0.01
taskset -c 0-"$final_core" ./bin/main 65/rz_6.txt -o 65/10_2/rz_6 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01
taskset -c 0-"$final_core" ./bin/main 65/rz_7.txt -o 65/10_2/rz_7 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01

taskset -c 0-"$final_core" ./bin/main 65/crz_2.txt -o 65/005/crz_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/crz_2.txt -o 65/10_2/crz_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01

taskset -c 0-"$final_core" ./bin/main 65/crk_2.txt -o 65/005/crk_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.05
taskset -c 0-"$final_core" ./bin/main 65/crk_2.txt -o 65/10_2/crk_2 --threads "$threads" --circuits 100 --time 3600 --save -eps 0.01

taskset -c 0-"$final_core" ./bin/main 65/ccrz_2.txt -o 65/10_3/ccrz_2 --threads "$threads" --circuits 100 --time 7200 --save -eps 0.001

# Section 6.6: Mode: Relative Phase Operators

taskset -c 0-"$final_core" ./bin/main 66/carry/carry_relative.txt -o 66/carry_relative --threads "$threads" --expand --save --circuits 100 -tc 4 -td 2 --time 43200 -cg composite_rccx
taskset -c 0-"$final_core" ./bin/main 66/carry/carry_relative.txt -o 66/carry_relative_cnot --threads "$threads" -gs CliffordTCNOT --expand --save --circuits 1000 --time 43200 --depth-gates cx