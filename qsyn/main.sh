i=0;
threads=64;
mkdir -p output
for file in benchmarks/*   
do  
    filename=$(basename "$file");
    name=${filename%.json};
    echo "$name";
    docker run qsyn taskset -c "$i" python3 qsyn/run_single.py --benchmark "$name" --mode Ours > output/"$name".txt &
    i=$((i+1));
    if [ "$i" -eq "$threads" ]
    then
        wait
        i=0;
    fi
done