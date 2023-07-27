for numModalities in 1; do #{2..6}; do
    for nsamples in 1; do #{1..5}; do
        for seed in 1; do #{1..3}; do
            prefix=${seed}'_'${nsamples}'_'${numModalities};
            echo $prefix;

            modeArr=""
            for ((i=0; i<numModalities; i++)); do
                modeArr+="5 "
            done

            python3 /scratch/data/nsdong2/projectPACTION/newpaction/src/simulation.py \
                -s ${seed} -n ${nsamples} -t 0 -p 0.2 \
                -m $modeArr \
                -o /scratch/data/nsdong2/projectPACTION/simulations/multipaction/${prefix}_
        done
    done
done
