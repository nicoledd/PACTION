for numModes in {3..6}; do
    for nsamples in {1..5}; do
        for seed in {1..3}; do
            prefix=${seed}'_'${nsamples}'_'${numModes};
            echo $prefix;
            PACTIONPREFIX='/scratch/data/nsdong2/projectPACTION/'

            TREES=$PACTIONPREFIX"simulations/multipaction/"$prefix"_tree_mode0.tsv"
            PROPS=$PACTIONPREFIX"simulations/multipaction/"$prefix"_props_mode0.out"
            for ((i=1; i<numModes; i++)); do
                TREES+=" None"
                PROPS+=" "$PACTIONPREFIX"simulations/multipaction/"$prefix"_props_mode"$i".out"
            done

            #echo $TREES
            #echo $PROPS

            python3 /scratch/data/nsdong2/projectPACTION/newpaction/src/paction.py \
                --seed ${seed} --samples ${nsamples} --noise 0 \
                -p $PROPS \
                -t $TREES \
                -o /scratch/data/nsdong2/projectPACTION/results/multipaction/${prefix}_
        done
    done
done
