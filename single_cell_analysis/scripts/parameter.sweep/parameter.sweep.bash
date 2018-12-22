# Qiwen Hu 2018
# submit jobs for parameter sweep

#learning rate
learning_rate=(0.0005 0.001 0.002)

#batch size
batch_size=(50 100 200)

#epoch
epoch=(25 50 100 200)

#depth
depth=(2 3)

#first layer dimension
first_layer=(100 250 500)

#
for f in `cat data.list.txt`;
do
 for l in ${learning_rate[@]};
 do
  for b in ${batch_size[@]};
  do
   for e in ${epoch[@]};
    do
      for d in ${depth[@]};
      do
       for c in ${first_layer[@]};
        do
          bsub -q gpu -R "select[ngpus>0] rusage [ngpus_shared=2]" python vae_singlecell.v4.py -f $f -l $l -b $b -e $e -d $d -c $c
	 done
      done
    done
  done
 done
done
