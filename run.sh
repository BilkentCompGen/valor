##########################################################################################################
############################### EDIT FROM HERE ###########################################################
# pools should be named as $bamfile$pool$suffix for $pool=1 to $count                                   ##
# path of bam files  -- it is advised to remove duplicate read prior to running this algorithm          ##
bamfile=/mnt/compgen/homes/marzie/InVersionDetection/mapping/NA12878pool                                ##
# total number of pools                                                                                 ##
count=288                                                                                               ##
suffix=.bam                                                                                             ##
                                                                                                        ##
# read size information (mean +- 3std)                                                                  ##
FRAG_MAX=1000                                                                                           ##
FRAG_MIN=200                                                                                            ##
                                                                                                        ##
# output directory -- make sure there is enough space on the directory                                  ##
outputDir="allReady"                                                                                    ##
# min coverage filtering for inferred clones                                                            ##
coverageThreshold=0.4                                                                                   ##
                                                                                                        ##
##########################################################################################################
############################### EDIT UP TO HERE ##########################################################

# make the temp and allReady directory
mkdir $outputDir
mkdir $outputDir/temp/
mkdir $outputDir/Clusters/
mkdir $outputDir/Pools/

# compile the java src files
# -- don't forget to update the config file
javac -cp src src/InversionDetection.java
javac -cp src src/InferClones.java
mv src/*.class J


# run for each pool
for pool in $(seq 1 $count)
do
        echo "pool"$pool
        # bamToBed
        bedtools bamtobed -bedpe -i $bamfile$pool$suffix > $outputDir/temp/pool$pool.bed
        #separate the reads according to signature
        if [ "$pool" = 1 ]
        then
                awk '{if ($1==$4 && $9=="-" && $10=="-") {print "chr" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $9 "\t" $10}}' $outputDir/temp/pool$pool.bed > $outputDir/mmReads.tsv
                awk '{if ($1==$4 && $9=="+" && $10=="+") {print "chr" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $9 "\t" $10}}' $outputDir/temp/pool$pool.bed > $outputDir/ppReads.tsv
        else
                awk '{if ($1==$4 && $9=="-" && $10=="-") {print "chr" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $9 "\t" $10}}' $outputDir/temp/pool$pool.bed >> $outputDir/mmReads.tsv
                awk '{if ($1==$4 && $9=="+" && $10=="+") {print "chr"  $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $9 "\t" $10}}' $outputDir/temp/pool$pool.bed >> $outputDir/ppReads.tsv
        fi
        awk -v MAXREAD=$FRAG_MAX -v MINREAD=$FRAG_MIN '{if ($1==$4 && $9=="+" && $10=="-" && ($6-$2)>=MINREAD && ($6-$2)<=MAXREAD) {print "chr" $1 "\t" $2 "\t" $6}}' $outputDir/temp/pool$pool.bed > $outputDir/temp/pool$pool.tsv
        #make inferred clones
        java -cp J/ InferClones $outputDir/temp/pool$pool.tsv $outputDir/temp/pool$pool.region.tsv $pool
        if [ "$pool" = 1 ]
        then
                awk '{print}' $outputDir/temp/pool$pool.region.tsv > $outputDir/allRegions.tsv
        else
                awk '{print}' $outputDir/temp/pool$pool.region.tsv >> $outputDir/allRegions.tsv
        fi
        #rm $outputDir/temp/pool$pool.tsv
        rm $outputDir/temp/pool$pool.region.tsv
        rm $outputDir/temp/pool$pool.bed
done

#################################################################################################
# it is advised to filter out regions with low coverage and then proceed
# if you don't want to just comment this part
echo "filtering low coverage clones"
for pool in $(seq 1 $count)
do
        # bamToBed
        bedtools bamtobed -i $bamfile$pool$suffix | awk '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4}' > $outputDir/temp/pool$pool.bed
        # regions of this pool
        awk -v pool=$pool '{if ($4==pool) print}' $outputDir/allRegions.tsv | awk '{print $1 "\t" $2 "\t"  $3}' > $outputDir/temp/allRegions$pool.tsv
        # get coverage and filter
        coverageBed -a $outputDir/temp/pool$pool.bed -b $outputDir/temp/allRegions$pool.tsv | awk -v pool=$pool -v threshold=$coverageThreshold '{if ($7 >= threshold) print $1 "\t" $2 "\t" $3 "\t" pool "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > $outputDir/temp/pool$pool.region.filtered.tsv
        if [ "$pool" = 1 ]
       then
                awk '{print}' $outputDir/temp/pool$pool.region.filtered.tsv > $outputDir/allRegions.filtered.tsv
        else
                awk '{print}' $outputDir/temp/pool$pool.region.filtered.tsv >> $outputDir/allRegions.filtered.tsv
        fi
        rm $outputDir/temp/allRegions$pool.tsv
        rm $outputDir/temp/pool$pool.region.filtered.tsv
        rm $outputDir/temp/pool$pool.region.tsv
        rm $outputDir/temp/pool$pool.bed
done
#mv $outputDir/allRegions.filtered.tsv $outputDir/allRegions.tsv
#################################################################################################

# get rid of temp directory
rmdir $outputDir/temp
# it is advised to filter out regions with low coverage and then proceed
java -cp J/ InversionDetection $outputDir/allRegions.tsv all $outputDir/ppReads.tsv $outputDir/mmReads.tsv $outputDir/Clusters/ $outputDir/Pools/

echo "It is advised to removed inversions laying on gaps and duplicated regions"
# remove overlaps by gaps and duplicates

