#
# Step #1
# Process individual samples
#
export QuickIsoSeq=/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeqPro
export PATH=$QuickIsoSeq:$PATH

run-isoseq.sh allIDs.txt run.config


#
# Step #2
# Summarization   Note: run this step after Step #1 is done
#
export QuickIsoSeq=/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeqPro
export PATH=$QuickIsoSeq:$PATH

#make sure your R envorionment is available. Below is speicifc for Pfizer's environment
ml RHEL6-apps R/3.2.3   

merge-isoseq.sh allIDs.txt  run.config &> Results.log
#bsub -q short -app large "merge-isoseq.sh allIDs.txt  run.config &> Results.log"
