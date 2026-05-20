# cd /projects/caeg/people/ngm902/apps/repos/rocs/results/microbial/sourcetracker
# mamba activate /home/ngm902/micromamba/envs/st2


#!/bin/bash
sourcetracker2 gibbs -i st_sp_median_wviruses.biom -m st_sp_median_wviruses.map  -o st_sp_median_wviruses --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments
# sourcetracker2 gibbs -i st_sp_median_wviruses.biom -m st_sp_median_wviruses.map  -o st_sp_median_wviruses --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments --restarts 100 --draws_per_restart 5 --diagnostics

# sourcetracker2 gibbs -i ....biom -m ....map -o ... --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments --restarts 100 --draws_per_restart 5 --diagnostics

# sourcetracker2 gibbs -i st-biome-class-gm_v3.biom -m st-biome-class-gm_v3.map -o st-biome-class-gm_v3 --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments
#sourcetracker2 gibbs -i st-biome-class-gm.biom -m st-biome-class-gm.map -o st-biome-class-gm2 --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments

#sourcetracker2 gibbs -i st-subclass-all.biom -m st-subclass-all.map -o st-subclass-all_basic2 --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments
#sourcetracker2 gibbs -i st-subclass-all_collapsed.biom -m st-subclass-all_collapsed.map -o st-subclass-all_basic_collapsed --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments

