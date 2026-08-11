# Pre-F–K null production readiness

## Recommended first pilot

Run 20 globally distinct null realizations on the same first 30 files used by
the corrected December pilot. Each realization includes both surrogate
families and both signed branches. This is a computational pilot, not a final
significance test: 20 realizations limit the smallest empirical p value to
1/21, approximately 0.048.

```bash
mkdir -p /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_full_pipeline_null_v2_pilot_n30_r20

sbatch --array=0-19 \
  --export=ALL,NULL_DATE=2024-12-20,NULL_FILE_START=0,NULL_FILE_COUNT=30,NULL_INDEX_OFFSET=0,NULLS_PER_TASK=1,NULL_OUTPUT_DIR=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_full_pipeline_null_v2_pilot_n30_r20 \
  awd_clean/ambient_fk_full_pipeline_null_v2_array.sbatch
```

The two null methods remain at their launcher defaults. Avoid passing their
comma-separated value through `sbatch --export`, because Slurm also uses commas
to delimit exported variables.

After every array task completes, aggregate with an exact completeness check:

```bash
/home/users/nberrios/miniconda3/envs/das/bin/python \
  awd_clean/aggregate_ambient_fk_full_pipeline_null_v2_checked.py \
  --input-dir /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_full_pipeline_null_v2_pilot_n30_r20 \
  --expected-start 0 --expected-count 20
```

The checked aggregator rejects missing or duplicate realization IDs, partial
JSON/NPZ pairs, inconsistent geometry, different observed metrics, and
incompatible methods or seeds. It concatenates only null distributions. The
recomputed observed result is equality-checked and retained once, not counted
once per array task.

## Expected cost

The one-file smoke with one realization from both null families required about
9 seconds, including the observed calculation. Linear scaling predicts about
4.5 minutes per 30-file task. Allow 5–10 minutes for filesystem and cluster
variability. Twenty tasks represent about 100–200 node-minutes, or roughly
7–14 allocated CPU-hours at four CPUs per task. Queue time is separate. The
launcher requests four hours as a conservative ceiling.

Do not expand directly to full seasonal days until this pilot shows sensible
observed and surrogate distributions. If it does, the next defensible step is
a predeclared held-out day with at least 99 realizations, which gives empirical
p-value resolution of 0.01.
