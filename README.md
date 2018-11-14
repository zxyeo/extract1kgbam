# extract1kgbam
Subset 1kG BAM with list of genomic regions (human only)

# Run
input path of bams \
input target files containing chromosome:start-stop \
input working directory

./extract1kgBAM.py \
--bamlist {{bam_path_list}} \
--target {{targetregion_list}} \
--workdir {{output_path}}

{{targetregion_list}} examples: \
1:12345-23456 \
1:44444-44500 \
3:55555-60100