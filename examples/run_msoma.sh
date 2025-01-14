#Inputs
bam="synthetic_data/chrM_seed1_ART.bam"
bed="synthetic_data/chrM.bed"
fasta="synthetic_data/chrM.fa"

#Params
minMQ=20
requireFlags=2
excludeFlags=3844
maxIndel=0
mismatchFrac=0.1
softclipFrac=0.5
seqLength=75
ntrim=5
minBQ=20
maxAltAllele=1
minDepth=10


msoma count \
    --bed $bed \
    --fasta $fasta \
    --min-MQ $minMQ \
    --require-flags $requireFlags \
    --exclude-flags $excludeFlags \
    --max-indel $maxIndel \
    --mismatch-frac $mismatchFrac \
    --softclip-frac $softclipFrac \
    --seq-length $seqLength \
    --ntrim $ntrim \
    --min-BQ $minBQ \
    --min-depth $minDepth \
    --max-alt-allele $maxAltAllele \
    --output "example.counts.gz" \
    $bam


msoma mle \
    --output "example.pvals" \
    --ab "example.ab" \
    --min-depth $minDepth \
    example.counts.gz
