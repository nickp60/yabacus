#!/bin/bash
set -uoe pipefail
IFS=$'\n\t'

#
FLANK=1000
if [ -z "$1" ]
then
    SEED=17
else
    SEED="$1"
fi
OUTDIR="./simulatedGenomeResults_${SEED}/"
mkdir ${OUTDIR}

##############################################################################
toyref="BA000007.2"
#ref="NC_000913.3"
# get the reference Sakai genome
if [ -e "${toyref}.fasta" ]
then
echo "using local copy of reference"
else
get_genomes.py -q ${toyref} -o ./
fi

# annotate regions
ribo scan ./${toyref}.fasta -o ${OUTDIR}toyGenome/scan/
# Cluster regions
ribo select ${OUTDIR}toyGenome/scan/scannedScaffolds.gb -o ${OUTDIR}toyGenome/select/
# extract regions with 5kb flanking
ribo snag ${OUTDIR}toyGenome/scan/scannedScaffolds.gb  ${OUTDIR}toyGenome/select/riboSelect_grouped_loci.txt -o ${OUTDIR}toyGenome/snag/ -l 5000 --just_extract

# only bother with 2 rDNAs, # 1 and 2

for rdna in {3..7}
do
    rm ${OUTDIR}toyGenome/snag/flanking_regions_output/*_region_${rdna}_riboSnag_flanking_regions.fasta
done

# combine the extracted regions into a toy genome
python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIR}toyGenome/snag/ \*_riboSnag.fasta -o ${OUTDIR}toyGenome/coli_genome/
# generate reads from the toy genome simulating a MiSeq V3 run
~/bin/pirs-2.0.2/pirs simulate -m 300 -l 100 -x 30 -v 10 -o ${OUTDIR}/toyGenome/reads -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz --compress -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat ${OUTDIR}/toyGenome/coli_genome/concatenated_seq.fasta -S ${SEED}


# annotate the toy genome for mauve visualization
ribo scan ${OUTDIR}toyGenome/coli_genome/concatenated_seq.fasta -o ${OUTDIR}toyGenome/coli_genome/scan/

##############################################################################

ref="NC_000913.3"
# get the genome
echo "Downloading reference: ${ref}"
mkdir ${OUTDIR}ref
if [ -e ${ref}.fasta ] ;
then
echo "using existing version of $ref"
else
get_genomes.py -q $ref -o ./
fi

# anotate the rDNAs

ribo scan ./${ref}.fasta -o ${OUTDIR}ref/scan/
# cluster
ribo select ${OUTDIR}ref/scan/scannedScaffolds.gb  -o ${OUTDIR}ref/select/
# run riboSeed
ribo seed -r ${OUTDIR}ref/scan/scannedScaffolds.gb  -o ${OUTDIR}ref/seed/ ${OUTDIR}ref/select/riboSelect_grouped_loci.txt -F ${OUTDIR}toyGenome/reads_100_300_1.fq.gz -R ${OUTDIR}toyGenome/reads_100_300_2.fq.gz -i 3 -c 2 -v 1 -l ${FLANK}
