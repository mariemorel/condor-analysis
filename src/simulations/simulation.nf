params.outdir="results_msa_simu"
params.scale=1

process RerootTree {
publishDir "${params.outdir}", mode: 'copy'

label 'gotree'

input:
path tree
path outgroup

output:
path "*_reroot_keep.treefile"
path "*_reroot.treefile"

script:
"""
gotree reroot outgroup -i ${tree} -l ${outgroup} | gotree unroot | gotree rename -a --tips=false --internal > ${tree.baseName}_reroot_keep.treefile
gotree reroot outgroup -i ${tree} -l ${outgroup} -r | gotree rename -a --tips=false --internal > ${tree.baseName}_reroot.treefile
"""
}

process MSARemoveOutgroup {

label 'goalign'

input:
path msa
path outgroup

output:
path "${msa.baseName}_nooutgroup.fasta"

script:
"""
goalign subset -i $msa -f $outgroup -r -o ${msa.baseName}_nooutgroup.fasta
"""
}

process ReplaceConvByGaps {
publishDir "${params.outdir}", mode: 'copy'

label 'goalign'

input:
path msa
path convlist

output:
path "${msa.baseName}_gaps.fasta"

script:
"""
grep -v "^#" $convlist | awk '{print \$1 "\t" \$2 "\t-"}' > convgaps.txt
goalign replace --posfile convgaps.txt -i $msa > ${msa.baseName}_gaps.fasta
"""
}

process AliLen {
label 'goalign'

input:
path msa

output:
stdout

script:
"""
printf \$(goalign stats length -i $msa)
"""
}

process PrepareAlignToPastML {
publishDir "${params.outdir}", mode: 'copy'

label 'python'

input:
path msa
val msalen
path matrices
path rates

output:
path "*pastml_input.tsv.gz"
path "parameter_*"
path "*pastml_matrix"

script:
"""
convert_align_to_pastml.py $msa ${msa.baseName}_

# Extract frequencies
freqs=`grep -i JTT -A 20 ${matrices} | tail -n 1 | sed 's/;//'`
AA="A R N D C Q E G H I L K M F P S T W Y V"
paste <(tr ' ' '\\n' <<< \${AA[*]}) <(tr ' ' '\\n' <<< \${freqs[*]}) >  frequencies.txt

# Extract matrix
matrix_pastml_format.py JTT $matrices

# PastML parameter files
for i in {1..$msalen}
do 
	R=`sed -n "\${i}"p ${rates}` 
	echo -e 'parameter\tvalue' > parameter_\$((\${i}-1))
	cat frequencies.txt >> parameter_\$((\${i}-1))
	echo -e "scaling_factor\t\${R}" >> parameter_\$((\${i}-1))
done
"""
}

process AncestralRootSequencePastML{
publishDir "${params.outdir}",mode:'copy'

label 'pastml'

input:
path msatsv
val msalen
path parameters
path matrix
path tree

output:
path "${tree.baseName}_root.fasta"

//rate sed numerotation from 1 ($line -1)
//input numerotation from 0
//create parameter file for pastml including rate per site (scaling factor) for each site and frequencies of amino acids
//run pastml 
//remove some temp files
//retrieve marginal proba for root (sed -n 2p)
script:
"""
ancestralPASTML.sh "$parameters" $matrix $msalen $msatsv $tree ${task.cpus} ${tree.baseName}
"""
}

process AncestralRootSequenceRAxML {
publishDir "${params.outdir}", mode: 'copy'

label 'raxml'

input:
path msa
path tree

output:
path "${tree.baseName}_root.fasta"

script:
"""
### Try to simulate with JTT Instead
# We try something else: reconstructing the ancestral sequence with raxml + simulating directly with outgroup
raxml-ng --ancestral --msa ${msa} --tree ${tree} --model JTT+R3
grep "N000000001" *.raxml.ancestralStates | awk '{print ">" \$1 "\\n" \$2}' > ${tree.baseName}_root.fasta
"""
}


process AddGapEnds {
publishDir "${params.outdir}", mode: 'copy'

label 'goalign'

input:
path root

output:
path "${root.baseName}_gaps.fasta"

script:
"""
# We add the "gaps" back to the root sequence
goalign replace -e -s '^' -n "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL" -i ${root} --unaligned \
	| goalign replace -e -s "\$" -n "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL" --unaligned > ${root.baseName}_gaps.fasta
"""
}

process EstimateRates {
publishDir "${params.outdir}", mode: 'copy'

label 'iqtree'

input:
path msa
path tree

output:
path "*_rates.txt"

script:
"""
# We estimate the rates per sites
iqtree -m JTT+R3 -nt ${task.cpus} -s $msa -te $tree -wsr
#iqtree -s $msa -t $tree -m JTT+R3 -n 0 --rate
grep -v "^#"  *.rate | tail -n+2| cut -f 2 | sed 's/1.00000/0.00000/g' > ${msa}_rates.txt
"""
}

process ReoptimizeTree {
publishDir "${params.outdir}", mode: 'copy'

label 'iqtree'

input:
path msa
path tree

output:
path "${msa}.treefile"

script:
"""
# We reestimate branch lengths
iqtree -m JTT+R3 -nt ${task.cpus} -s $msa -te $tree -wsr
"""
}


process SimulateMarie {

publishDir "${params.outdir}", mode: 'copy'

label 'python'

input:
path tree
path root
path rates
path matrices

output:
path "simulated_alignment_jtt.fasta"
path "simulated_alignment_ancestral.fasta"

script:
"""
echo \$(grep -v "^>" $root | tr -d '\\n') > root
freqs=`grep -i jtt -A 20 $matrices | tail -n 1 | sed 's/;//'`
AA="A R N D C Q E G H I L K M F P S T W Y V"
paste <(tr ' ' '\n' <<< \${AA[*]}) <(tr ' ' '\n' <<< \${freqs[*]}) >  frequencies.txt
matrix_pastml_format.py JTT $matrices
simple_simulator.py root $tree out $rates frequencies.txt JTTsimulator_matrix.model > simulated_alignment_jtt.fasta
"""
}

process CountEEMs{
publishDir "${params.outdir}", mode: 'copy'

label 'gotreedev'

input:
path tree
path align

output:
path 'simulated_alignment_eems.tsv'

script:
"""
gotree compute mutations --eems -i $tree -a $align > simulated_alignment_eems.tsv
"""
}

process ScaleTree {
publishDir "${params.outdir}", mode: 'copy'

label 'gotree'

input:
path tree
val scale

output:
path "${tree.baseName}_scale.treefile"

script:
"""
gotree brlen scale -f $scale -i $tree -o ${tree.baseName}_scale.treefile
"""
}

process AddConvergentMutations {
publishDir "${params.outdir}", mode: 'copy'

label 'goalign'

input:
path msa
path positions

output:
path "${msa.baseName}_convergent.fasta"

script:
"""
# Then we add the convergent mutations
goalign replace --posfile $positions -i $msa > tmp
# We set the GAPS as they were in the initial alignment
goalign mask -i tmp --start 0 --length 453 --replace GAP | goalign mask --start 908 --length 10000 --replace GAP > ${msa.baseName}_convergent.fasta
"""
}


workflow{
    tree=file(params.tree)
    msa=file(params.msa)
    outgroup=file(params.outgroup)
    convpos=file(params.convpos)
    matrices=file(params.matrices)
    treescale=params.scale
    outdir=file(params.outdir)

    msagapskeep=ReplaceConvByGaps(msa,convpos)
    msagaps = MSARemoveOutgroup(msagapskeep,outgroup)

    optim=ReoptimizeTree(msagapskeep,tree)
    rerootout=RerootTree(optim,outgroup)
    reroot = rerootout[1]
    rerootfull = rerootout[0]

    rates=EstimateRates(msagaps,reroot)
    
    msalen = AliLen(msagaps)
    pastmlprep=PrepareAlignToPastML(msagaps, msalen, matrices, rates)
    msapastml=pastmlprep[0]
    paramspastml=pastmlprep[1]
    matrixpastml=pastmlprep[2]
    root=AncestralRootSequencePastML(msapastml,msalen, paramspastml, matrixpastml, reroot)
    
    scaledtree = ScaleTree(rerootfull, treescale)
    simumsa=SimulateMarie(scaledtree, root, rates, matrices)
    simueems=CountEEMs(scaledtree, simumsa[1])
    convmsa=AddConvergentMutations(simumsa[0],convpos)
}
