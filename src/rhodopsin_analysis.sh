goalign subset -f spp_tokeep.txt -i final_alignment.translated.fullrhodopsin.fasta > to_keep.fullrhodo.fasta
goalign subset -f spp_toremove.txt -r -i to_keep.fullrhodo.fasta > rhodo_woutgroup.fasta
goalign clean sites -c 0.9 -i rhodo_woutgroup.fasta > clean_rhodo_woutgroup.fasta

iqtree -m MFP -s rhodo_woutgroup.fasta
