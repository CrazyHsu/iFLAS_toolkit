############ wangbo un-hybrid
python iflas.py preproc -c iFLAS_config_wangbo.cfg
python iflas.py mapping -c iFLAS_config_wangbo.cfg
python iflas.py filter -c iFLAS_config_wangbo.cfg
python iflas.py collapse -c iFLAS_config_wangbo.cfg
python iflas.py find_as -c iFLAS_config_wangbo.cfg
python iflas.py visual_as -c iFLAS_config_wangbo.cfg -g Zm00001d027231
python iflas.py rank_as -c iFLAS_config_wangbo.cfg
python iflas.py diff_as -c iFLAS_config_wangbo.cfg -d wangbo_compCond.txt
python iflas.py go -c iFLAS_config_wangbo.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3 
python iflas.py report -c iFLAS_config_wangbo.cfg

############ wanglab
python iflas.py preproc -c iFLAS_config_wanglab.cfg
python iflas.py mapping -c iFLAS_config_wanglab.cfg
python iflas.py filter -c iFLAS_config_wanglab.cfg
python iflas.py collapse -c iFLAS_config_wanglab.cfg
python iflas.py find_as -c iFLAS_config_wanglab.cfg
python iflas.py visual_as -c iFLAS_config_wanglab.cfg -g Zm00001d027231
python iflas.py rank_as -c iFLAS_config_wanglab.cfg
python iflas.py diff_as -c iFLAS_config_wanglab.cfg -d wanglab_compCond.txt
python iflas.py go -c iFLAS_config_wanglab.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3 
python iflas.py report -c iFLAS_config_wanglab.cfg

############ ont
python iflas.py preproc -c iFLAS_config_ont.cfg
python iflas.py mapping -c iFLAS_config_ont.cfg
python iflas.py filter -c iFLAS_config_ont.cfg
python iflas.py collapse -c iFLAS_config_ont.cfg
python iflas.py find_as -c iFLAS_config_ont.cfg
python iflas.py visual_as -c iFLAS_config_ont.cfg -g Zm00001d027231
python iflas.py rank_as -c iFLAS_config_ont.cfg
python iflas.py diff_as -c iFLAS_config_ont.cfg -d ont_compCond.txt
python iflas.py palen_as -c iFLAS_config_ont.cfg
python iflas.py go -c iFLAS_config_ont.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3 
python iflas.py report -c iFLAS_config_ont.cfg

############ wangbo hybrid
python iflas.py preproc -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py mapping -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py filter -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py collapse -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py find_as -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py visual_as -c iFLAS_config_wangbo_hybrid.cfg -g Zm00001d027231
python iflas.py rank_as -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py diff_as -c iFLAS_config_wangbo_hybrid.cfg -d wangbo_hybrid_compCond.txt
python iflas.py allelic_as -c iFLAS_config_wangbo_hybrid.cfg
python iflas.py go -c iFLAS_config_wangbo_hybrid.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2goId.txt -s sample1,sample2,sample3
python iflas.py report -c iFLAS_config_wangbo_hybrid.cfg