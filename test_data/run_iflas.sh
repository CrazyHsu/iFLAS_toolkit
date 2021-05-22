############ wangbo un-hybrid
python iflas.py preproc -c test_data/iFLAS_config_wangbo.cfg
python iflas.py mapping -c test_data/iFLAS_config_wangbo.cfg
python iflas.py filter -c test_data/iFLAS_config_wangbo.cfg
python iflas.py collapse -c test_data/iFLAS_config_wangbo.cfg
python iflas.py find_as -c test_data/iFLAS_config_wangbo.cfg
python iflas.py visual_as -c test_data/iFLAS_config_wangbo.cfg -g Zm00001d050245
python iflas.py rank_as -c test_data/iFLAS_config_wangbo.cfg
python iflas.py diff_as -c test_data/iFLAS_config_wangbo.cfg -d test_data/wangbo_compCond.txt
python iflas.py go -c test_data/iFLAS_config_wangbo.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -c test_data/iFLAS_config_wangbo.cfg

############ wangbo un-hybrid merged
python iflas.py preproc -c test_data/iFLAS_config_wangbo_merge.cfg
python iflas.py mapping -c test_data/iFLAS_config_wangbo_merge.cfg
python iflas.py filter -c test_data/iFLAS_config_wangbo_merge.cfg
python iflas.py collapse -c test_data/iFLAS_config_wangbo_merge.cfg
python iflas.py find_as -c test_data/iFLAS_config_wangbo_merge.cfg
python iflas.py visual_as -c test_data/iFLAS_config_wangbo_merge.cfg -g Zm00001d050245
python iflas.py rank_as -c test_data/iFLAS_config_wangbo_merge.cfg
python iflas.py diff_as -c test_data/iFLAS_config_wangbo_merge.cfg -d test_data/wangbo_compCond.txt
python iflas.py go -c test_data/iFLAS_config_wangbo_merge.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -c test_data/iFLAS_config_wangbo_merge.cfg

############ wanglab
python iflas.py preproc -c test_data/iFLAS_config_wanglab.cfg
python iflas.py mapping -c test_data/iFLAS_config_wanglab.cfg
python iflas.py filter -c test_data/iFLAS_config_wanglab.cfg
python iflas.py collapse -c test_data/iFLAS_config_wanglab.cfg
python iflas.py find_as -c test_data/iFLAS_config_wanglab.cfg
python iflas.py visual_as -c test_data/iFLAS_config_wanglab.cfg -g Zm00001d050245
python iflas.py rank_as -c test_data/iFLAS_config_wanglab.cfg
python iflas.py diff_as -c test_data/iFLAS_config_wanglab.cfg -d test_data/wanglab_compCond.txt
python iflas.py go -c test_data/iFLAS_config_wanglab.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -c test_data/iFLAS_config_wanglab.cfg

############ ont
python iflas.py preproc -c test_data/iFLAS_config_ont.cfg
python iflas.py mapping -c test_data/iFLAS_config_ont.cfg
python iflas.py filter -c test_data/iFLAS_config_ont.cfg
python iflas.py collapse -c test_data/iFLAS_config_ont.cfg
python iflas.py find_as -c test_data/iFLAS_config_ont.cfg
python iflas.py visual_as -c test_data/iFLAS_config_ont.cfg -g Zm00001d050245
python iflas.py rank_as -c test_data/iFLAS_config_ont.cfg
python iflas.py diff_as -c test_data/iFLAS_config_ont.cfg -d test_data/ont_compCond.txt
python iflas.py palen_as -c test_data/iFLAS_config_ont.cfg
python iflas.py go -c test_data/iFLAS_config_ont.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -c test_data/iFLAS_config_ont.cfg

############ wangbo hybrid
python iflas.py preproc -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py mapping -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py filter -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py collapse -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py find_as -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py visual_as -c test_data/iFLAS_config_wangbo_hybrid.cfg -g Zm00001d050245
python iflas.py rank_as -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py diff_as -c test_data/iFLAS_config_wangbo_hybrid.cfg -d test_data/wangbo_hybrid_compCond.txt
python iflas.py allelic_as -c test_data/iFLAS_config_wangbo_hybrid.cfg
python iflas.py go -c test_data/iFLAS_config_wangbo_hybrid.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2goId.txt -s sample1,sample2,sample3
python iflas.py report -c test_data/iFLAS_config_wangbo_hybrid.cfg