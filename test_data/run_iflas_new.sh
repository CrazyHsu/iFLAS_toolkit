############ wangbo un-hybrid
python iflas.py preproc -cfg test_data/iFLAS_config_wangbo_new.cfg
python iflas.py mapping -cfg test_data/iFLAS_config_wangbo_new.cfg
python iflas.py collapse -cfg test_data/iFLAS_config_wangbo_new.cfg
python iflas.py refine -cfg test_data/iFLAS_config_wangbo_new.cfg
python iflas.py find_as -cfg test_data/iFLAS_config_wangbo_new.cfg
python iflas.py visual_as -cfg test_data/iFLAS_config_wangbo_new.cfg -g Zm00001d050245
python iflas.py rank_as -cfg test_data/iFLAS_config_wangbo_new.cfg
python iflas.py diff_as -cfg test_data/iFLAS_config_wangbo_new.cfg -d test_data/wangbo_compCond.txt
python iflas.py go -cfg test_data/iFLAS_config_wangbo_new.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -cfg test_data/iFLAS_config_wangbo_new.cfg

############ wangbo un-hybrid merged
python iflas.py preproc -cfg test_data/iFLAS_config_wangbo_merge_new.cfg
python iflas.py mapping -cfg test_data/iFLAS_config_wangbo_merge_new.cfg
python iflas.py collapse -cfg test_data/iFLAS_config_wangbo_merge_new.cfg
python iflas.py refine -cfg test_data/iFLAS_config_wangbo_merge_new.cfg
python iflas.py find_as -cfg test_data/iFLAS_config_wangbo_merge_new.cfg
python iflas.py visual_as -cfg test_data/iFLAS_config_wangbo_merge_new.cfg -g Zm00001d050245
python iflas.py rank_as -cfg test_data/iFLAS_config_wangbo_merge_new.cfg
python iflas.py diff_as -cfg test_data/iFLAS_config_wangbo_merge_new.cfg -d test_data/wangbo_compCond.txt
python iflas.py go -cfg test_data/iFLAS_config_wangbo_merge_new.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -cfg test_data/iFLAS_config_wangbo_merge_new.cfg

############ wanglab
python iflas.py preproc -cfg test_data/iFLAS_config_wanglab_new.cfg
python iflas.py mapping -cfg test_data/iFLAS_config_wanglab_new.cfg
python iflas.py collapse -cfg test_data/iFLAS_config_wanglab_new.cfg
python iflas.py refine -cfg test_data/iFLAS_config_wanglab_new.cfg
python iflas.py find_as -cfg test_data/iFLAS_config_wanglab_new.cfg
python iflas.py visual_as -cfg test_data/iFLAS_config_wanglab_new.cfg -g Zm00001d050245
python iflas.py rank_as -cfg test_data/iFLAS_config_wanglab_new.cfg
python iflas.py diff_as -cfg test_data/iFLAS_config_wanglab_new.cfg -d test_data/wanglab_compCond.txt
python iflas.py go -cfg test_data/iFLAS_config_wanglab_new.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -cfg test_data/iFLAS_config_wanglab_new.cfg

############ ont
python iflas.py preproc -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py mapping -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py collapse -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py refine -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py find_as -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py visual_as -cfg test_data/iFLAS_config_ont_new.cfg -g Zm00001d050245
python iflas.py rank_as -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py diff_as -cfg test_data/iFLAS_config_ont_new.cfg -d test_data/ont_compCond.txt
python iflas.py palen_as -cfg test_data/iFLAS_config_ont_new.cfg
python iflas.py go -cfg test_data/iFLAS_config_ont_new.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2go.txt -s sample1,sample2,sample3
python iflas.py report -cfg test_data/iFLAS_config_ont_new.cfg

############ wangbo hybrid
python iflas.py preproc -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py mapping -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py collapse -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py refine -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py find_as -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py visual_as -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg -g Zm00001d050245
python iflas.py rank_as -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py diff_as -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg -d test_data/wangbo_hybrid_compCond.txt
python iflas.py allelic_as -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg
python iflas.py go -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg gene2goId.txt -s sample1,sample2,sample3
python iflas.py report -cfg test_data/iFLAS_config_wangbo_hybrid_new.cfg