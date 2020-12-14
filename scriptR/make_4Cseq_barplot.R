require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)

bless80 <- "/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed" %>% read_bed()


#BOXPLOT SUR LES PROFILES 4C-seq
replaceName <- list(
    "4C-Seq_C_Legube_1_fastq_after_trimming"="C_mOHT"           
    ,"4C-Seq_C_Legube_2_fastq_after_trimming"="C_pOHT"           
    ,"4C-Seq_C_Legube_3_fastq_after_trimming"="C_pOHTpDNAPKi"           
    ,"4C-Seq_C_Legube_4_fastq_after_trimming"="C_pOHTpATMi"           
    ,"4C-seq_E-Legube_1_S1_all_R1_001_cutadapt"="E_mOHT"         
    ,"4C-seq_E-Legube_2_S2_all_R1_001_cutadapt"="E_pOHT24h"         
    ,"4C-Seq-D-Legube-Banque1_S1_all_R1_001_cutadapt"="D_mOHT"   
    ,"4C-Seq-D-Legube-Banque2_S2_all_R1_001_cutadapt"="D_pOHT"   
    ,"4C-Seq-D-Legube-Banque3_S3_all_R1_001_cutadapt"="D_pOHTpDNAPKi"  
    ,"4C-Seq-D-Legube-Banque4_S4_all_R1_001_cutadapt"="D_pOHTpATMi"   
    ,"4C-Seq-Dbis-Legube-Banque1_S1_all_R1_001_cutadapt"="Dbis_mOHT"
    ,"4C-Seq-Dbis-Legube-Banque2_S2_all_R1_001_cutadapt"="Dbis_pOHT"
    ,"4C-Seq-Dbis-Legube-Banque3_S3_all_R1_001_cutadapt"="Dbis_pOHTpDNAPKi"
    ,"4C-Seq-Dbis-Legube-Banque4_S4_all_R1_001_cutadapt"="Dbis_pOHTpATMi"
    ,"4C-seq-F-Legube_1_fastq_after_trimming"  ="F_G1mOHT"         
    ,"4C-seq-F-Legube_2_fastq_after_trimming" ="F_G1pOHT"          
    ,"4C-seq-F-Legube_3_fastq_after_trimming" ="F_G2mOHT"          
    ,"4C-seq-F-Legube_4_fastq_after_trimming" ="F_G2pOHT"          
    ,"4C-seq-F-Legube_5_fastq_after_trimming" ="F_SmOHT"          
    ,"4C-seq-F-Legube_6_fastq_after_trimming"  ="F_SpOHT"         
    ,"4Cseq-LBCMCP-index1_S1_all_R1_001_cutadapt"="A_mOHT"       
    ,"4Cseq-LBCMCP-index2_S2_all_R1_001_cutadapt"="A_pOHT"       
    ,"LEGU-7_1_fastq_after_trimming"="G_siCTRLmOHT"                    
    ,"LEGU-7_2_fastq_after_trimming" ="G_siCTRLpOHT"                   
    ,"LEGU-7_3_fastq_after_trimming"="G_siSCC1mOHT"                    
    ,"LEGU-7_4_fastq_after_trimming"="G_siSCC1pOHT"                     
    ,"LEGU-7_5_fastq_after_trimming"="G_siCTCFmOHT"                    
    ,"LEGU-7_6_fastq_after_trimming"="G_siCTCFpOHT"                    
    ,"LEGU-8_1_fastq_after_trimming"="H_siCTRLmOHT"                    
    ,"LEGU-8_2_fastq_after_trimming"="H_siCTRLpOHT"                    
    ,"LEGU-8_3_fastq_after_trimming"="H_siSCC1mOHT"                    
    ,"LEGU-8_4_fastq_after_trimming"="H_siSCC1pOHT"                     
    ,"LEGU-8_5_fastq_after_trimming"="H_siCTCFmOHT"                    
    ,"LEGU-8_6_fastq_after_trimming"="H_siCTCFpOHT" 
)

#Using bamcompare

replaceName.compare <- list(
    "4C-Seq_A_pvsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr17_57578597-57579677.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr1_80484021-80485207.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr1_88985251-88986220.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr1_89659070-89659811.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr1_90042784-90045873.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr21_32863772-32864618.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_A_pvsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
    ,"4C-Seq_C_pATMivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
    ,"4C-Seq_C_pATMivsmOHT_Legube_chr17_57168614-57169531_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
    ,"4C-Seq_C_pATMivsmOHT_Legube_chr1_89455867-89456712_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
    ,"4C-Seq_C_pATMivsmOHT_Legube_chr21_33251469-33252587_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
    ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
    ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr17_57168614-57169531_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
    ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr1_89455867-89456712_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
    ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr21_33251469-33252587_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
    ,"4C-Seq_C_pvsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_C_pvsmOHT"
    ,"4C-Seq_C_pvsmOHT_Legube_chr17_57168614-57169531_pvsmOHT"="4C-Seq_C_pvsmOHT"
    ,"4C-Seq_C_pvsmOHT_Legube_chr1_89455867-89456712_pvsmOHT"="4C-Seq_C_pvsmOHT"
    ,"4C-Seq_C_pvsmOHT_Legube_chr21_33251469-33252587_pvsmOHT"="4C-Seq_C_pvsmOHT"
    ,"4C-Seq_Dbis_pATMivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_Dbis_pATMivsmOHT"
    ,"4C-Seq_Dbis_pATMivsmOHT_Legube_chr20_30946314-30947710.bw_pvsmOHT"="4C-Seq_Dbis_pATMivsmOHT"
    ,"4C-Seq_Dbis_pDNAPKivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_Dbis_pDNAPKivsmOHT"
    ,"4C-Seq_Dbis_pDNAPKivsmOHT_Legube_chr20_30946314-30947710.bw_pvsmOHT"="4C-Seq_Dbis_pDNAPKivsmOHT"
    ,"4C-Seq_Dbis_pvsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_Dbis_pvsmOHT"
    ,"4C-Seq_Dbis_pvsmOHT_Legube_chr20_30946314-30947710.bw_pvsmOHT"="4C-Seq_Dbis_pvsmOHT"
    ,"4C-Seq_D_pATMivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
    ,"4C-Seq_D_pATMivsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
    ,"4C-Seq_D_pATMivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
    ,"4C-Seq_D_pATMivsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
    ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
    ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
    ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
    ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
    ,"4C-Seq_D_pvsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
    ,"4C-Seq_D_pvsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
    ,"4C-Seq_D_pvsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
    ,"4C-Seq_D_pvsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
    ,"4C-seq_E-Legube_p24vsmOHT_chr17_45759770-45760603_pvsmOHT"="4C-seq_E_p24vsmOHT"
    ,"4C-seq_E-Legube_p24vsmOHT_chr17_57168614-57169531_pvsmOHT"="4C-seq_E_p24vsmOHT"
    ,"4C-seq_E-Legube_p24vsmOHT_chr1_89455867-89456712_pvsmOHT"="4C-seq_E_p24vsmOHT"
    ,"4C-seq_E-Legube_p24vsmOHT_chr1_90292860-90295480_pvsmOHT"="4C-seq_E_p24vsmOHT"
    ,"4C-seq_E-Legube_p24vsmOHT_chr20_30946314-30947710_pvsmOHT"="4C-seq_E_p24vsmOHT"
    ,"4C-seq_E-Legube_p24vsmOHT_chr21_33251469-33252587_pvsmOHT"="4C-seq_E_p24vsmOHT"
    ,"4C-seq-F-Legube_pvsmOHT_G1_chr17_45759770-45760603_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
    ,"4C-seq-F-Legube_pvsmOHT_G1_chr17_57168614-57169531_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
    ,"4C-seq-F-Legube_pvsmOHT_G1_chr1_89455867-89456712_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
    ,"4C-seq-F-Legube_pvsmOHT_G1_chr20_30946314-30947710_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
    ,"4C-seq-F-Legube_pvsmOHT_G1_chr21_33251469-33252587_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
    ,"4C-seq-F-Legube_pvsmOHT_G2_chr17_45759770-45760603_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
    ,"4C-seq-F-Legube_pvsmOHT_G2_chr17_57168614-57169531_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
    ,"4C-seq-F-Legube_pvsmOHT_G2_chr1_89455867-89456712_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
    ,"4C-seq-F-Legube_pvsmOHT_G2_chr20_30946314-30947710_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
    ,"4C-seq-F-Legube_pvsmOHT_G2_chr21_33251469-33252587_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
    ,"4C-seq-F-Legube_pvsmOHT_S_chr17_45759770-45760603_pvsmOHT"="4C-seq_F_pvsmOHT_S"
    ,"4C-seq-F-Legube_pvsmOHT_S_chr17_57168614-57169531_pvsmOHT"="4C-seq_F_pvsmOHT_S"
    ,"4C-seq-F-Legube_pvsmOHT_S_chr1_89455867-89456712_pvsmOHT"="4C-seq_F_pvsmOHT_S"
    ,"4C-seq-F-Legube_pvsmOHT_S_chr20_30946314-30947710_pvsmOHT"="4C-seq_F_pvsmOHT_S"
    ,"4C-seq-F-Legube_pvsmOHT_S_chr21_33251469-33252587_pvsmOHT"="4C-seq_F_pvsmOHT_S"
    ,"LEGU-G_pvsmOHT_siCTCF_chr17_45759770-45760603_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
    ,"LEGU-G_pvsmOHT_siCTCF_chr17_57168614-57169531_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
    ,"LEGU-G_pvsmOHT_siCTCF_chr1_89455867-89456712_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
    ,"LEGU-G_pvsmOHT_siCTCF_chr20_30946314-30947710_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
    ,"LEGU-G_pvsmOHT_siCTCF_chr21_33251469-33252587_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
    ,"LEGU-G_pvsmOHT_siCTRL_chr17_45759770-45760603_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
    ,"LEGU-G_pvsmOHT_siCTRL_chr17_57168614-57169531_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
    ,"LEGU-G_pvsmOHT_siCTRL_chr1_89455867-89456712_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
    ,"LEGU-G_pvsmOHT_siCTRL_chr20_30946314-30947710_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
    ,"LEGU-G_pvsmOHT_siCTRL_chr21_33251469-33252587_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
    ,"LEGU-G_pvsmOHT_siSCC1_chr17_45759770-45760603_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
    ,"LEGU-G_pvsmOHT_siSCC1_chr17_57168614-57169531_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
    ,"LEGU-G_pvsmOHT_siSCC1_chr1_89455867-89456712_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
    ,"LEGU-G_pvsmOHT_siSCC1_chr20_30946314-30947710_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
    ,"LEGU-G_pvsmOHT_siSCC1_chr21_33251469-33252587_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
    ,"LEGU-H_pvsmOHT_siCTCF_chr17_45759770-45760603_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
    ,"LEGU-H_pvsmOHT_siCTCF_chr17_57168614-57169531_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
    ,"LEGU-H_pvsmOHT_siCTCF_chr1_89455867-89456712_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
    ,"LEGU-H_pvsmOHT_siCTCF_chr20_30946314-30947710_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
    ,"LEGU-H_pvsmOHT_siCTCF_chr21_33251469-33252587_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
    ,"LEGU-H_pvsmOHT_siCTRL_chr17_45759770-45760603_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
    ,"LEGU-H_pvsmOHT_siCTRL_chr17_57168614-57169531_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
    ,"LEGU-H_pvsmOHT_siCTRL_chr1_89455867-89456712_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
    ,"LEGU-H_pvsmOHT_siCTRL_chr20_30946314-30947710_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
    ,"LEGU-H_pvsmOHT_siCTRL_chr21_33251469-33252587_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
    ,"LEGU-H_pvsmOHT_siSCC1_chr17_45759770-45760603_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
    ,"LEGU-H_pvsmOHT_siSCC1_chr17_57168614-57169531_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
    ,"LEGU-H_pvsmOHT_siSCC1_chr1_89455867-89456712_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
    ,"LEGU-H_pvsmOHT_siSCC1_chr20_30946314-30947710_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
    ,"LEGU-H_pvsmOHT_siSCC1_chr21_33251469-33252587_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
    ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
    ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
    ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
    ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
    ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
    ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
    ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
    ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
    ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
    ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
    ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
    ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
    ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
    ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
    ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
    ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
    ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
    ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
    ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
    ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
    #Derniers plots pour les dernières données (31 oct 19)
    
    ,"4C-Seq-I_pvsmOHT_G1_chr17_45759770-45760603_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
    ,"4C-Seq-I_pvsmOHT_G1_chr17_57168614-57169531_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
    ,"4C-Seq-I_pvsmOHT_G1_chr1_89455867-89456712_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
    ,"4C-Seq-I_pvsmOHT_G1_chr20_30946314-30947710_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
    ,"4C-Seq-I_pvsmOHT_G1_chr21_33251469-33252587_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
    ,"4C-Seq-I_pvsmOHT_G2_chr17_45759770-45760603_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
    ,"4C-Seq-I_pvsmOHT_G2_chr17_57168614-57169531_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
    ,"4C-Seq-I_pvsmOHT_G2_chr1_89455867-89456712_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
    ,"4C-Seq-I_pvsmOHT_G2_chr20_30946314-30947710_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
    ,"4C-Seq-I_pvsmOHT_G2_chr21_33251469-33252587_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
    ,"4C-Seq-I_pvsmOHT_S_chr17_45759770-45760603_pvsmOHT"="4C-seq_I_pvsmOHT_S"
    ,"4C-Seq-I_pvsmOHT_S_chr17_57168614-57169531_pvsmOHT"="4C-seq_I_pvsmOHT_S"
    ,"4C-Seq-I_pvsmOHT_S_chr1_89455867-89456712_pvsmOHT"="4C-seq_I_pvsmOHT_S"
    ,"4C-Seq-I_pvsmOHT_S_chr20_30946314-30947710_pvsmOHT"="4C-seq_I_pvsmOHT_S"
    ,"4C-Seq-I_pvsmOHT_S_chr21_33251469-33252587_pvsmOHT"="4C-seq_I_pvsmOHT_S"
    ,"4C-Seq-K_pvsmOHT_siCTCF_chr17_45759770-45760603_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
    ,"4C-Seq-K_pvsmOHT_siCTCF_chr17_57168614-57169531_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
    ,"4C-Seq-K_pvsmOHT_siCTCF_chr1_89455867-89456712_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
    ,"4C-Seq-K_pvsmOHT_siCTCF_chr20_30946314-30947710_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
    ,"4C-Seq-K_pvsmOHT_siCTCF_chr21_33251469-33252587_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
    ,"4C-Seq-K_pvsmOHT_siCTRL_chr17_45759770-45760603_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
    ,"4C-Seq-K_pvsmOHT_siCTRL_chr17_57168614-57169531_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
    ,"4C-Seq-K_pvsmOHT_siCTRL_chr1_89455867-89456712_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
    ,"4C-Seq-K_pvsmOHT_siCTRL_chr20_30946314-30947710_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
    ,"4C-Seq-K_pvsmOHT_siCTRL_chr21_33251469-33252587_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
    ,"4C-Seq-K_pvsmOHT_siSCC1_chr17_45759770-45760603_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
    ,"4C-Seq-K_pvsmOHT_siSCC1_chr17_57168614-57169531_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
    ,"4C-Seq-K_pvsmOHT_siSCC1_chr1_89455867-89456712_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
    ,"4C-Seq-K_pvsmOHT_siSCC1_chr20_30946314-30947710_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
    ,"4C-Seq-K_pvsmOHT_siSCC1_chr21_33251469-33252587_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
)

path4C.1 <- "/home/rochevin/Documents/PROJET_INGE/4CSeq_PROCESS_VINCENT/AFTER_EVA_PROCESSING/4C_ALN_TRANS/BAMCOMPARE"
path4C.2 <- "/home/rochevin/Documents/PROJET_INGE/4CSeq_PROCESS_VINCENT/AFTER_EVA_PROCESSING/4C_ALN_TRANS/BAMCOMPARE_GH"

wigs.4C <- c(list.files(path4C.1,pattern="_normalized.bw",full.names = T),list.files(path4C.2,pattern="_normalized.bw",full.names = T))

to.plot.4C.compare <- lapply(c(50000,100000,500000,1000000),function(bin){
    message(bin)
    lapply(wigs.4C,function(wig){
        one.w <- import.bw(wig,as="RleList")
        wig <- basename(wig)
        File <- replaceName.compare[[str_remove(wig,"_normalized.bw")]]
        vp <- str_extract(wig,"chr[0-9A-Z]+_[0-9]+-[0-9]+") %>% str_replace("_",":")
        x <- GRanges(vp)%>% anchor_center() %>% mutate(width = bin) 
        x$name <- vp
        Get1val(File,one.w,x) %>% mutate(viewpoint = vp) %>% mutate(binsize = bin)
    }) %>% bind_rows()
}) %>% bind_rows()
to.plot.4C.compare <- to.plot.4C.compare %>%
    mutate(Manip = str_extract(wig,"4C-[sS]eq_[A-Za-z0-9]+_")) %>%
    mutate(Condition = str_remove(wig,"4C-[sS]eq_[A-Za-z0-9]+_"))

bin=1000000
#PLOT1
subdf <- to.plot.4C.compare %>% filter(Condition %in% c("pvsmOHT","pATMivsmOHT","pDNAPKivsmOHT")) %>% mutate(Condition = fct_relevel(as.factor(Condition),c("pvsmOHT","pATMivsmOHT","pDNAPKivsmOHT")))

subsubdf <- subdf  %>% filter(binsize==bin,
                              viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),
                              Manip != "4C-Seq_A_")
subsubdf <- subsubdf %>% mutate(Group = ifelse(viewpoint == "chr17:45759770-45760603","CTRL","VIEWPOINT")) %>% filter(Group != "CTRL")
p1 <- subsubdf %>%
    ggplot(aes(x=Condition,y=value,fill = Condition)) +  
    geom_boxplot() + 
    facet_wrap(~Group,ncol=1) +
    theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + ggtitle(bin)
print(p1)
#PLOT2
subdf <- to.plot.4C.compare %>%
    filter(Condition %in% c("pvsmOHT_siCTRL","pvsmOHT_siSCC1","pvsmOHT_siCTCF")) %>%
    mutate(Condition = as.factor(Condition)) %>%
    mutate(Condition = fct_relevel(Condition,c("pvsmOHT_siCTRL","pvsmOHT_siSCC1","pvsmOHT_siCTCF")))
subsubdf <- subdf  %>% filter(binsize==bin) %>% filter(Manip != "4C-seq_K_")
subsubdf <- subsubdf %>% mutate(Group = ifelse(viewpoint == "chr17:45759770-45760603","CTRL","VIEWPOINT")) %>% filter(Group != "CTRL")
p2 <- subsubdf %>% 
    ggplot(aes(x=Condition,y=value,fill = Condition)) +  
    geom_boxplot() + 
    facet_wrap(~Group,ncol=1) +
    theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + ggtitle(bin)
print(p2)

#MANIP A/C/D/Dbis mean barplot/boxplot
subdf <- to.plot.4C.compare  %>% filter(binsize==bin,
                                        Manip %in% c("4C-Seq_A_","4C-Seq_C_","4C-Seq_D_","4C-Seq_Dbis_"),
                                        viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),
                                        Condition =="pvsmOHT")
subdf <- subdf %>% mutate(Group = ifelse(viewpoint == "chr17:45759770-45760603","CTRL","VIEWPOINT"))
p3 <- subdf %>%
    ggplot(aes(x=Condition,y=value,fill = Condition)) +  
    geom_boxplot() + 
    facet_wrap(~Group) +
    theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + ggtitle(bin)
print(p3)

p4 <- subdf %>%
    ggplot(aes(x=Group,y=value,fill = Group)) +  
    stat_summary(fun.y = mean, geom = "bar") + 
    geom_point() +
    stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2) +
    theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + ggtitle(bin)
print(p4)
