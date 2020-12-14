from snakemake.utils import R
THREADS = 3
RAW_PATH = "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/4C-seq/4CSeq_PROCESS_VINCENT/RAW/"
DEMULTIPLEX_PATH = "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/4C-seq/4CSeq_PROCESS_VINCENT/DEMULTIPLEXING/"
OUT_PATH = "../4C_ALN_TRANS"
EXTENSION = ".fastq.gz"
GENOME="/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/GENOMES/female.hg19/female.hg19.fa"

CONDITIONS = {
	"mvspOHT":{
		"sample_name":["4C-Seq_C_Legube_1_fastq_after_trimming",
		"4C-Seq_C_Legube_2_fastq_after_trimming",
		"4C-Seq-D-Legube-Banque1_S1_all_R1_001_cutadapt",
		"4C-Seq-D-Legube-Banque2_S2_all_R1_001_cutadapt"],
		"condition":["mOHT_C","pOHT_C","mOHT_D","pOHT_D"]
	}
}
#SCRIPTS
demultiplex_script="demultiplex.py"
make_frags="make_frags.R"
create_matrix_count="create_matrix_count.R"
makeDE="DE_processing_snakemake.R"
generate_bedGraph_smoothed="smoothData_adapted.R"
generate_bedGraph_norm="normFrags.R"
#PROCESS CONFIG
PRIMERS_PATH="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/4C-seq/4CSeq_PROCESS_VINCENT/RAW/PRIMERS/"
PRIMERS_FILES=[
"primer_file_demultiplexing_mapping_A.FA",
"primer_file_demultiplexing_mapping_C.FA",
"primer_file_demultiplexing_mapping_D.FA",
"primer_file_demultiplexing_mapping_E.FA",
"primer_file_demultiplexing_mapping_Dbis.FA",
"primer_file_demultiplexing_mapping_LEGU_7.FA",
"primer_file_demultiplexing_mapping_LEGU_8.FA",
"primer_file_demultiplexing_mapping_F.FA",
"primer_file_demultiplexing_mapping_I.FA",
"primer_file_demultiplexing_mapping_J.FA",
"primer_file_demultiplexing_mapping_K.FA",
"primer_file_demultiplexing_mapping_L.FA",
"primer_file_demultiplexing_mapping_M.FA"]
PRIMERS = {}
for f in PRIMERS_FILES:
	with open(PRIMERS_PATH+f) as file:
		PRIMERS[f] = [line.strip().replace(">","") for line in file if line.startswith(">")]
EXCLUDES = {}
EXCLUDES_PATH="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/4C-seq/4CSeq_PROCESS_VINCENT/RAW/EXCLUDES/"
for f in PRIMERS_FILES:
	with open(EXCLUDES_PATH+f) as file:
		EXCLUDES[f] = [line.strip() for line in file]


INFO_FASTQ = {
# "4C-Seq_C_Legube_1_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_C.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_C.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_C.FA"]
# }
# "4C-Seq_C_Legube_2_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_C.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_C.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_C.FA"]	
# }
# "4C-Seq_C_Legube_3_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_C.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_C.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_C.FA"]	
# }
# "4C-Seq_C_Legube_4_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_C.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_C.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_C.FA"]	
# }
# "4C-Seq-Dbis-Legube-Banque1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_C.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_C.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_C.FA"]	
# }
# "4C-Seq-Dbis-Legube-Banque2_S2_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]	
# }
# "4C-Seq-Dbis-Legube-Banque3_S3_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]		
# }
# "4C-Seq-Dbis-Legube-Banque4_S4_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]		
# }
# "4C-Seq-D-Legube-Banque1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]		
# }
# "4C-Seq-D-Legube-Banque2_S2_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]		
# }
# "4C-Seq-D-Legube-Banque3_S3_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]		
# }
# "4C-Seq-D-Legube-Banque4_S4_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_D.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_D.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_D.FA"]		
# }
# "4C-seq_E-Legube_1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_E.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_E.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_E.FA"]		
# }
# "4C-seq_E-Legube_2_S2_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_E.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_E.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_E.FA"]		
# }
# "4C-seq-F-Legube_1_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_F.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_F.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_F.FA"]	
# }
# "4C-seq-F-Legube_2_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_F.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_F.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_F.FA"]		
# }
# "4C-seq-F-Legube_3_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_F.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_F.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_F.FA"]		
# }
# "4C-seq-F-Legube_4_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_F.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_F.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_F.FA"]		
# }
# "4C-seq-F-Legube_5_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_F.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_F.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_F.FA"]		
# }
# "4C-seq-F-Legube_6_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_F.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_F.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_F.FA"]		
# }
# "4Cseq-LBCMCP-index1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_A.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_A.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_A.FA"]		
# }
# "4Cseq-LBCMCP-index2_S2_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_A.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_A.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_A.FA"]		
# }
# "LEGU-7_1_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_G.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_G.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_G.FA"]	
# }
# "LEGU-7_2_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_G.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_G.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_G.FA"]		
# }
# "LEGU-7_3_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_G.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_G.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_G.FA"]		
# }
# "LEGU-7_4_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_G.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_G.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_G.FA"]		
# }
# "LEGU-7_5_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_G.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_G.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_G.FA"]		
# }
# "LEGU-7_6_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_G.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_G.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_G.FA"]		
# }
# "LEGU-8_1_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_H.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_H.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_H.FA"]		
# }
# "LEGU-8_2_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_H.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_H.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_H.FA"]		
# }
# "LEGU-8_3_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_H.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_H.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_H.FA"]	
# }
# "LEGU-8_4_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_H.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_H.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_H.FA"]	
# }
# "LEGU-8_5_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_H.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_H.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_H.FA"]	
# }
# "LEGU-8_6_fastq_after_trimming":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_H.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_H.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_H.FA"]	
# },

# "4C-Seq-K_6_S6_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# 	},
# "4C-Seq-K_5_S5_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]		
# 	},
# "4C-Seq-K_4_S4_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]		
# 	},
# "4C-Seq-K_3_S3_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]		
# 	},
# "4C-Seq-K_2_S2_all_R1_001_cutadapt":{
			# "primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
			# "viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
			# "to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]		
# 	},
# "4C-Seq-K_1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]		
# 	},
# "4C_seq_J_index4_S4_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_J.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_J.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_J.FA"]		
# 	},
# "4C_seq_J_index3_S3_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_J.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_J.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_J.FA"]		
# 	},
# "4C_seq_J_index2_S2_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_J.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_J.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_J.FA"]		
# 	},
# "4C_seq_J_index1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_J.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_J.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_J.FA"]		
# 	},
# "4C_seq_I_index6_S6_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_I.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_I.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_I.FA"]		
# 	},
# "4C_seq_I_index5_S5_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_I.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_I.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_I.FA"]			
# 	},
# "4C_seq_I_index4_S4_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_I.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_I.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_I.FA"]			
# 	},
# "4C_seq_I_index3_S3_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_I.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_I.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_I.FA"]			
# 	},
# "4C_seq_I_index2_S2_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_I.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_I.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_I.FA"]			
# 	},
# "4C_seq_I_index1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_I.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_I.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_I.FA"]			
# 	},
# 
# "4C-Seq-L_4_S4_all_R1_001_cutadapt":{
			# "primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_L.FA",
			# "viewpoint":PRIMERS["primer_file_demultiplexing_mapping_L.FA"],
			# "to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_L.FA"]			
	# },
# "4C-Seq-L_3_S3_all_R1_001_cutadapt":{
			# "primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_L.FA",
			# "viewpoint":PRIMERS["primer_file_demultiplexing_mapping_L.FA"],
			# "to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_L.FA"]			
	# },
# "4C-Seq-L_2_S2_all_R1_001_cutadapt":{
			# "primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_L.FA",
			# "viewpoint":PRIMERS["primer_file_demultiplexing_mapping_L.FA"],
			# "to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_L.FA"]			
	# },
# "4C-Seq-L_1_S1_all_R1_001_cutadapt":{
			# "primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_L.FA",
			# "viewpoint":PRIMERS["primer_file_demultiplexing_mapping_L.FA"],
			# "to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_L.FA"]			
	# },


# "4Cseq_K_newseq_6_S9_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# },
# "4Cseq_K_newseq_5_S8_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# },
# "4Cseq_K_newseq_4_S7_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# },
# "4Cseq_K_newseq_3_S6_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# },
# "4Cseq_K_newseq_2_S5_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# },
# "4Cseq_K_newseq_1_S4_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_K.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_K.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_K.FA"]
# }
# "4Cseq_M1_S1_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_M.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_M.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_M.FA"]
# },
"4Cseq_M2_S2_all_R1_001_cutadapt":{
			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_M.FA",
			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_M.FA"],
			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_M.FA"]
},
# "4Cseq_M3_S3_all_R1_001_cutadapt":{
# 			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_M.FA",
# 			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_M.FA"],
# 			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_M.FA"]
# },
"4Cseq_M4_S4_all_R1_001_cutadapt":{
			"primers":PRIMERS_PATH+"primer_file_demultiplexing_mapping_M.FA",
			"viewpoint":PRIMERS["primer_file_demultiplexing_mapping_M.FA"],
			"to.exclude":EXCLUDES["primer_file_demultiplexing_mapping_M.FA"]
}



}

FILES =  {z+"_"+INFO_FASTQ[z]["viewpoint"][i].replace(":","_"):{"viewpoint":INFO_FASTQ[z]["viewpoint"][i],"to.exclude":INFO_FASTQ[z]["to.exclude"][i]} for z in INFO_FASTQ.keys() for i in range(0,len(INFO_FASTQ[z]["viewpoint"]))}

FILE_BY_VIEWPOINT = {}
for f in FILES.keys():
	viewpoint = FILES[f]["viewpoint"].replace(":","_")
	if viewpoint in FILE_BY_VIEWPOINT.keys():
		FILE_BY_VIEWPOINT[viewpoint].append(f)
	else:
		FILE_BY_VIEWPOINT[viewpoint] = [f]

# def getPrimerFile(wildcards):
# 	return(INFO_FASTQ[wildcards.sample]["primers"])

def getPrimers(wildcards):
	return(INFO_FASTQ[wildcards.prefix]["viewpoint"])


def getExclude(wildcards):
	return(FILES[wildcards.sample]["to.exclude"])

def getVP(wildcards):
	return(FILES[wildcards.sample]["viewpoint"])

def getFileForViewpoint(wildcards):
	return(expand(OUT_PATH+"/COV/{prefix}.count",prefix=FILE_BY_VIEWPOINT[wildcards.viewpoint]))

def getConditionSampleName(wildcards):
	return(CONDITIONS[wildcards.cond]["sample_name"])

def getConditionConditionName(wildcards):
	return(CONDITIONS[wildcards.cond]["condition"])

rule all:
	input:
		#expand(OUT_PATH+"/DE/{bin}/{viewpoint}_{bin}_{extend}_DE_{pval}_{cond}.tsv",viewpoint=['chr1:89455867-89456712', 'chr17:57168614-57169531', 'chr21:33251469-33252587', 'chr17:45759770-45760603'],bin=BINS,extend=EXTENSION_VIEW,pval=PVALS,cond=CONDITIONS.keys())
		#expand(OUT_PATH+"/COUNTS/{viewpoint}_{extend}.tsv",viewpoint=FILE_BY_VIEWPOINT.keys(),extend=EXTENSION_VIEW,pval=PVALS,cond=CONDITIONS.keys())
		expand(OUT_PATH+"/COUNTS/{viewpoint}.tsv",viewpoint=FILE_BY_VIEWPOINT.keys()),
		expand(OUT_PATH+"/BAM/{prefix}.sorted.bai",prefix=FILES.keys()),
		#expand(OUT_PATH+"/BIGWIG_SMOOTHED/{bin}/{prefix}_normalized_{bin}.bw",prefix=FILES.keys(),bin = [50000,20000,10000]),
		expand(OUT_PATH+"/GAELLE_SMOOTHED/{bin}/{prefix}_normalized.{bin}.{sw}.bw",prefix=FILES.keys(),bin = [100000,50000,20000,5000,2000],sw = 200)

rule demultiplex:
	input:
		expand(RAW_PATH+"{prefix}"+EXTENSION,prefix=INFO_FASTQ.keys())
	output:
		expand(DEMULTIPLEX_PATH+"{sample}"+EXTENSION,sample=FILES.keys())
	params:
		script=demultiplex_script,
		outdir=DEMULTIPLEX_PATH
	run:
		for f in INFO_FASTQ.keys():
			my_input=RAW_PATH+f+EXTENSION
			if False in [os.path.exists(DEMULTIPLEX_PATH+f+"_"+vp+EXTENSION) for vp in INFO_FASTQ[f]["viewpoint"]]:
				barcode=INFO_FASTQ[f]["primers"]
				shell("python2 {params.script} --fastq "+my_input+" --barcode "+barcode+" --out {params.outdir}")
		shell("rename 's/:/_/' {params.outdir}/*.fastq.gz")

rule bwa_mem:
	input:
		fastq=DEMULTIPLEX_PATH+"{sample}"+EXTENSION,
		genome=GENOME
	output:
		bam=OUT_PATH+"/BAM/{sample}.bam"
	threads:THREADS
	shell:
		"bwa mem -M -t {threads} "
		"{input.genome} "
		"{input.fastq} | "
		"samtools view -q 20 -@ {threads} -Sb -h -t {input.genome} -o {output.bam} -"


rule samtools_sort:
	input:
		bam=OUT_PATH+"/BAM/{sample}.bam"
	output:
		bam=OUT_PATH+"/BAM/{sample}.sorted.bam"
	threads: THREADS
	message: "##RUNNING : samtools sort {input.bam}"
	shell:
		"samtools sort -@ {threads} "
		"-o {output} "
		"{input.bam}"

rule samtools_index:
	input:
		bam=OUT_PATH+"/BAM/{sample}.sorted.bam"
	output:
		bam=OUT_PATH+"/BAM/{sample}.sorted.bai"
	threads: THREADS
	message: "##RUNNING : samtools index {input}"
	shell:
		"samtools index -@ {threads} {input} {output}"

rule create_frags:
	output:
		good_bed=OUT_PATH+"/FRAGMENTS/frag_MboI_NlaIII_hg19_with_motif_excluded.bed",
		bad_bed=OUT_PATH+"/FRAGMENTS/frag_MboI_or_NlaIII_hg19_with_motif_excluded.bed"
	params:
		dir = os.getcwd()
	script:
		make_frags

rule make_coverage:
	input:
		bam=OUT_PATH+"/BAM/{sample}.bam",
		fragment=OUT_PATH+"/FRAGMENTS/frag_MboI_NlaIII_hg19_with_motif_excluded.bed"
	output:
		OUT_PATH+"/COV/{sample}.count"
	shell:
		"bedtools coverage "
		"-F 0.9 -counts -a {input.fragment} "
		"-b {input.bam} "
		"> {output}"

rule create_matrix_count:
	input:
		getFileForViewpoint
	output:
		OUT_PATH+"/COUNTS/{viewpoint}.tsv"
	params:
		dir = os.getcwd(),
		inpath=OUT_PATH+"/COV",
		viewpoint="{viewpoint}"
	script:
		create_matrix_count

rule generate_bedGraph_smoothed:
	input:
		OUT_PATH+"/COV/{sample}.count"
	output:
		OUT_PATH+"/BEDGRAPH/{sample}.bedGraph"
	params:
		dir = os.getcwd(),
		nFragsPerWin=11,
		curName="{sample}",
		regToExclude=getExclude,
		scoreCol=7
	script:
		generate_bedGraph_smoothed


rule generate_bedGraph_normalized:
	input:
		OUT_PATH+"/BEDGRAPH/{sample}.bedGraph"
	output:
		OUT_PATH+"/BEDGRAPH_NORMALIZED/{sample}_normalized.bedGraph"
	params:
		dir = os.getcwd(),
		curName="{sample}",
		viewpoint=getVP,
		extSize= 2000000,
		regToExclude=getExclude
	script:
		generate_bedGraph_norm

rule bedGraphToBigWig:
	input:
		chromsizes = "/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/GENOMES/female.hg19/hg19.chrom.sizes",
		bedgraph = OUT_PATH+"/BEDGRAPH_NORMALIZED/{sample}_normalized.bedGraph"
	output:
		OUT_PATH+"/BIGWIG_NORMALIZED/{sample}_normalized.bw"
	shell:
		"bedGraphToBigWig {input.bedgraph} {input.chromsizes} {output}"

rule SmoothBWusingdeeptools:
	input:
		OUT_PATH+"/BIGWIG_NORMALIZED/{sample}_normalized.bw"
	output:
		bedgraph = temp(OUT_PATH+"/BIGWIG_SMOOTHED/{bin}/{sample}_normalized_{bin}_withna.bedGraph"),
		npz = temp(OUT_PATH+"/BIGWIG_SMOOTHED/{sample}_normalized_{bin}_withna.npz")
	params:
		bins = "{bin}"
	shell:
		"multiBigwigSummary bins -b {input} -out {output.npz} -bs={params.bins} --outRawCounts {output.bedgraph}"

rule removeNA:
	input:bedgraph = OUT_PATH+"/BIGWIG_SMOOTHED/{bin}/{sample}_normalized_{bin}_withna.bedGraph"
	output:temp(OUT_PATH+"/BIGWIG_SMOOTHED/{bin}/{sample}_normalized_{bin}_withoutna.bedGraph")
	shell : "grep -v 'nan' {input} | tail -n+2 > {output}"

rule bedGraphToBigWigSmoothed:
	input:
		chromsizes = "/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/GENOMES/female.hg19/hg19.chrom.sizes",
		bedgraph = OUT_PATH+"/BIGWIG_SMOOTHED/{bin}/{sample}_normalized_{bin}_withoutna.bedGraph"
	output:
		OUT_PATH+"/BIGWIG_SMOOTHED/{bin}/{sample}_normalized_{bin}.bw"
	shell:
		"bedGraphToBigWig {input.bedgraph} {input.chromsizes} {output}"

rule favorite_gaelle_smooth:
	input:
		OUT_PATH+"/BIGWIG_NORMALIZED/{sample}_normalized.bw"
	output:
		OUT_PATH+"/GAELLE_SMOOTHED/{bin}/{sample}_normalized.{bin}.{sw}.bw"
	params:
		dir = os.getcwd(),
		smooth_windows = "{bin}",
		smooth_subwindows = "{sw}"
	script:
		"favorite_gaelle_smooth.R"
