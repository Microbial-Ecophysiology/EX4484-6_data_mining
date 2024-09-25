#rule remove_phix:
#	output:
#		trim1_R1 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim1_R1.fastq.gz",
#		trim1_R2 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim1_R2.fastq.gz",
#		stats1 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_stats1.txt"
#	params:
#		download = config["adir"] + "/" + config["study"] + "/downloaded",
#		sample = "{sample}"
#	threads: config["parallel_threads"]
#	conda: config["wdir"] + "/envs/bbmap.yaml"
#	log: config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim1.log"
#	shell:
#		"""
#		bbduk.sh in={params.download}/{params.sample}_1.fastq.gz in2={params.download}/{params.sample}_2.fastq.gz out={output.trim1_R1} out2={output.trim1_R2} ref=phix k=28 stats={output.stats1} threads={threads} fastawrap=300 > {log} 2>&1
#		"""

#rule remove_adapter_right:
#	input:
#		trim1_R1 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim1_R1.fastq.gz",
#		trim1_R2 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim1_R2.fastq.gz"
#	output:
#		trim2_R1 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim2_R1.fastq.gz",
#		trim2_R2 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim2_R2.fastq.gz"
#	params:
#		adapter = config["adapter"]   
#	threads: config["parallel_threads"]
#	conda: config["wdir"] + "/envs/bbmap.yaml"
#	log: config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim2.log"
#	shell:
#		"""
#		bbduk.sh in1={input.trim1_R1} in2={input.trim1_R2} out1={output.trim2_R1} out2={output.trim2_R2} ref={params.adapter} ktrim=r k=23 mink=11 tpe tbo threads={threads} > {log} 2>&1
#		"""

rule quality_trimming:
#	input:
#		trim2_R1 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim2_R1.fastq.gz",
#		trim2_R2 = config["adir"] + "/" + config["study"] + "/trimmed_reads/{sample}_trim2_R2.fastq.gz"
	output:
		clean_R1 = config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_1.fastq.gz",
		clean_R2 = config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_2.fastq.gz"
	params:
		minlen = config["min_read_length"],
		sample = "{sample}",
		download = config["adir"] + "/" + config["study"] + "/downloaded"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/bbmap.yaml"
	log: config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_trim3.log"
	shell:
		"""
		bbduk.sh in={params.download}/{params.sample}_1.fastq.gz in2={params.download}/{params.sample}_2.fastq.gz out1={output.clean_R1} out2={output.clean_R2} qtrim=rl trimq=20 minlen={params.minlen} maq=10 threads={threads} > {log} 2>&1
		"""

# bbduk.sh in1={input.trim2_R1} in2={input.trim2_R2} out1={output.clean_R1} out2={output.clean_R2} qtrim=rl trimq=20 minlen={params.minlen} maq=10 threads={threads} > {log} 2>&1
