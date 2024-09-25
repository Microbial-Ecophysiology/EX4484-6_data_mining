rule assembly_megahit:
	input:
		clean_R1 = expand(config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_1.fastq.gz", sample = SAMPLES),
		clean_R2 = expand(config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_2.fastq.gz", sample = SAMPLES)
	output:
		assembly = config["adir"] + "/" + config["study"] + "/assembly/final.contigs.fa"
	params:
		outdir = config["adir"] + "/" + config["study"] + "/assembly",
		kmer = config["k_list"]
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/megahit.yaml"
	log: config["adir"] + "/" + config["study"] + "/assembly/coassembly.log"
	shell:
		"""
		echo "{input.clean_R1}" | tr ' ' ',' > {params.outdir}/R1.csv
		echo "{input.clean_R2}" | tr ' ' ',' > {params.outdir}/R2.csv
		megahit --k-list {params.kmer} -1 $(<{params.outdir}/R1.csv) -2 $(<{params.outdir}/R2.csv) -t {threads} -m 0.98 --min-contig-len 2500 --continue -o {params.outdir}/tmp > {log} 2>&1
		mv {params.outdir}/tmp/* {params.outdir}
		"""

#rule indexing_bowtie2:
#	input:
#		assembly = config["adir"] + "/" + config["study"] + "/assembly/final.contigs.fa"
#	output:
#		contigs = config["adir"] + "/" + config["study"] + "/mapping/contigs.1.bt2"
#	params:
#		outdir = config["adir"] + "/" + config["study"] + "/mapping"
#	threads: config["parallel_threads"]
#	conda: config["wdir"] + "/envs/metawrap.yaml"
#	log: config["adir"] + "/" + config["study"] + "/mapping/indexing.log"
#	shell:
#		"""
#		bowtie2-build --threads {threads} {input.assembly} {params.outdir}/contigs >> {log} 2>&1
#		"""

rule mapping_bowtie2:
	input:
		clean_R1 = config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_1.fastq.gz",
		clean_R2 = config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_2.fastq.gz",
		contigs = config["adir"] + "/" + config["study"] + "/mapping/contigs.1.bt2"
	output:
		bam = config["adir"] + "/" + config["study"] + "/mapping/{sample}.bam",
		bai = config["adir"] + "/" + config["study"] + "/mapping/{sample}.bam.bai"
	params:
		outdir = config["adir"] + "/" + config["study"] + "/mapping",
		sample = "{sample}"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	# log: config["adir"] + "/" + config["study"] + "/mapping/{sample}_mapping.log"
	shell:
		"""
		bowtie2 --threads {threads} -x {params.outdir}/contigs -1 {input.clean_R1} -2 {input.clean_R2} -S {params.outdir}/{params.sample}.sam
		samtools sort -o {output.bam} {params.outdir}/{params.sample}.sam
		samtools index {output.bam}
		rm {params.outdir}/{params.sample}.sam
		"""

rule binning_concoct:
	input:
		assembly = config["adir"] + "/" + config["study"] + "/assembly/final.contigs.fa",
		bam = expand(config["adir"] + "/" + config["study"] + "/mapping/{sample}.bam", sample = SAMPLES)
	output:
		bed_out = config["adir"] + "/" + config["study"] + "/binning/concoct/contigs_10K.bed",
		comp_out = config["adir"] + "/" + config["study"] + "/binning/concoct/contigs_10K.fa",
		cov_table = config["adir"] + "/" + config["study"] + "/binning/concoct/coverage_table.tsv",
		done = config["adir"] + "/" + config["study"] + "/binning/concoct/done"
	params:
		mapdir = config["adir"] + "/" + config["study"] + "/mapping",
		outdir = config["adir"] + "/" + config["study"] + "/binning/concoct"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	log: config["adir"] +  "/" + config["study"] + "/binning/concoct.log"
	shell:
		"""
		cut_up_fasta.py {input.assembly} -c 10000 -o 0 --merge_last -b {output.bed_out} > {output.comp_out} >> {log} 2>&1
		concoct_coverage_table.py {output.bed_out} {params.mapdir}/*.bam > {output.cov_table} >> {log} 2>&1
		concoct --composition_file {output.comp_out} --coverage_file {output.cov_table} -b {params.outdir}/concoct_output -l 2500 -t {threads} >> {log} 2>&1
		merge_cutup_clustering.py {params.outdir}/concoct_output/clustering_gt2500.csv > {params.outdir}/concoct_output/clustering_merged.csv >> {log} 2>&1
		extract_fasta_bins.py {input.assembly} {params.outdir}/concoct_output/clustering_merged.csv --output_path {params.outdir}/concoct_output/fasta_bins >> {log} 2>&1
		touch {output.done}
		"""

rule binning_metabat:
	input:
		assembly = config["adir"] + "/" + config["study"] + "/assembly/final.contigs.fa",
		bam = expand(config["adir"] + "/" + config["study"] + "/mapping/{sample}.bam", sample = SAMPLES)
	output:
		depth = config["adir"] + "/" + config["study"] + "/mapping/depth.txt",
		done = config["adir"] + "/" + config["study"] + "/binning/metabat/done"
	params:
		mapdir = config["adir"] + "/" + config["study"] + "/mapping",
		outdir = config["adir"] + "/" + config["study"] + "/binning/metabat"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	log: config["adir"] + "/" + config["study"] + "/binning/metabat.log"
	shell:
		"""
		jgi_summarize_bam_contig_depths --outputDepth {output.depth} {params.mapdir}/*.bam
		metabat2 -i {input.assembly} -a {output.depth} -o {params.outdir}/metabat_out -m 2500 -v -t {threads}
		touch {output.done}
   		"""
