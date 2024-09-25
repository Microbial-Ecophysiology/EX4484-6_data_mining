rule bin_refinement:
	input:
		concoct_done = config["adir"] + "/" + config["study"] + "/binning/concoct/done",
		metabat_done = config["adir"] + "/" + config["study"] + "/binning/metabat/done"
	output:
		checkm = config["adir"] + "/" + config["study"] + "/refined_bins/metawrap_" + str(config["completeness"]) + "_" + str(config["contamination"]) + "_bins.stats"
	params:
		script_dir = config["metawrap_scripts"],
		indir = config["adir"] + "/" + config["study"] + "/binning",
		outdir = config["adir"] + "/" + config["study"] + "/refined_bins",
		comp = config["completeness"],
		cont = config["contamination"]
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	log: config["adir"] + "/" + config["study"] + "refined_bins/refinement.log"
	shell:
		"""
		metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.indir}/concoct/concoct_output/fasta_bins -B {params.indir}/metabat/metabat_out -c {params.comp} -x {params.cont} -m 370 > {log} 2>&1
		"""

checkpoint gtdbtk:
	input:
		checkm = config["adir"] + "/" + config["study"] + "/refined_bins/metawrap_" + str(config["completeness"]) + "_" + str(config["contamination"]) + "_bins.stats"
	output:
		ar53 = config["adir"] + "/" + config["study"] + "/gtdbtk_out/gtdbtk.ar53.summary.tsv",
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins/done",
		selectdir = directory(config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins")
	params:
		indir = config["adir"] + "/" + config["study"] + "/refined_bins/metawrap_" + str(config["completeness"]) + "_" + str(config["contamination"]) + "_bins",
		outdir = config["adir"] + "/" + config["study"] +  "/gtdbtk_out"
	threads: config["parallel_threads"] 
	conda: config["wdir"] + "/envs/gtdbtk.yaml"
	log: config["adir"] + config["adir"] + "/" + config["study"] + "/gtdbtk_out/gtdbtk.log"
	shell:
		"""
		gtdbtk classify_wf --genome_dir {params.indir} --out_dir {params.outdir} --cpus {threads} -x fa >> {log} 2>&1
		grep 'Thermoplasmatota' {output.ar53} | cut -f1 > {params.outdir}/Thermoplasmatota_MAGs.txt
		while read line
		do
		  cp {params.indir}/${{line}}.fa {output.selectdir} 
		done < {params.outdir}/Thermoplasmatota_MAGs.txt
		touch {output.done}
		"""

rule cat_clean_reads:
	input:
		clean_R1 = config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_1.fastq.gz",
		clean_R2 = config["adir"] + "/" + config["study"] + "/clean_reads/{sample}_2.fastq.gz"
	output:
		concat_R1 = config["adir"] + "/" + config["study"] + "/reassembled_bins/combined_reads_1.fastq",
		concat_R2 = config["adir"] + "/" + config["study"] + "/reassembled_bins/combined_reads_2.fastq"
	shell:
		"""
		zcat {input.clean_R1} > {output.concat_R1}
		zcat {input.clean_R2} > {output.concat_R2}
		"""

rule bin_reassembly_index:
	input:
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins/done" 
	output:
		allbins = config["adir"] + "/" + config["study"] + "/reassembled_bins/bin_index/all_bins.fa",
		binindex = config["adir"] + "/" + config["study"] + "/reassembled_bins/bin_index/all_bins.fa.sa"
	params:
		selectdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins",
		bindir = config["adir"] + "/" + config["study"] + "/reassembled_bins/original_bins"
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.bindir}
		ls -1 {params.selectdir} *.fa | while read line
		do
		  BIN=$(echo $line | xargs basename | sed 's/\\.fa//')
		  sed "/^>/s/^>/>${{BIN}}_/" $line > {params.bindir}/${{BIN}}.fa
		done 
		cat {params.bindir}/*.fa > {output.allbins}
		bwa index {output.allbins}
		"""

# remove supplementary and secondary alignments from bam file based on sam flag
rule bin_reassembly_mapping:
	input:
		concat_R1 = config["adir"] + "/" + config["study"] + "/reassembled_bins/combined_reads_1.fastq",
		concat_R2 = config["adir"] + "/" + config["study"] + "/reassembled_bins/combined_reads_2.fastq",
		allbins = config["adir"] + "/" + config["study"] + "/reassembled_bins/bin_index/all_bins.fa",
		binindex = config["adir"] + "/" + config["study"] + "/reassembled_bins/bin_index/all_bins.fa.sa"
	output:
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly/done"
	params:
		script_dir = config["metawrap_scripts"],
		indir = config["adir"] + "/" + config["study"] + "/reassembled_bins/original_bins",
		outdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.outdir}
		# ulimit -n 10000
		bwa mem -t {threads} {input.allbins} {input.concat_R1} {input.concat_R2} | samtools view -F 256 - | {params.script_dir}/filter_reads_for_bin_reassembly.py {params.indir} {params.outdir} 2 5
		touch {output.done}
		"""


rule bin_reassembly_strict:
	input:
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly/done",
		fa = config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins/{tbin}.fa"
	output:
		fa = config["adir"] + "/" + config["study"] + "/reassembled_bins/{tbin}.strict.fa"
	params:
		script_dir = config["metawrap_scripts"],
		readdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly",
		tmpdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/tmp_spades",
		spadesdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassemblies",
		bin = "{cbin}"
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.spadesdir}
		# if spades fails, just copy the original assembly (this should avoid breaking the snakemake workflow)
		spades.py -t {threads} -m 150 --tmp {params.tmpdir}/tmp_{params.bin}.strict --careful --untrusted-contigs {input.fa} -1 {params.readdir}/{params.bin}.strict_1.fastq -2 {params.readdir}/{params.bin}.strict_2.fastq -o {params.spadesdir}/{params.bin}.strict || sed '/^>/s/^[^_]*_/>/' {input.fa} > {params.spadesdir}/{params.bin}.strict/scaffolds.fasta
		if [[ -s {params.spadesdir}/{params.bin}.strict/scaffolds.fasta ]]
		then
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.strict/scaffolds.fasta > {output.fa}
		else
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.strict/contigs.fasta > {output.fa}
		fi
		"""

rule bin_reassembly_permissive:
	input:
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly/done",
		fa = config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins/{tbin}.fa"
	output:
		fa = config["adir"] + "/" + config["study"] + "/reassembled_bins/{tbin}.permissive.fa"
	params:
		script_dir = config["metawrap_scripts"],
		readdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly",
		tmpdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/tmp_spades",
		spadesdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassemblies",
		bin = "{cbin}"
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.spadesdir}
		# if spades fails, just copy the original assembly (this should avoid breaking the snakemake workflow)
		spades.py -t {threads} -m 150 --tmp {params.tmpdir}/tmp_{params.bin}.permissive --careful --untrusted-contigs {input.fa} -1 {params.readdir}/{params.bin}.permissive_1.fastq -2 {params.readdir}/{params.bin}.permissive_2.fastq -o {params.spadesdir}/{params.bin}.permissive || sed '/^>/s/^[^_]*_/>/' {input.fa} > {params.spadesdir}/{params.bin}.permissive/scaffolds.fasta
		if [[ -s {params.spadesdir}/{params.bin}.permissive/scaffolds.fasta ]]
		then
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.permissive/scaffolds.fasta > {output.fa}
		else
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.permissive/contigs.fasta > {output.fa}
		fi
		"""

rule bin_reassembly_orig:
	input:
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/reads_for_reassembly/done",
		fa = config["adir"] + "/" + config["study"] + "/reassembled_bins/selected_bins/{tbin}.fa"
	output:
		fa = config["adir"] + "/" + config["study"] + "/reassembled_bins/{tbin}.orig.fa"
	shell:
		"""
		cp {input.fa} {output.fa}
		"""

# I am using the wildcard name cbin here (candidate bin), to not cause any conflicts with the wildcard bin used later
def aggregate_strict(wildcards):
	checkpoint_output = checkpoints.gtdbtk.get().output['selectdir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	tbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/" + config["study"] + "/reassembled_bins/{tbin}.strict.fa", tbin = tbin)
	
def aggregate_permissive(wildcards):
	checkpoint_output = checkpoints.gtdbtk.get().output['selectdir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	tbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/" + config["study"] + "/reassembled_bins/{tbin}.permissive.fa", tbin = tbin)
	
def aggregate_orig(wildcards):
	checkpoint_output = checkpoints.gtdbtk.get().output['selectdir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	tbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/" + config["study"] + "/reassembled_bins/{tbin}.orig.fa", tbin = tbin)


rule bin_reassembly_checkm:
	input:
		strict = aggregate_strict,
		permissive = aggregate_permissive,
		orig = aggregate_orig
	output:
		stats = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassembled_bins.stats"
	params:
		script_dir = config["metawrap_scripts"],
		tmpdir = config["adir"] + "/tmp",
		checkmdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassembled_bins.checkm",
		indir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassembled_bins"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.tmpdir}
		checkm lineage_wf -x fa {params.indir} {params.checkmdir} -t {threads} --tmpdir {params.tmpdir} --pplacer_threads 10
		{params.script_dir}/summarize_checkm.py {params.checkmdir}/storage/bin_stats_ext.tsv | (read -r; printf "%s\\n" "$REPLY"; sort) > {output.stats}
		"""

# settings for parameters a and b for calculating score to select best bin have been empirically chosen
# they should not be changed independently from minimum completeness and maximum contamination settings
rule bin_reassembly_best:
	input:
		stats = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassembled_bins.stats"
	output:
		best_bin_list = config["adir"] + "/" + config["study"] + "/reassembled_bins/best_bins.list"
	params:
		script_dir = config["metawrap_scripts"],
		comp = config["completeness"],
		cont = config["contamination"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		{params.script_dir}/choose_best_bin.py {input.stats} {params.comp} {params.cont} > {output.best_bin_list}
		"""

rule collect_final_bins:
	input:
		best_bin_list = config["adir"] + "/" + config["study"] + "/reassembled_bins/best_bins.list"
	output:
		done = config["adir"] + "/" + config["study"] + "/reassembled_bins/bin_final/done",
	params:
		outdir = config["adir"] + "/" + config["study"] + "/reassembled_bins/bin_final",
		indir = config["adir"] + "/" + config["study"] + "/reassembled_bins/reassembled_bins"
	shell:
		"""
		mkdir -p {params.outdir}
		cat {input.best_bin_list} | while read line
		do 
		  sed -e "/^>/s/^>/>${{line}}_/" -e '/^>/s/[-\\.]/_/g' {params.indir}/${{line}}.fa > {params.outdir}/${{line}}.fa
		done
		touch {output.done}
		"""
