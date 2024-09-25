rule download_study:
	input: 
		link = config["adir"] + config["analysis"] + "/links/{study}.txt",
		link_log = config["adir"] + config["analysis"] + "/links/links.log"
	params:
		outdir = config["adir"] + config["analysis"] + "/downloaded"
	conda: config["wdir"] + "/envs/aria.yaml"
	log: config["adir"] + config["analysis"] + "logfiles/MGs_download_{study}.log"
	shell:
		"""
		aria2c -i {input.link} -c -l {input.link_log} --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads=100 >> {log} 2>&1
		"""
