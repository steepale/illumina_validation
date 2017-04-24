import os
import sys
import re

# infile (open independently)
infile = sys.argv[1]

# outfile (open once and close once)
outfile = open(sys.argv[2], 'w')

# reference files
orthologue_file = "/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/ensembl/chicken_human_orthologues_full_annotation.txt"
cosmic_cgc_file = "/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/cosmic/cosmic_CGC_gene_list_2017_01_03.tsv"

# Create a dictionary of mutated genes to tally distinct sample counts
samples_by_gene = {}
for inline in open(infile):
	if inline[0] != '#':
		inline = inline.rstrip()
		incol = inline.split('\t')
		ingeneid = incol[7]
		insample = incol[9]
		insample_list = insample.split('|')
		if ingeneid not in samples_by_gene.keys():
			samples_by_gene[ingeneid] = set()
			samples_by_gene[ingeneid] |= set(insample_list)
		if ingeneid in samples_by_gene.keys():
			samples_by_gene[ingeneid] |= set(insample_list)

# Create a series of dictionaries of high confidence orthologues:
# Create a dictionary of high confidence orthologues of chicken ensembl gene id 2 human ensembl gene id
chicken2human_ortho_ensembl_gene_id = {}
# Create a dictionary of high confidence orthologues of chicken ensembl gene id 2 human ensembl gene name
chicken_ensembl_gene_id2human_ensembl_gene_name_ortho = {}
for o_line in open(orthologue_file):
	if o_line[0] != '#':
		o_line = o_line.rstrip()
		o_col = o_line.split('\t')
		c_e_gene_id = o_col[0]
		c_e_gene_name = o_col[1]
		c_e_gene_desc = o_col[2]
		h_e_gene_id = o_col[3]
		h_e_gene_name = o_col[4]
		homology = o_col[5]
		id_gene_h2c = o_col[6]
		id_gene_c2h = o_col[7]
		goc_score = o_col[8]
		wga_cov = o_col[9]
		ortho_score = o_col[10]
		# Fill dictionary of high confidence orthologues of chicken ensembl gene id 2 human ensembl gene id
		if c_e_gene_id not in chicken2human_ortho_ensembl_gene_id.keys() and ortho_score == '1':
			chicken2human_ortho_ensembl_gene_id[c_e_gene_id] = set()
			chicken2human_ortho_ensembl_gene_id[c_e_gene_id].add(h_e_gene_id)
		if c_e_gene_id in chicken2human_ortho_ensembl_gene_id.keys() and ortho_score == '1':
			chicken2human_ortho_ensembl_gene_id[c_e_gene_id].add(h_e_gene_id)
		# Fill dictionary of high confidence orthologues of chicken ensembl gene id 2 human ensembl gene name
		if c_e_gene_id not in chicken_ensembl_gene_id2human_ensembl_gene_name_ortho.keys() and ortho_score == '1':
			chicken_ensembl_gene_id2human_ensembl_gene_name_ortho[c_e_gene_id] = set()
			chicken_ensembl_gene_id2human_ensembl_gene_name_ortho[c_e_gene_id].add(h_e_gene_name)
		if c_e_gene_id in chicken_ensembl_gene_id2human_ensembl_gene_name_ortho.keys() and ortho_score == '1':
			chicken_ensembl_gene_id2human_ensembl_gene_name_ortho[c_e_gene_id].add(h_e_gene_name)


# Create a set of mutated genes in cancer gene consensus
cgc_gene_set = set()
for cgc_line in open(cosmic_cgc_file):
	cgc_line = cgc_line.rstrip()
	cgc_col = cgc_line.split('\t')
	cgc_gene = cgc_col[0]
	cgc_gene_set.add(cgc_gene)

# Create a set of final variants for output that pass filters
final_var_set = set()

# Write a header for the outfile
outfile.write('#CHROM'+'\t'+'POS'+'\t'+'REF'+'\t'+'ALT'+'\t'+'MUT'+'\t'+'IMPACT'+'\t'+'SYMBOL'+'\t'+'GENE_ID'+'\t'+'ORTHOLOGUE'+'\t'+'TSN_VAR'+'\t'+'TSN_GENE'+'\t'+'SAMPLE'+'\t'+'VAC'+'\t'+'VAF'+'\t'+'NUM_TOOLS'+'\t'+'CGC_STATUS'+'\t'+'FILTER'+'\n')

# Iterate through the infile and collect additional information for each filter
for inline in open(infile):
	if inline[0] != '#':
		inline = inline.rstrip()
		incol = inline.split('\t')
		inchr = incol[0]
		inpos = incol[1]
		inref = incol[2]
		inalt = incol[3]
		var_id = inchr+inpos+inref+inalt
		inmut = incol[4]
		inimpact = incol[5]
		insymbol = incol[6]
		ingeneid = incol[7]
		intsn_by_var = incol[8]
		intsn_by_gene = str(len(samples_by_gene[ingeneid]))
		insample = incol[9].replace('|', ';')
		invac = incol[10].replace('|', ';')
		invaf = incol[11].replace('|', ';')
		# Annotate validation site ID
		# Annotate gene if it has a high confidence orthologue
		if ingeneid in chicken_ensembl_gene_id2human_ensembl_gene_name_ortho.keys():
			ortho_status = 'Yes'
			ortho_set = chicken_ensembl_gene_id2human_ensembl_gene_name_ortho[ingeneid]
			orthologue = ';'.join(map(str,ortho_set))
		if ingeneid not in chicken_ensembl_gene_id2human_ensembl_gene_name_ortho.keys():
			ortho_status = 'No'
			orthologue = 'NA'
		# Annotate gene if found in cancer gene consensus
		if ingeneid in chicken_ensembl_gene_id2human_ensembl_gene_name_ortho.keys():
			for human_ensembl_gene_name in chicken_ensembl_gene_id2human_ensembl_gene_name_ortho[ingeneid]:
				if human_ensembl_gene_name in cgc_gene_set:
					cgc_status = 'Yes'
				elif human_ensembl_gene_name not in cgc_gene_set:
					cgc_status = 'No'
		# Annotate variant by consenus of variant callers (need to extract information from vep files)
		callers = []
		vep_samples = set(incol[9].split('|'))
		# If the variant is an snv
		if len(inref) == 1 and len(inalt) == 1:
			var_type = 'snv'
		elif len(inref) != 1 or len(inalt) != 1:
			var_type = 'indel'
		for vep_sample in vep_samples:
			vep_file = "./data/somaticseq_vcf/"+vep_sample+"_somaticseq_"+var_type+"_vep.vcf"
			for vep_line in open(vep_file):
				if vep_line[0] != '#':
					vep_line = vep_line.rstrip()
					vep_cols = vep_line.split('\t')
					vep_chr = vep_cols[0]
					vep_pos = vep_cols[1]
					vep_ref = vep_cols[3]
					vep_alt = vep_cols[4]
					vep_var = vep_chr+vep_pos+vep_ref+vep_alt
					vep_info = vep_cols[7]
					if vep_var == var_id:
						if re.search('NUM_TOOLS', vep_info.split(';')[1]):
							vep_tool_num = vep_info.split(';')[1][-1]
						elif re.search('NUM_TOOLS', vep_info.split(';')[2]):
							vep_tool_num = vep_info.split(';')[2][-1]
			callers.append(vep_tool_num)
			alg_num = ';'.join(map(str,callers))
		# Filter 1: Variants that occur in more than one tumor sample (filter out undesirable genes, e.g. olfactory receptors)
		if int(intsn_by_var) > 1 and ortho_status == 'Yes' and orthologue[0:2] != 'OR' and var_id not in final_var_set:
			filter_status = '1'
			final_var_set.add(var_id)
			outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
		# Filter 2: Variants in genes, which are mutated with different variants in multiple samples (filter out undesirable genes, e.g. olfactory receptors)
		if int(intsn_by_gene) > 1 and ortho_status == 'Yes' and orthologue[0:2] != 'OR' and var_id not in final_var_set:
			filter_status = '2'
			final_var_set.add(var_id)
			outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
		# Filter 3: Variants that are found in COSMIC's cancer gene consensus (filter out undesirable genes, e.g. olfactory receptors)
		if cgc_status == 'Yes' and ortho_status == 'Yes' and orthologue[0:2] != 'OR' and var_id not in final_var_set:
			filter_status = '3'
			final_var_set.add(var_id)
			outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+ insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
		# Filter 4: Variants with a consensus of gene callers greater than 1 and variant allele count greater than one
		for call_num in alg_num.split(';'):
			for vac in invac.split(';'):
				if int(call_num) > 1 and int(vac) > 1 and ortho_status == 'Yes' and orthologue[0:2] != 'OR' and var_id not in final_var_set:
					filter_status = '4'
					final_var_set.add(var_id)
					outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+ insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
		# Filter 5: Variants with a mutation that is associated with a high impact (annotaion via VEP) and variant allele count greater than 1
		if inimpact == 'HIGH' and int(vac) > 1 and ortho_status == 'Yes' and orthologue[0:2] != 'OR' and var_id not in final_var_set:
			filter_status = '5'
			final_var_set.add(var_id)
			outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+ insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
		# Filter 6: Variants with a consensus of gene callers greater than 1
		for call_num in alg_num.split(';'):
			for vac in invac.split(';'):
				if int(call_num) > 1 and ortho_status == 'Yes' and orthologue[0:2] != 'OR' and var_id not in final_var_set:
					filter_status = '6'
					final_var_set.add(var_id)
					outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+ insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
		# Filter 7: The remaining variants
		if var_id not in final_var_set:
			filter_status = '7'
			final_var_set.add(var_id)
			outfile.write(inchr+'\t'+inpos+'\t'+inref+'\t'+inalt+'\t'+inmut+'\t'+inimpact+'\t'+insymbol+'\t'+ingeneid+'\t'+orthologue+'\t'+intsn_by_var+'\t'+intsn_by_gene+'\t'+ insample+'\t'+invac+'\t'+invaf+'\t'+alg_num+'\t'+cgc_status+'\t'+filter_status+'\n')
# Close outfile
outfile.close
