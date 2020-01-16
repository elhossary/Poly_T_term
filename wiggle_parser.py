import pandas as pd
from io import StringIO


class WiggleParser:

	def __init__(self, file_path):
		self.file_path = file_path

	def parse(self):
		chrom_accession = ""
		chrom_wig_vals_dict = {}
		temp_wig_vals = ""
		for line in open(self.file_path).readlines():
			if line[0].isnumeric():
				temp_wig_vals += line
			else:
				if "chrom=" in line:
					if temp_wig_vals == "":
						chrom_accession = line.split('chrom=')[-1].split(' ')[0]
						continue
					print(f"\tLoading coverage for: {chrom_accession}")
					chrom_wig_vals_dict[chrom_accession] = pd.read_csv(StringIO(temp_wig_vals), sep=" ", header=None)
					temp_wig_vals = ""
					chrom_accession = line.split('chrom=')[-1].split(' ')[0]
				if "name=" in line:
					wig_name = line.split('name=')[-1].replace('\n', '')
					print(f"Parsing file: {wig_name}")
		print(f"\tLoading coverage for: {chrom_accession}")
		chrom_wig_vals_dict[chrom_accession] = pd.read_csv(StringIO(temp_wig_vals), sep=" ", header=None)

		return chrom_wig_vals_dict
