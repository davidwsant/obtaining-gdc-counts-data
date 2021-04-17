#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from argparse import ArgumentParser
import sys
import glob
import os
import tarfile

args = ArgumentParser('./merge_counts.py', description="""This program has been designed to merge
count files obtained from the Genomic Data Commons (https://portal.gdc.cancer.gov/). Please
specify the 'tar.gz' file that contains all of the count files using the -c parameter and the
samplesheet using the -s parater. The unpacked clinical cart file can optionally be saved using
the -u parameter.
Example usage: python merge_counts.py -s gdc_sample_sheet.2021-04-13.tsv -c clinical.cart.2021-04-13.tar.gz
-u Unpacked_GDC """)

args.add_argument(
	'-c',
	'--cart_file',
	help="""This is the file that contains count data from Genomic Data Commons. The file name should end
	in 'tar.gz'. Unless otherwise specified, Python will unzip this file and put the resulting 
        subdirectories in a folder titled 'GDC_unpack_tmp' that will be deleted when the files are 
        finished merging.""",
	default=None,
)

args.add_argument(
	'-s',
	'--samplesheet',
	help="""This is the name of the samplesheet that you downloaded from GDC. The file name should start
	with 'gdc_sample_sheet' followed by a date and end in '.tsv'. """,
	default = None
)

args.add_argument(
	'-u',
	'--unpacked_folder',
	help="""Use this option to specify the name of a folder to store the unpacked 'tar.gz' file
	if you would like it saved. Individual count files will remain zipped. """,
	default = None
)

args.add_argument(
	'-o',
	'--output_prefix',
	help="""You can use this option to specify a prefix that you would like added to the output files.
	By default, no prefix will be added and the output file names will reflect the sample type. """,
	default = None
)


args = args.parse_args()


def print_error_message():
	print()
	print("\tWelcome to merge_counts.py.")
	print("\tThis program has been designed to merge count files obtained from the Genomic Data Commons ")
	print("\t(https://portal.gdc.cancer.gov/). Please specify the 'tar.gz' file that contains all of the count ")
	print("\tfiles using the -c parameter and the samplesheet using the -s parater. The unpacked clinical cart file ")
	print("\tcan optionally be saved using the -u parameter.")
	print()
	print("\tExample usage: python merge_counts.py -s gdc_sample_sheet.2021-04-13.tsv")
	print("\t-c clinical.cart.2021-04-13.tar.gz -u Unpacked_GDC")
	print()

cart_file = args.cart_file
samplesheet = args.samplesheet
unpacked_folder = args.unpacked_folder
output_prefix = args.output_prefix


if not cart_file:
	files = glob.glob('*.tar.gz')
	print_error_message()
	print("\tYou have not specified a gzipped tar file containing the count files that you wish to merge.")
	print("\tPlease specify this file using the -c option.")
	print()
	print("\tPossible files from your current working directory are: ")
	for file in files:
		print('\t'+file)
	print()
	sys.exit(1) # Exit with a status of 1.

if not cart_file.endswith('tar.gz'):
	files = glob.glob('*.tar.gz')
	print_error_message()
	print("\tThe cart file that you have specified does not appear to be a gzipped tar file.")
	print("\tPlease specify this file using the -c option.")
	print()
	print("\tPossible files from your current working directory are: ")
	for file in files:
		print('\t'+file)
	print()
	sys.exit(1) # Exit with a status of 1.

if not samplesheet:
	files = glob.glob('*.tsv')
	print_error_message()
	print("\tYou have not specified a samplesheet containing information about the individual count files.")
	print("\tPlease specify this file using the -s option.")
	print()
	print("\tPossible files from your current working directory are: ")
	for file in files:
		print('\t'+file)
	print()
	sys.exit(1) # Exit with a status of 1.

if not samplesheet.endswith('tsv'):
	files = glob.glob('*.tsv')
	print_error_message()
	print("\tThe samplesheet file that you have specified does not appear to be a tsv file provided by GDC.")
	print("\tPlease specify this file using the -s option.")
	print()
	print("\tPossible files from your current working directory are: ")
	for file in files:
		print('\t'+file)
	print()
	sys.exit(1) # Exit with a status of 1.

unpack_name = 'GDC_unpack_tmp'
if unpacked_folder:
	unpack_name = unpacked_folder

os.mkdir(unpack_name)
os.chdir(unpack_name)
tar = tarfile.open('../'+cart_file, "r:gz")
tar.extractall()
tar.close()
os.chdir('../')

samplesheet_df = pd.read_csv(samplesheet, sep = "\t")

primary_df = pd.DataFrame()
metastaitc_df = pd.DataFrame()
solid_normal_df = pd.DataFrame()
other_type_df = pd.DataFrame()

def merge_dataframes(df, tmp_df):
	if len(df) == 0:
		df = pd.concat([df, tmp_df])
	else:
		df = df.merge(tmp_df, left_index = True, right_index=True)
	return df

sample_names = {'primary': {},
	'metastatic': {},
	'solid_normal': {},
	'other': {}}

for i, row in samplesheet_df.iterrows():
	file_name = row['File Name']
	file_path = unpack_name+'/'+row['File ID']+'/'+file_name
	tmp_df = pd.read_csv(file_path, sep="\t", index_col = 0, header = None)
	sample_name = row['Case ID']
	sample_type = row['Sample Type']

	if sample_type == 'Primary Tumor':
		if sample_name not in sample_names['primary']:
			tmp_df.columns = [sample_name]
			sample_names['primary'][sample_name] = 1
		else:
			sample_names['primary'][sample_name] += 1
			count = sample_names['primary'][sample_name]
			new_name = sample_name+'_'+str(count)
			tmp_df.columns = [new_name]
			print(sample_name +' has more than one Primary Tumor file in this directory. ')
			print('The file '+file_name+' will be given the title '+new_name+'.')
		primary_df = merge_dataframes(primary_df, tmp_df)

	elif sample_type == 'Metastatic':
		if sample_name not in sample_names['metastatic']:
			tmp_df.columns = [sample_name]
			sample_names['metastatic'][sample_name] = 1
		else:
			sample_names['metastatic'][sample_name] += 1
			count = sample_names['metastatic'][sample_name]
			new_name = sample_name+'_'+str(count)
			tmp_df.columns = [new_name]
			print(sample_name +' has more than one Metastatic Tumor file in this directory. ')
			print('The file '+file_name+' will be given the title '+new_name+'.')
		metastaitc_df = merge_dataframes(metastaitc_df, tmp_df)

	elif sample_type == 'Solid Tissue Normal':
		if sample_name not in sample_names['solid_normal']:
			tmp_df.columns = [sample_name]
			sample_names['solid_normal'][sample_name] = 1
		else:
			sample_names['solid_normal'][sample_name] += 1
			count = sample_names['solid_normal'][sample_name]
			new_name = sample_name+'_'+str(count)
			tmp_df.columns = [new_name]
			print(sample_name +' has more than one Solid Tissue Normal file in this directory. ')
			print('The file '+file_name+' will be given the title '+new_name+'.')
		solid_normal_df = merge_dataframes(solid_normal_df, tmp_df)

	else:
		if sample_name not in sample_names['other']:
			tmp_df.columns = [sample_name]
			sample_names['other'][sample_name] = 1
		else:
			sample_names['other'][sample_name] += 1
			count = sample_names['other'][sample_name]
			new_name = sample_name+'_'+str(count)
			tmp_df.columns = [new_name]
			print(sample_name +' has more than one file in this directory that did not fit into a TCGA classification.')
			print('The file '+file_name+' will be given the title '+new_name+'.')
		other_type_df = merge_dataframes(other_type_df, tmp_df)

# All dataframes have been created, now just to save out the files with the correct prefix

primary_name = "Primary_Tumors_counts.csv"
metastatic_name = "Metastatic_Tumors_counts.csv"
solid_name = "Solid_Tissue_Normal_counts.csv"
other_name = "Non_TCGA_Tissue_Type_counts.csv"
if output_prefix:
	primary_name = output_prefix+'_'+primary_name
	metastatic_name = output_prefix+'_'+metastatic_name
	solid_name = output_prefix+'_'+solid_name
	other_name = output_prefix+'_'+other_name

if len(primary_df) > 0:
	primary_df.to_csv(primary_name)
if len(metastaitc_df) > 0:
	metastaitc_df.to_csv(metastatic_name)
if len(solid_normal_df) > 0:
	solid_normal_df.to_csv(solid_name)
if len(other_type_df) > 0:
	other_type_df.to_csv(other_name)

if not unpacked_folder:
	import shutil
	try:
		shutil.rmtree('GDC_unpack_tmp')
	except OSError as e:
		print('\tCount files have been merged, but unable to delete "GDC_unpack_tmp" folder.')
		print("\tError: %s : %s" % (dir_path, e.strerror))
