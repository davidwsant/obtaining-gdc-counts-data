# merge-augmented-counts

This program has been designed to merge count files obtained from the Genomic Data Commons
(https://portal.gdc.cancer.gov/). Files containing unstranded counts and TPM values for each
tissue type will be created. Please specify the 'tar.gz' file that contains all of the
count files using the -c parameter and the samplesheet using the -s parameter. The unpacked
clinical cart file can optionally be saved using the -u parameter.

Example usage:
```
python merge_augmented_counts.py -s gdc_sample_sheet.2021-04-13.tsv -c gdc_download_20210414_202057.715477.tar.gz -u Unpacked_GDC
```

optional arguments:

  -h, --help

                        show this help message and exit

  -c CART_FILE, --cart_file CART_FILE

                        This is the file that contains count data from Genomic Data Commons. The file name should end
                        in 'tar.gz'. Unless otherwise specified, Python will unzip this file and put the resulting
                        subdirectories in a folder titled 'GDC_unpack_tmp' that will be deleted when the files are
                        finished merging.

  -s SAMPLESHEET, --samplesheet SAMPLESHEET

                        This is the name of the samplesheet that you downloaded from GDC. The file name should start
                        with 'gdc_sample_sheet' followed by a date and end in '.tsv'.

  -u UNPACKED_FOLDER, --unpacked_folder UNPACKED_FOLDER

                        Use this option to specify the name of a folder to store the unpacked 'tar.gz' file if you
                        would like it saved. Individual count files will remain zipped.

  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX

                        You can use this option to specify a prefix that you would like added to the output files. By
                        default, no prefix will be added and the output file names will reflect the sample type.
