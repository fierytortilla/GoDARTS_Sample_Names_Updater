from os import path
import subprocess
import argparse
from numpy import nan
import pandas as pd
import re
from gzip import open as gzopen
from tqdm import tqdm
import glob
import warnings

# old_vcf_reader_path = "/Users/dravis/Documents/Dundee/AMR/Example_VCFs/chr21_3samples_head100k_dummy.vcf"
# new_vcf_path = "output_test.vcf.gz"
# GDlinker_prochi_df = pd.read_csv('/Users/dravis/Documents/Dundee/AMR/GoDARTS_Data_Converter/GDlinker_PROCHI.csv', index_col=0)

#please have bgzip and tabix preinstalled in order for this code to work
#bgzip and tabix can be downloaded from http://www.htslib.org/download/


def update_prochi_sample_names(old_vcf_reader, prochi_linker_filepath, new_vcf_writer_filepath):
    current_directory = path.dirname(path.realpath(__file__))
    new_vcf_writer_filepath = new_vcf_writer_filepath.rstrip(".vcf.gz")
    non_corresponding_vcf_writer_name = new_vcf_writer_filepath+"_non_corresponding_samples_output.vcf"
    prochi_linker_dict = pd.read_csv(prochi_linker_filepath, index_col=0)
    corresponding_samples_dict= {}
    non_corresponding_samples_dict = {}
    vcf_header= []
    new_vcf_writer = open(new_vcf_writer_filepath + ".vcf", 'w')
    for line in tqdm(old_vcf_reader):
        if line.startswith('##'):
            vcf_header.append(line)
            new_vcf_writer.write(line)
        elif line.startswith('#CHROM'):
            index = 0
            for sample in line.rstrip('\n').split('\t'):
                if index < 9:
                    index += 1
                else:
                    if sample in prochi_linker_dict['NewProchi_GD'] and prochi_linker_dict.get('NewProchi_GD')[sample] is not nan:
                        line = re.sub(sample, prochi_linker_dict.get('NewProchi_GD')[sample], line)
                        corresponding_samples_dict[index] = sample
                        index += 1
                    else:
                        line = re.sub('\t'+sample, '', line)
                        non_corresponding_samples_dict[index] = sample
                        index += 1
            vcf_header.append(line)
            new_vcf_writer.write(line)
        else:
            line = line.split("\t")
            corresponding_sample_vcf_line = [col for col in line if line.index(col) not in non_corresponding_samples_dict.keys()]
            new_vcf_writer.write("\t".join(corresponding_sample_vcf_line))
            if len(non_corresponding_samples_dict.keys()) > 0:
                if path.exists(non_corresponding_vcf_writer_name)==False:
                    vcf_header[-1] = vcf_header[-1].split("\t")
                    vcf_header[-1] = vcf_header[-1][0:9]
                    vcf_header[-1] += non_corresponding_samples_dict.values()
                    vcf_header[-1] = "\t".join(vcf_header[-1]) + '\n'
                    non_corresponding_vcf_writer = open(non_corresponding_vcf_writer_name, 'w')
                    non_corresponding_vcf_writer.writelines(list(vcf_header))
                non_corresponding_sample_line = [col for col in line if line.index(col) < 9 or line.index(col) in non_corresponding_samples_dict.keys()]
                non_corresponding_vcf_writer.write("\t".join(non_corresponding_sample_line))
    new_vcf_writer.close()
    subprocess.run(["bgzip", "{0}.vcf".format(new_vcf_writer_filepath)], cwd = path.dirname(new_vcf_writer_filepath))
    subprocess.run(["tabix", "-p", "vcf", "{0}.vcf.gz".format(new_vcf_writer_filepath)], cwd = path.dirname(new_vcf_writer_filepath))
    if path.exists(non_corresponding_vcf_writer_name)== True:
        non_corresponding_vcf_writer.close()
        subprocess.run(["bgzip", "{0}.vcf".format(non_corresponding_vcf_writer_name)], cwd= new_vcf_writer_filepath)
        subprocess.run(["tabix", "-p", "vcf", "{0}.vcf.gz".format(non_corresponding_vcf_writer_name)], cwd= new_vcf_writer_filepath)
    return None


if __name__=="__main__":
    GoDARTSSampleNameUpdater = argparse.ArgumentParser(prog="GoDARTSSampleNameUpdater")
    GoDARTSSampleNameUpdater.description = "converts the format of the older .vcf.gz files from the GoDARTS project into the newer format."

    GoDARTSSampleNameUpdater.add_argument("-linker", "--godarts_prochi_linker", type=str, 
    required= True, help= "Please indicate the file path to the location of the GoDARTSlinker_PROCHI.csv file.")

    GoDARTSSampleNameUpdater.add_argument("-i", "--input_old_vcf_file", type=str,
    required= True, help= "Please indicate the file path to the location of the older .vcf.gz file.")

    GoDARTSSampleNameUpdater.add_argument("-o", "--output_new_vcf_file", type=str, nargs="?", const= "default",
    required= True, help= "Please indicate the full path of the new .vcf.gz file you wish to create. By default, the output directory will be where the input file is")

    args= GoDARTSSampleNameUpdater.parse_args()

    old_vcf_filepath= args.input_old_vcf_file
    gd_linker= args.godarts_prochi_linker
    new_vcf_filepath = args.output_new_vcf_file

    if old_vcf_filepath.endswith('vcf.gz'):
        if new_vcf_filepath == "default":
            new_vcf_filepath = glob.glob(old_vcf_filepath)[0].rstrip(".vcf.gz")+ "_updated_sample_names.vcf.gz"
        else:
            new_vcf_filepath = new_vcf_filepath
        with gzopen(old_vcf_filepath, 'rt') as old_vcf_reader:
            update_prochi_sample_names(old_vcf_reader, gd_linker, new_vcf_filepath)
    else:
        warnings.warn("Please give a file with the file extension .vcf.gz.")


