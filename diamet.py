import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

def run_megahit_and_diamond():
    fastq_file = "undetermined_reads.fastq.gz"
    if not os.path.isfile(fastq_file):
        print(f"Error: {fastq_file} not found in the working directory.")
        exit(1)

    # Run Megahit
    megahit_command = f"megahit -r {fastq_file} -o contigs"
    subprocess.run(megahit_command, shell=True, check=True)

    # Run Diamond blastx for contigs
    output_file = "undetermined_contigs_diamet.tsv"
    diamond_command = "/rv_home/iapich/DiaMet/diamond blastx -d /data/diamond/swissprot.dmnd " \
                      f"-q contigs/final.contigs.fa -o {output_file} -f 6 qseqid qlen length sscinames sskingdoms"
    subprocess.run(diamond_command, shell=True, check=True)

    if os.path.isfile(output_file):
        df = pd.read_csv(output_file, sep='\t', header=None,
                         names=['qseqid', 'qlen', 'length', 'sscinames', 'sskingdoms'])

        df_filtered = df.drop_duplicates(subset='qseqid', keep='first')
        df_filtered.to_csv(output_file, sep='\t', header=None, index=False)

        print(f"Duplicates removed from {output_file}.")

        contigs_directory = "contigs"
        if os.path.isdir(contigs_directory):
            subprocess.run(f"rm -r {contigs_directory}", shell=True)
            print(f"Removed {contigs_directory}.")
        else:
            print(f"Error: {contigs_directory} not found.")

    else:
        print(f"Error: {output_file} not found.")
        exit(1)

    print("Megahit and Diamond blastx execution completed.")

def run_diamond_blastx():
    command = "/rv_home/iapich/DiaMet/diamond blastx -d /data/diamond/swissprot.dmnd -q undetermined_reads.fastq.gz -o undetermined_reads_diamet.tsv -f 6 qseqid qlen length sscinames sskingdoms --ultra-sensitive"
    os.system(command)

def remove_duplicate_rows(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)
    df = df.drop_duplicates(subset=0, keep='first')
    df.to_csv(file_path, sep='\t', index=False, header=False)

def custom_colors(entries):
    color_dict = {
        'Eukaryota': (20 / 255, 54 / 255, 66 / 255),
        'Bacteria': (16 / 255, 138 / 255, 140 / 255),
        'Viruses': (168 / 255, 31 / 255, 27 / 255),
        'Archaea': (236 / 255, 153 / 255, 41 / 255)
    }
    return [color_dict[entry] for entry in entries]

def plot_column_5(file_path, total_sequences, pdf_pages):
    df = pd.read_csv(file_path, sep='\t', header=None)
    counts = df[4].value_counts()
    
    desired_order = ['Eukaryota', 'Bacteria', 'Archaea', 'Viruses']
    
    unique_entries = sorted(counts.index.tolist(), key=lambda x: desired_order.index(x))

    colors = custom_colors(unique_entries)

    plt.figure(figsize=(8, 5))

    plt.gca().set_axisbelow(True)
    plt.grid(which='major', axis='y', linestyle='-', linewidth=1.2, alpha=0.5)
    plt.grid(which='minor', axis='y', linestyle='-', linewidth=0.5, alpha=0.2)

    bars = plt.bar(range(len(unique_entries)), counts[unique_entries], color=colors)
    
    plt.ylabel('Count', fontsize=12)
    plt.yscale('log')
    plt.title('Taxonomic classification of undetermined reads', fontsize=14, fontweight='bold', loc='left', y=1.05)
    
    percentage_classified = len(df[0].unique()) / total_sequences * 100
    subtitle = f'{percentage_classified:.2f}% of undetermined reads could be classified on protein level'
    plt.suptitle(subtitle, y=0.92, x=0.125, fontsize=10, ha='left')
    
    legend_labels = [entry for entry in unique_entries]
    legend_colors = custom_colors(unique_entries)
    legend_patches = [Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]
    
    legend = plt.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')

    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(True)
    plt.gca().spines['left'].set_visible(False)

    plt.ylim(1, counts.max() * 1.2)

    pdf_pages.savefig(bbox_inches='tight')
    
    df[df[4] == 'Viruses'][[3]].value_counts().reset_index().to_csv('undetermined_reads_diamet_viral.csv', index=False, header=['virus_species', 'count'])

if __name__ == "__main__":
    pdf_path = 'undetermined_reads_diamet.pdf'

    with PdfPages(pdf_path) as pdf_pages:
        run_megahit_and_diamond()
        
        run_diamond_blastx()

        result_file = "undetermined_reads_diamet.tsv"
        remove_duplicate_rows(result_file)

        os.system("seqkit fq2fa undetermined_reads.fastq.gz > undetermined_reads.fasta")

        with open('undetermined_reads.fasta') as fasta_file:
            total_sequences = sum(1 for line in fasta_file if line.startswith('>'))

        os.remove('undetermined_reads.fasta')

        plot_column_5(result_file, total_sequences, pdf_pages)

        print(f"Script executed successfully.")
