import os

output_dir = '/Users/pierrj/fungap_runs/guy11_debug/braker_patch_test'

prefix = 'SRR8842990'

get_anno_script = '/Users/pierrj/anaconda3/envs/fungap/bin/getAnnoFastaFromJoingenes.py'
adjusted_assembly = '/Users/pierrj/fungap_runs/guy11/fungap_out/maker_out/masked_assembly.fasta.adjusted'
translation_table = '1'

command_merge_transcripts = '/usr/local/BRAKER/scripts/merge_transcript_sets.pl {} {} > {}'.format(
    os.path.join(output_dir, prefix, 'augustus.hints.gtf'), os.path.join(output_dir, prefix, 'GeneMark-ET', 'genemark.gtf'), os.path.join(output_dir, prefix, 'braker.gtf') 
)
os.system(command_merge_transcripts)

command_gtf2gff = '/usr/local/maker/bin/genemark_gtf2gff3 {} > {}'.format(
    os.path.join(output_dir, prefix, 'braker.gtf'), os.path.join(output_dir, prefix, 'braker.gff3')
)
os.system(command_gtf2gff)

# Change file name
command2 = 'mv {} {}'.format(
    os.path.join(output_dir, prefix, 'braker.gff3'), os.path.join(output_dir, prefix, 'braker.gff3')
)
os.system(command2)



command3 = '{} -g {} -o {} -t {} -3 {}'.format(
    get_anno_script, adjusted_assembly,
    prefix, translation_table, os.path.join(output_dir, prefix, 'braker.gtf')
)
os.system(command3)

command4 = 'mv {} {}'.format(
    '{}.aa'.format(prefix),
    os.path.join(output_dir, prefix, 'braker_{}.faa'.format(prefix))
)
os.system(command4)

command3 = '{} -g {} -o {} -t {} -3 {}'.format(
    get_anno_script, adjusted_assembly,
    os.path.splitext(gff3_braker)[0], translation_table, os.path.join(output_dir, prefix, 'braker.gtf')
)
os.system(command3)

command4 = 'mv {} {}'.format(
    '{}.aa'.format(os.path.splitext(gff3_braker)[0]),
    os.path.join(output_dir, prefix, 'braker_{}.faa'.format(prefix))
)
os.system(command4)