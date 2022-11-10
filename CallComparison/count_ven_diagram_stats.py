import os
import glob

def count_lines(file_path):
    return sum(1 for line in open(file_path,'r'))

def count_lines_ignore_double_GATKs(file_path):
    # only count the line if there is not a GATK overlapped with another GATK
    return sum( line.count('GATK') < 2 for line in open(file_path,'r'))

"""
expand('work/ThreeWay/{sample}.3way.{cnv_type}.filtered.bed',sample=SAMPLES,cnv_type=TYPES),
expand('work/ThreeWayBreakpoint/{sample}.3way.{cnv_type}.filtered.bedpe',sample=SAMPLES,cnv_type=TYPES),
expand('work/{sample}.{caller}.{cnv_type}.filtered.bed',sample=SAMPLES,caller=CALLERS,cnv_type=TYPES),
expand('work/BreakPoint/{slop}/{sample}.savvy_x_gatk.{cnv_type}.breakpoint.filtered.bedpe',slop='5000',sample=SAMPLES,cnv_type=TYPES),
expand('work/BreakPoint/{slop}/{sample}.savvy_x_cnvkit.{cnv_type}.breakpoint.filtered.bedpe',slop='5000',sample=SAMPLES,cnv_type=TYPES),
expand('work/BreakPoint/{slop}/{sample}.gatk_x_cnvkit.{cnv_type}.breakpoint.filtered.bedpe',slop='5000',sample=SAMPLES,cnv_type=TYPES),
expand('work/Reciprocal/Results/savvy_x_gatk.{cnv_type}.reciprocal.filtered.bed',cnv_type=TYPES),
expand('work/Reciprocal/Results/savvy_x_cnvkit.{cnv_type}.reciprocal.filtered.bed',cnv_type=TYPES),
expand('work/Reciprocal/Results/gatk_x_cnvkit.{cnv_type}.reciprocal.filtered.bed',cnv_type=TYPES)
"""
# BED
# DUPLICATIONS
# count triples
def count_beds(cnv_type):
    files = glob.glob('work/ThreeWay/*.3way.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
    triple_count = sum(count_lines(f) for f in files)
    print('triples',triple_count)
    # count savvy x cnvkit
    cross_stats = {}
    for cross in ['savvy_x_cnvkit','savvy_x_gatk','gatk_x_cnvkit']:
        files = glob.glob('work/Reciprocal/Results/{cross}.{cnv_type}.reciprocal.filtered.bed'.format(cross=cross, cnv_type=cnv_type))
        cross_stats[cross] = sum(count_lines(f) for f in files) - triple_count
        print(cross, cross_stats[cross])
    # singles
    # gatk
    files = glob.glob('work/*.gatk.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
    gatk_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['savvy_x_gatk'] - cross_stats['gatk_x_cnvkit']
    print('GATK',gatk_count)
    # savvy
    files = glob.glob('work/*.savvy.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
    savvy_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['savvy_x_gatk'] - cross_stats['savvy_x_cnvkit']
    print('savvy',savvy_count)
    # cnvkit
    files = glob.glob('work/*.cnvkit.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
    cnvkit_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['gatk_x_cnvkit'] - cross_stats['savvy_x_cnvkit']
    print('CNVkit',cnvkit_count)

print('RECIP')
print('DUP')
count_beds('DUP')
print()
print('DEL')
count_beds('DEL')

def count_bedpe(cnv_type):
    # expand('work/ThreeWayBreakpoint/{sample}.3way.{cnv_type}.filtered.bedpe',sample=SAMPLES,cnv_type=TYPES),
    files = glob.glob('work/BufferThreeWay/*.3way.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
    triple_count = sum(count_lines_ignore_double_GATKs(f) for f in files)
    print('triples',triple_count)
    # count savvy x cnvkit
    cross_stats = {}
    algos = ['Savvy','CNVkit','GATK']
    for algo1 in algos:
        files = glob.glob('work/BufferDoublesWay/*.2way.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
        for algo2 in algos:
            if algo1 == algo2:
                continue
            cross = algo1 + '-' + algo2
            if cross == 'Savvy-CNVkit':
                print('Debugging', sum([ sum([algo1 in l and algo2 in l for l in open(f,'r')]) for f in files]))
            cross_stats[cross] = sum([ sum([algo1 in l and algo2 in l for l in open(f,'r')]) for f in files]) - triple_count
            print(cross, cross_stats[cross])
    # singles
    # work/BEDPE/{sample}.{caller}.{cnv_type}.filtered.bedpe
    # gatk
    files = glob.glob('work/BEDPE/*.gatk.{cnv_type}.filtered.bedpe'.format(cnv_type=cnv_type))
    gatk_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['Savvy-GATK'] - cross_stats['GATK-CNVkit']
    print('GATK',gatk_count)
    # savvy
    files = glob.glob('work/BEDPE/*.savvy.{cnv_type}.filtered.bedpe'.format(cnv_type=cnv_type))
    savvy_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['Savvy-GATK'] - cross_stats['Savvy-CNVkit']
    print('Savvy',savvy_count)
    # cnvkit
    files = glob.glob('work/BEDPE/*.cnvkit.{cnv_type}.filtered.bedpe'.format(cnv_type=cnv_type))
    cnvkit_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['GATK-CNVkit'] - cross_stats['Savvy-CNVkit']
    print('CNVkit',cnvkit_count)

def both_ol_types(cnv_type):
    #work/BothOverlaps/{sample}.3way.{cnv_type}.filtered.bed
    #work/BothOverlaps/{sample}.2way.{cnv_type}.filtered.bed
    files = glob.glob('work/BothOverlaps/*.3way.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
    triple_count = sum(count_lines(f) for f in files)
    print('Triple',triple_count)
    cross_stats = {}
    algos = ['Savvy','CNVkit','GATK']
    for algo1 in algos:
        files = glob.glob('work/BothOverlaps/*.2way.{cnv_type}.filtered.bed'.format(cnv_type=cnv_type))
        for algo2 in algos:
            if algo1 == algo2:
                continue
            cross = algo1 + '-' + algo2
            cross_stats[cross] = sum([ sum([algo1 in l and algo2 in l for l in open(f,'r')]) for f in files]) - triple_count
            print(cross, cross_stats[cross])
    # singles
    # work/BEDPE/{sample}.{caller}.{cnv_type}.filtered.bedpe
    # gatk
    files = glob.glob('work/BEDPE/*.gatk.{cnv_type}.filtered.bedpe'.format(cnv_type=cnv_type))
    gatk_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['Savvy-GATK'] - cross_stats['GATK-CNVkit']
    print('GATK',gatk_count)
    # savvy
    files = glob.glob('work/BEDPE/*.savvy.{cnv_type}.filtered.bedpe'.format(cnv_type=cnv_type))
    savvy_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['Savvy-GATK'] - cross_stats['Savvy-CNVkit']
    print('Savvy',savvy_count)
    # cnvkit
    files = glob.glob('work/BEDPE/*.cnvkit.{cnv_type}.filtered.bedpe'.format(cnv_type=cnv_type))
    cnvkit_count = sum(count_lines(f) for f in files) - triple_count - cross_stats['GATK-CNVkit'] - cross_stats['Savvy-CNVkit']
    print('CNVkit',cnvkit_count)


print()
print('Buffer')
print('DUP')
count_bedpe('DUP')
print()
print('DEL')
count_bedpe('DEL')
print()
print('Both Recip & Buffer')
print('DUP')
both_ol_types('DUP')
print()
print('DEL')
both_ol_types('DEL')
