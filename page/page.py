""" PAGE stands for post-alignment genotyping of single cell gene expression 
and genotyping of transcriptomes libraries  """

import chileup
import pysam 
import pandas as pd
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os 


def getHileup(samples, positions, path, run):
    """ 
    Gets a pileup of reads at a particular position from bams output by 10X Cellranger count pipeline.

    Args: 
    Samples -- list of sample names (should match cellranger count output directory names )
    Positions -- a dataframe containing gene, variant, chr, and pos
    Path -- path to the cell ranger count project folder - ends in "/"
    Run -- A name for the whole project or experiment 

    Returns a file called hileup_verbose_dups* with each read as a row showing the allele found there
     """
    for s in samples:
        for p in positions.variant: 
            bampath = path + s + "/outs/possorted_genome_bam.bam"
            bam = pysam.AlignmentFile(bampath, "rb")
            snv = positions[positions.variant == p]
            
            # config for chileup 
            config = chileup.Config(tags=["CB", "UB"], track_read_names=True,
            track_base_qualities=True, track_mapping_qualities=True,
            exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY ,#| pysam.FDUP,
            min_base_quality=10, min_mapping_quality=10)
            
            # do the pileup 
            h = chileup.pileup(bam, snv.chr[0], snv.pos[0], config)
            
            # reformat the results - prepare to export 
            h.bases.decode("UTF-8") # Decode from byte
            tups = list(zip(h.tags,h.bases.decode("UTF-8"), h.read_names, h.bqs, h.mqs))
            bases = pd.DataFrame(tups, columns=['Barcodes','base', 'read_name','base_quality', 'map_quality'])
            bases.Barcodes = bases.Barcodes.str.decode("UTF-8")
            bases.read_name = bases.read_name.str.decode("UTF-8")
            bases[['cellbarcode','umi']] = bases.Barcodes.str.split("/", expand=True) 
            keep0 = bases.loc[bases['cellbarcode']!= "."  ] # remove reads with no cell barcode or umi tags
            keep = keep0.loc[keep0['umi']!= "."  ]
            final = keep[['cellbarcode', 'umi', 'base','read_name','base_quality', 'map_quality']]
            
            # Save output file as csv 
            os.mkdir(path + "genotype")
            filename = (path + "genotype/hileup_verbose_dups_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv")
            final.to_csv(filename, index = False)


def getHileup_indels(samples, positions, path, run):
    """ 
    Identifies reads with an indel. See docs for getHileup. 
    
    """
    for s in samples:
        for p in positions.variant: 
            bampath = path + s + "/outs/possorted_genome_bam.bam"
            bam = pysam.AlignmentFile(bampath, "rb")
            snv = positions[positions.variant == p]

            config = chileup.Config(tags=["CB", "UB"], track_read_names=True,
            track_base_qualities=True, track_mapping_qualities=True,
            exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY ,#| pysam.FDUP,
            min_base_quality=10, min_mapping_quality=10)

            h = chileup.pileup(bam, snv.chr[0], snv.pos[0], config)
            barcodes=[a.decode("UTF-8").split("/")[0] for a in h.tags]
            umi=[a.decode("UTF-8").split("/")[1] for a in h.tags]
            read_names = [a.decode("UTF-8") for a in h.read_names]
            tups = list(zip(barcodes,umi, read_names, h.bqs, h.mqs))
            inserts={k:v for k,v in h.insertions}
            data={}
            for i in range(len(h.tags)):
                if i in inserts:
                    data[i]=inserts[i]
                else:
                    data[i]='0'
            del h


            bases = pd.DataFrame(tups, columns=['cellbarcode','umi','read_name','base_quality', 'map_quality'])

            indels = pd.Series(data)
            bases['base']=indels
    #         keep0 = bases.loc[bases['cellbarcode']!= "."  ] # remove reads with no cell barcode or umi tags
    #         keep = keep0.loc[keep0['umi']!= "."  ]
            filename = (path + "genotype/hileup_verbose_dups_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv")    
            bases[(bases.cellbarcode != ".") & (bases.umi != ".")].to_csv(filename, index = False)
            del indels
            del bases
            del tups
            del inserts

def genotypeAssignment(samples, positions, WT_allele, ALT_allele, path, run):
    """ 
    A simplistic way to assign cell barcode MT or WT status. 

    If any MT reads, then MT. Else, WT. 

    WT_allele and ALT_allele are strings. "C" or "T" 
    """
    for s in samples:
        for p in positions.variant: 
            snv = positions[positions.variant == p]
            dat = pd.read_csv((path + "genotype/hileup_verbose_dups_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv"), header=0)
            
            variant = snv.variant[0]
            gene = snv.gene[0]
            # reformat the dataframe 
            datmerge = dat.groupby(['cellbarcode'])['base'].apply(''.join).reset_index() # Merge all basecalls per cell barcode 
            umis = dat.groupby(['cellbarcode'])['umi'].apply(','.join).reset_index() # Retain the UMIs from each cell barcode in separate table 
            datmerge.merge(umis) # Merge basecalls and UMIs back together 
            datmerge['Assignment'] = "NA"

            # Make WT and MT assignments
            datmerge.Assignment[datmerge.base.apply(lambda x: WT_allele in map(str.upper, x))] = 'WT'
            datmerge.Assignment[datmerge.base.apply(lambda x: ALT_allele in map(str.upper, x))] = 'MT'

            datmerge['cellbarcode'] = datmerge['cellbarcode'].str.split('-').str[0]
            datmerge.to_csv((path + "genotype/assignment_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv"))
            fig = plt.figure()
            dat_sort = datmerge.sort_values('Assignment', ascending=True)
    #         ax = dat_sort['Assignment'].value_counts()
    #                                         .plot(kind='bar',
    #                                     figsize=(8,5),
    #                                     title=( s + " " + gene + "_" + variant))
            fig = sns.countplot(x = 'Assignment',
                  data = datmerge,
                  order = ['WT', 'MT'])
            #fig.ylabel = "Number of Cells"
            fig.set(xlabel=' ', ylabel=' ')
            fig.set_title((s ))#+ " " + gene+ " " + variant ))
#             fig.set_ylim([0,10000])
    #         for tick in ax.get_xticklabels():
    #             tick.set_rotation(0)

            plt.tight_layout()
            fig = plt.gcf()
            fig.savefig(path + "genotype/assignment_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".png")
            print(s)
            print(print(datmerge['Assignment'].value_counts()))
    #         # For plot of base calls before assignment 
    #         fig = plt.figure()
    #         dat_sort = datmerge.sort_values('base', ascending=True)
    #         fig = sns.countplot(x = 'base',
    #               data = datmerge) #,

    #         fig.set(xlabel=' ', ylabel=' ')
    #         fig.set_title((s ))
    #         fig = plt.gcf()

def disambiguateBCassignment(samples, positions, WT_allele, ALT_allele, path, run):
    """ 
    A slightly more complex way to do barcode assignments and deals with ambiguous reads/UMIs. 
    
    Throws away UMIs that have 0.15 > x > 0.85 fractions of reads supporting the assignment. 

    Returns a file with cell barcode assignments and the UMIs supporting the assignment. 
    Also makes a plot showing fraction of MT and WT cells. 

    """
    for s in samples:
        for p in positions.variant: 
                snv = positions[positions.variant == p]
                dat = pd.read_csv((path + "genotype/hileup_verbose_dups_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv"), header=0)

                variant = snv.variant[0]
                gene = snv.gene[0]
                # reformat the dataframe 
                dat.base = dat.base.astype(str)
                dat.base.str.upper()
                merged = dat.groupby(['cellbarcode','umi'])['base'].apply(''.join).reset_index() # Merge all basecalls per cell barcode 
                merged['conflict'] = merged.base.str.contains(WT_allele, regex=False) & merged.base.str.contains(ALT_allele, regex=False)

                merged['wt_count'] = merged.base.str.count(WT_allele)
                merged['mt_count'] = merged.base.str.count(ALT_allele)

                merged['ratio'] = merged['mt_count'] / (merged['mt_count'] + merged['wt_count'])

                merged['umi_assignment'] = "N"

                # Make WT and MT assignment
                merged.umi_assignment[merged.ratio.apply(lambda x: x <= 0.15)] = "R" # Ref
                merged.umi_assignment[merged.ratio.apply(lambda x: x >= 0.85)] = "A" # Alt

                merged.to_csv((path + "genotype/UMIassignment_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv"))

                # plot QC
                df = merged[merged['conflict']==True]
                df['wt_count'] = df.base.str.count(WT_allele)
                df['mt_count'] = df.base.str.count(ALT_allele)
                df['ratio'] = df['mt_count'] / (df['mt_count'] + df['wt_count'])
                fig = plt.figure()
                plt.xlim(xmin=0, xmax = 1)
                df['ratio'].hist(bins=100, grid = False)
                plt.axvline(0.15, color='k', linestyle='dashed', linewidth=1)
                plt.axvline(0.85, color='k', linestyle='dashed', linewidth=1)
                print("In sample " + s + ", " + str(df.umi.nunique()) +" ("+ str(round(100* df.umi.nunique() / merged.umi.nunique(),1)) + "%)" + " out of " + str(merged.umi.nunique()) + " umi have conflicts")
                plt.tight_layout()
                fig = plt.gcf()
                fig.savefig(path + "genotype/disambiguate_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".png")


                barcodeAssignments = merged.groupby(['cellbarcode'])['umi_assignment'].apply(''.join).reset_index() # Merge all basecalls per cell barcode 
                umis =  merged.groupby(['cellbarcode'])['umi'].apply(','.join).reset_index() # Retain the UMIs from each cell barcode in separate table 
                barcodeAssignments = barcodeAssignments.merge(umis) # Merge basecalls and UMIs back together 

                barcodeAssignments['cell_assignment'] = "NA"
                barcodeAssignments['cell_assignment'][barcodeAssignments.umi_assignment.apply(lambda x: "R" in map(str.upper, x))] = 'WT'
                barcodeAssignments['cell_assignment'][barcodeAssignments.umi_assignment.apply(lambda x: "A" in map(str.upper, x))] = 'MT'

                barcodeAssignments['cellbarcode'] = barcodeAssignments['cellbarcode'].str.split('-').str[0]
                barcodeAssignments.to_csv((path + "genotype/BCassignment_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".csv"))
                fig = plt.figure()
                fig = sns.countplot(x = 'cell_assignment',
                      data = barcodeAssignments,
                      order = ['WT', 'MT'])
                #fig.ylabel = "Number of Cells"
                fig.set(xlabel=' ', ylabel=' ')
                fig.set_title((s ))#+ " " + gene+ " " + variant ))
#                 fig.set_ylim([0,5500])
        #         for tick in ax.get_xticklabels():
        #             tick.set_rotation(0)

                #plt.tight_layout()
                plt.tight_layout()
                fig = plt.gcf()
                fig.savefig(path + "genotype/BCassignment_" + run + "_" + s + "_" + snv.gene[0] + "_" + p + ".png")
                print(s)
                print(print(barcodeAssignments['cell_assignment'].value_counts()))


def page(samples, positions, WT_allele, ALT_allele, path, run):
    """ 
    Generates a read pileup and cell barcode assignment for a list of samples and positions. 

    samples is a list of sample names - should correspond with the sample names given to the cellranger 
    count output directory 
    positions is a pandas dataframe containing gene, variant, chr, and pos. You may need to subtract 1 from
    pos (0- or 1-based start -- double check )
    WT_allele is a string. ALT_allele is another string. 
    path is the directory containing subdirectories of 10x cellranger count sample directories. 
    run is a name for the project or experiment (string) 

    Example usage: 
    path = "/ifs/work/leukgen/home/ms4/projects/254/10XscRNAseq/dat/GOT/S172/trimmed_count/"
    samples = ['IDH2i-02-T1_DNMT3A-011-R882_E136','IDH2i-02-T2_DNMT3A-011-R882_E136']

    positions = pd.DataFrame({"gene": [ "DNMT3A"], "variant": [ "R882"], "chr" : ["2"], "pos" : [25457241]})
    run = "S172" 
    page(samples, positions, "C" , "T", path, run) """

    getHileup(samples, positions, path, run)
    disambiguateBCassignment(samples, positions, WT_allele, ALT_allele, path, run)
