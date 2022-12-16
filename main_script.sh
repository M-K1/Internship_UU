#!/usr/bin/env bash
################################################################
### Author: Matthijs Kon                                    ####
### Wrapper script for kmerGWAS, based on Justin's scripts. ####
### Date: 01-12-2022                                        ####
################################################################

#Functions

function show_usage (){
##############################################################################################################
### Function which prints this help message when an option given is not a possibility and exits the program###
##############################################################################################################

    printf -- "\n"
    printf -- "Usage: $0 [-options [parameter]]\n"
    printf -- "\n"
    printf -- "Required options:\n"
    printf -- "--steps/-s\t\t\t<1,2,3,4 or all>\tflag for running some or all sub-steps [default: all][no spaces!!!]\n"
    printf -- "--pfile/-pf\t\t\t<filename>\t\tfile containing the samples and their phenotypes\n"
    printf -- "--outdir/-o\t\t\t<directory>\t\tdirectory to which all output is written\n"
    printf -- "--reads/-r\t\t\t<directory>\t\tdirectory containing paired-end read files [necessary option for step 1 and 4]\n"
    printf -- "--kmers/-kd\t\t\t<directory>\t\tdirectory containing the sampled kmers [necessary option when a non-default directory is used for storing the results]\n"
    printf -- "--heap/-hd\t\t\t<directory>\t\tdirectory containing the kmer heap and table [necessary option when a non-default directory is used for storing the results]\n"
    printf -- "--gwas/-g\t\t\t<directory>\t\tdirectory containing the gwas output [necessary option when a non-default directory is used for storing the results]\n"
    printf -- "--phenotraits/-pt\t\t<filename>\t\tfile containing a list of phenotypes for analysis\n"
    printf -- "--voichek/-v\t\t\t<directory>\t\tdirectory containing the original kmerGWAS installation and additional scripts\n"
    printf -- "--justin/-j\t\t\t<directory>\t\tdirectory containing scripts for an improved kmerGWAS\n"
    printf -- "--blast_db/-bdb\t\t\t<prefix>\t\tprefix of the blastn database for comparing assembled contigs [necessary option for step 4]\n"
    printf -- "--pangenome/-pg\t\t\t<prefix>\t\tprefix of the pangenome bwa and minimap2 index for comparing assembled contigs [necessary option for step 4]\n"


    printf -- "\nOptional options:\n"
    printf -- "--pheno_id/-pid\t\t\t<string>\t\tstring containing the header of the column in which the samples are stored\n"
    printf -- "--omitid/-oid\t\t\t<filename>\t\tfile with Sample ID's to omit during analysis\n"
    printf -- "--klength/-kl\t\t\t<int>\t\t\tkmer length (must be odd) [default: 31]\n"
    printf -- "--nCanon/-nc\t\t\t<int>\t\t\tminimal threshold of kmer-counts to filter kmers for canonization[default: 2]\n"
    printf -- "--threads/-t\t\t\t<int>\t\t\tamount of threads to use for analysis [default: 20]\n"
    printf -- "--minorallelecount/-mac\t\t<int>\t\t\tminor allele count for filtering during heaping of kmers [default: 5]\n"
    printf -- "--mcso/-m\t\t\t<float>\t\t\tminimal coverage of specific orientation of kmer over all kmer presences for filtering during heaping of kmers [default: 0.2]\n"
    printf -- "--minorallelefrequency/-maf\t<float>\t\t\tthreshold minor allele frequency [default: 0.05]\n"
    printf -- "--nBest/-b\t\t\t<int>\t\t\tmaximum amount of kmers to return per phenotype [default: 10001]\n"
    printf -- "--help/-h\t\t\t\t\t\tshow this message\n"
    printf -- "\nDISCLAIMER\n"
    printf -- "Step 1 and 2 haven't been implemented since these steps are only performed once for an entire dataset"
    }

function parse_args () {
##################################################################
###Parse the arguments and save them as variables for later use###
##################################################################
    ###Parse all arguments given
    if [ -z "$1" ]; then
        show_usage
        exit 1
    fi
    while [ ! -z "$1" ]; do
        case "$1" in
            --pfile|-pf)
                ###
                shift
                pheno_file=$1
                ;;
            --steps|-s)
                ### Which steps to include
                shift
                steps=$1
                if [[ $steps =~ "all" ]]; then
                    steps="1,2,3,4"
                fi
                ;;
            --outdir|-o)
                ### Output directory
                shift
                out_dir=$1
                out_dir=${out_dir::-1}
#                 mkdir -p $out_dir
                ;;
            --reads|-r)
                ### Where the reads are stored
                shift
                reads_dir=$1
                reads_dir=${reads_dir::-1}
                ;;
            --kmers|-kd)
                ### Output directory of created kmers
                shift
                k_dir=$1
                k_dir=${k_dir::-1}
                ### Create dirs
                ;;
            --heap|-hd)
                ### Output directory heaping/preparing kmers
                shift
                heap_dir=$1
                heap_dir=${heap_dir::-1}
                ;;
            --gwas|-g)
                ### Raw output directory of kmer GWAS
                shift
                gwas_dir=$1
                gwas_dir=${gwas_dir::-1}
                ;;
            --phenoid|-pid)
                shift
                pheno_id=$1
                ;;
            --phenotraits|-pt)
                ### File containing a list of the traits that are to  be analysed (the header names in the phenotype file)
                shift
                pheno_traits=$1
                ;;
            --omitid|-oid)
                ### Omitted IDs for analysis
                shift
                omit_id=$1
                ;;
            --klength|-k)
                ### Kmer length
                shift
                k_length=$1
                ;;
            --nCanon|-nc)
                ### Minimal threshold for kmer-count
                shift
                nCanon=$1
                ;;
            --threads|-t)
                shift
                threads=$1
                ### Amount of threads
                ;;
            --minorallelecount|-mac)
                ### Minor allele count threshold
                shift
                mac=$1
                ;;
            --mcso|-m)
                ### Minimal coverage of specific orientation of kmer over all kmer presences
                shift
                mcso=$1
                ;;
            --minorallelefrequency|-maf)
                ### Threshold minor allele frequency
                shift
                maf=$1
                ;;
            --nBest|-b)
                ### Maximum amount of kmers to return per phenotype
                shift
                nBest=$1
                ;;
            --voichek|-v)
                ### Directory containing the original kmerGWAS installation and additional scripts
                shift
                voichek=$1
                voichek=${voichek::-1}
                ;;
            --justin|-j)
                ### Directory containing scripts for an improved kmerGWAS
                shift
                justin=$1
                justin=${justin::-1}
                ;;
            --blast_db|-bdb)
                ### Output directory
                shift
                blast_db=$1
                ;;
            --pangenome|-pg)
                ### Output directory
                shift
                pangenome=$1
                ;;
            --help|-h)
                show_usage
                exit 1
                ;;
            *)
                show_usage
                unfound_option=$1
                echo "Error: $unfound_option is not an option"
                exit 1
                ;;
        esac
    shift
    done
    }

function set_values (){
##########################################
###Set default values for unset options###
##########################################
    if [[ "$steps" =~ "1" ]]; then
        if [[ -z "$k_dir" ]]; then
            k_dir=${k_dir:-"$out_dir/kmers"}
            mkdir -p $k_dir/samples
            mkdir -p $k_dir/logs
        else
            mkdir -p $k_dir/samples
            mkdir -p $k_dir/logs
        fi
    fi

    if [[ "$steps" =~ "2" ]]; then
        if [[ -z "$heap_dir" ]]; then
            heap_dir=${heap_dir:-"$out_dir/heap"}
            mkdir -p $heap_dir/logs
        else
            mkdir -p $heap_dir/logs
        fi
    fi


    if [[ "$steps" =~ "3" ]]; then
        if [[ -z "$gwas_dir" ]]; then
            gwas_dir=${gwas_dir:-"$out_dir/gwas"}
            mkdir -p $gwas_dir
        else
            mkdir -p $gwas_dir
        fi
    fi

    if [[ "$steps" =~ "4" ]]; then
        if [[ -z "$gwas_dir" ]]; then
            gwas_dir=${gwas_dir:-"$out_dir/gwas"}
        fi
        if [[ -z "$heap_dir" ]]; then
            heap_dir=${heap_dir:-"$out_dir/heap"}
        fi
    fi
    omit_id=${omit_id:-''}
    k_length=${k_length:-31}
    nCanon=${nCanon:-2}
    threads=${threads:-20}
    mac=${mac:-5}
    mcso=${mcso:-0.2}
    maf=${maf:-0.05}
    nBest=${nBest:-10001}
    }

function check_empty_values () {
#################################################################################################
###Check if all required options are given (accounting for the step that is run). If not, exit###
#################################################################################################
    if [[ -z $pheno_file ]] || [[ -z $steps ]] || [[ -z $out_dir ]] || [[ -z $reads_dir && ("$steps" =~ 1 || "$steps" =~ 4) ]] || [[ -z $pheno_traits && ("$steps" =~ 3 || "$steps" =~ 4) ]] || [[ -z $voichek ]] || [[ -z $justin ]] || \
    [[ -z $blast_db && "$steps" =~ 4 ]] || [[ -z $pangenome && "$steps" =~ 4 ]]; then
        echo "At least 1 required parameter is empty, exiting"
        show_usage
        exit 1
    fi
    }

#############################################################################################################################################################
### From this point forward, the functions called are written by Justin (previous student). I will try to interpret what they do and edit them accordingly###
#############################################################################################################################################################



function cut_pheno_file () {
#######################################################################
###Cut the first column/column of choice to get all the sample names###
#######################################################################
    if [[ $pheno_id == '' ]]; then
        id_lst=$(cut -f 1 < $pheno_file)
    else
        id_lst=$(bash $justin/cut-tsv.sh $pheno_file $pheno_id)
    fi
    }


function sample-kmers () {
################################################################
###Generate (canonised/non-canonised) kmers from sample reads###
################################################################
    id_lst=$(echo $id_lst | cut -d ' ' -f 2-)
    echo "Initialising first step"
    # STEP 1. Sample kmers
    ## Create collection of true kmers per sample including their form (canon/non-canon/both)
    sample_index=1
    echo $id_lst
    for sample_id in $id_lst; do
        echo $sample_id
        sample_reads_lst=$(find -L $reads_dir -name "$sample_id*.fq.gz")
        echo $sample_reads_lst
        if [[ $(isin-lst $sample_id "$omit_ids") -eq 0 ]] || [[ "$sample_reads_lst" == '' ]]; then
            echo "Omitting sample #$sample_index : $sample_id"
            echo $sample_id >> $k_dir/logs/omitted.tsv
        else
            echo "Working on sample #$sample_index : $sample_id "

            sample_kmerdir="$k_dir/samples/$sample_id"
            sample_logdir="$k_dir/logs/$sample_id"
            mkdir -p $sample_kmerdir
            mkdir -p $sample_logdir

            fastq_lst_file="$sample_logdir/fastq_files.txt"
            echo $sample_reads_lst | tr ' ' '\n' > $fastq_lst_file


            has_cns=$(find $sample_kmerdir -name kmers_cns*)
            has_ncns=$(find $sample_kmerdir -name kmers_ncns*)
            has_stranded=$(find $sample_kmerdir -name kmers_stranded)

            if [[ "$has_stranded" == '' ]]; then

                if [[ "$has_cns" == '' ]]; then
                    echo-cmd tee-out $sample_logdir/kmers_cns.log \
                        "$voichek/external_programs/kmc_v3 \
                            -t $threads \
                            -k $k_length \
                            -ci $nCanon \
                            @$fastq_lst_file \
                            $sample_kmerdir/kmers_cns \
                            $sample_kmerdir"
                fi

                if [[ "$has_ncns" == '' ]]; then
                    echo-cmd tee-out $sample_logdir/kmers_ncns.log \
                        "$voichek/external_programs/kmc_v3 \
                            -t $threads \
                            -k $k_length \
                            -ci0 \
                            -b \
                            @$fastq_lst_file \
                            $sample_kmerdir/kmers_ncns \
                            $sample_kmerdir"
            fi

                echo-cmd tee-out $sample_logdir/kmers_stranded.log \
                    "$voichek/bin/kmers_add_strand_information \
                        -c $sample_kmerdir/kmers_cns \
                        -n $sample_kmerdir/kmers_ncns \
                        -k $k_length \
                        -o $sample_kmerdir/kmers_stranded"

            fi


            has_cns=$(find $sample_kmerdir -name kmers_cns*)
            has_ncns=$(find $sample_kmerdir -name kmers_ncns*)
#             has_stranded=$(find $sample_kmerdir -name kmers_stranded)

            if [[ ! "$has_stranded" == '' ]]; then
                if [[ ! "$has_cns" ==  '' ]]; then echo-cmd rm $sample_kmerdir/kmers_cns*; fi
                if [[ ! "$has_ncns" == '' ]]; then echo-cmd rm $sample_kmerdir/kmers_ncns*; fi
            fi

        fi

        ((sample_index++))
    done

    }

function heap-kmers () {
################################################################################################################
###Collect and filter best kmers from individual samples, and from it create a heap, table and kinship matrix###
################################################################################################################
    echo $id_lst | cut -d ' ' -f 1 > $heap_dir/logs/omitted.tsv
    id_lst=$(echo $id_lst | cut -d ' ' -f 2-)


    kmer_table="$heap_dir/kmer_table"
    kmer_kinship="$kmer_table.kinship"

    if [[ ! -f $heap_dir/logs/kmer_files.txt ]]; then
        sample_index=1
        for sample_id in $id_lst; do
            sample_kmers_file="$kmers_dir/samples/$sample_id/kmers_stranded"
            if [[ ! -f $sample_kmers_file ]] || [[ $(isin-lst $sample_id "$omit_ids") -eq 0 ]]; then
                echo "Omitting sample #$sample_index : $sample_id"
                echo $sample_id >> $heap_dir/logs/omitted.tsv
            else
                echo "$sample_kmers_file $sample_id" | tr ' ' '\t' >> $heap_dir/logs/kmer_files.txt
            fi
            ((sample_index++))
        done
    fi

    if [[ "$(find $heap_dir -name kmer_heap*)" == '' ]]; then
        echo-cmd \
            "$voichek/bin/list_kmers_found_in_multiple_samples \
                -l $heap_dir/logs/kmer_files.txt \
                -k $k_length \
                --mac $mac \
                -p $mcso \
                -o $heap_dir/kmer_heap"
    fi

    if [[ "$(find $heap_dir -name kmer_table*)" == '' ]]; then
        echo-cmd \
            "$voichek/bin/build_kmers_table \
                -l $heap_dir/logs/kmer_files.txt \
                -k $k_length \
                -a $heap_dir/kmer_heap \
                -o $heap_dir/kmer_table"
    fi

    if [[ "$(find $heap_dir -name kmer_kinship*)" == '' ]]; then
        echo-cmd \
            "$voichek/bin/emma_kinship_kmers \
                -t $heap_dir/kmer_table \
                -k $k_length \
                --maf $maf \
                    > $kmer_kinship"
    fi

    }

function kmer_gwas () {
##################################################################################################################################
###perform kmer-based GWAS (GRAMMAR gamma score; mod. chi-squared test of relation between genotype to (co-variated) phenotype)###
##################################################################################################################################
    for trait in ${traits[@]}; do
        echo "Starting GWAS on trait: '$trait'"
        trait_dir="$gwas_dir/$trait"
        trait_pheno_file="$trait_dir/logs/phenotypes.tsv"
        trait_omit_ids="$trait_dir/logs/omitted.tsv"
        mkdir -p $trait_dir
        mkdir -p $trait_dir/logs

        if [[ $pheno_id == '' ]]; then pheno_id=$(head -n 1 | cut -f 1 < $pheno_file); fi
        trait_phenos_head=$(bash $justin/cut-tsv.sh $pheno_file $pheno_id $trait | head -n 1 | tr ' ' '_')
        echo $trait_phenos_head | tr ' ' '\t'> $trait_pheno_file
        echo $trait_phenos_head | cut -d ' ' -f 1 > $trait_omit_ids

        sample_index=0;sample_amt=0
        for pheno_line in $(bash $justin/cut-tsv.sh $pheno_file $pheno_id $trait | tr ' ' '_' | tr '\t' ';'); do
            if [[ $sample_index -eq 0 ]]; then ((sample_index++)); continue; fi
            pheno_line=$(echo $pheno_line | tr ',' '.')
            sample_id=$(echo $pheno_line | cut -d ';' -f 1)
            sample_pheno=$(echo $pheno_line | cut -d ';' -f 2)
            if [[ $(isin-lst $sample_id $omit_ids) -eq 0 || "${sample_pheno,,}" == "na" || "${sample_pheno,,}" == "n/a" ]]; then
                echo "Omitting sample #$sample_index : $sample_id"
                echo $sample_id >> $trait_omit_ids
            elif [[ $sample_pheno == '' ]]; then
                echo "Omitting sample #$sample_index : $sample_id"
                echo $sample_id >> $trait_omit_ids
            else
                echo $pheno_line | tr ';' '\t' >> $trait_pheno_file
                ((sample_amt++))
            fi
            ((sample_index++))
        done

        trait_phenos=$(tail -n +3 $trait_pheno_file | cut -f 2)
        trait_phenos_unique=$(echo $trait_phenos | tr ' ' '\n' | sort | uniq)

        if [[ $sample_amt -le $(((($sample_index-$(echo $omit_ids | wc -w)))/2)) ]]; then
            err_msg="Skipped gwas of trait $trait: less phenotype values present than half the total amount of samples included"
            echo $err_msg | tee $trait_dir/logs/error.log
            continue
        elif [[ $(echo $trait_phenos_unique | wc -w) -lt 2 ]]; then
            err_msg="Skipped gwas of trait $trait: less than two unique phenotype values present"
            echo $err_msg | tee $trait_dir/logs/error.log
            continue
        elif [[ $(echo $trait_phenos_unique | wc -w) -eq 2 ]]; then
            value_types=$(for trait_val in $trait_phenos_unique; do [[ $trait_val =~ ^[0-9]+([.][0-9]+)*$ ]] && echo 'n' || echo 'x'; done)
            if [[ ! $(echo $value_types | tr 'x' ' ' | wc -w) == $(echo $value_types | wc -w) ]]; then
                echo "Found exactly two unique categorical phenotype values in gwas of trait $trait: converting to binary numerical"
                unique_i=0
                for unique_val in $trait_phenos_unique; do
                    cat $trait_pheno_file | python3 -c "import sys; open('$trait_pheno_file', 'w').write(sys.stdin.read().replace('$unique_val', '$unique_i'))"
                    ((unique_i++))
                done
            fi
        fi
#         elif [[ $(echo $trait_phenos_unique | wc -w) -ge 3 ]]; then
#             value_types=$(for trait_val in $trait_phenos_unique; do [[ $trait_val =~ ^[0-9]+([.][0-9]+)*$ ]] && echo 'n' || echo 'x'; done)
#             if [[ ! $(echo $value_types | tr 'x' ' ' | wc -w) == $(echo $value_types | wc -w) ]]; then
#                 err_msg="Skipped gwas of trait $trait: more than two unique categorical phenotype values present (advise: convert to linear numerical)"
#                 echo $err_msg | tee $trait_dir/logs/error.log
#                 continue
#             fi
#         fi

#         if [[ $(echo $trait_phenos_unique | wc -w) -ge 2 ]]; then
#             value_types=$(for trait_val in $trait_phenos_unique; do [[ $trait_val =~ ^[0-9]+[.][0-9]+$ ]] && echo 'f' || echo 'x'; done)
#             if [[ ! $(echo $value_types | tr 'x' ' ' | wc -w) == $(echo $value_types | wc -w) ]]; then
#                 trait_phenos_unique_counts=$(echo $trait_phenos | tr ' ' '\n' | sort | uniq -c | tr ' ' '_')
#                 poor_pheno_variation=-1
#                 echo "trait phenos unique counts = $trait_phenos_unique_counts"
#                 for unique_pheno_count in $trait_phenos_unique_counts; do
#                     unique_pheno_count=$(echo $unique_pheno_count | tr '_' ' ' | rev | cut -d ' ' -f 2 | rev)
#                     if [[ $unique_pheno_count -le $(($sample_amt/20)) ]]; then
#                         echo "$unique_pheno_count and $sample_amt"
#                         poor_pheno_variation=0
#                         break
#                     fi
#                 done
#                 if [[ $poor_pheno_variation -eq 0 ]]; then
#                     err_msg="Skipped gwas of trait $trait: non-continuous phenotype value category present that covers less than 5% of all included phenotype values"
#                     echo $err_msg | tee $trait_dir/logs/error.log
#                     continue
#                 fi
#             fi
#         fi

        echo-cmd tee-out $trait_dir/logs/gwas.log \
            "python2.7 $voichek/kmers_gwas.py \
                --pheno $trait_pheno_file \
                --kmers_table $heap_dir/kmer_table \
                --gemma_path $voichek/external_programs/gemma_0_98_1 \
                -l $k_length \
                -p $threads \
                -k $nBest \
                --outdir $trait_dir/results"


    done
    }


####################################################################################
### From this point forward, the functions are once again written by Matthijs Kon###
####################################################################################

function parse_traits () {
###########################
###Parse the traits file###
###########################
    traits=()
    while read line || [ -n "$line" ]; do
        traits+=("$line")
    done < $pheno_traits

    }

function create_fasta_sig_kmers () {
##########################################################################################
###Create a fasta file from the significant kmers in the output of the voichek pipeline###
##########################################################################################
    for trait in ${traits[@]}; do
        mkdir -p $gwas_dir/$trait/postprocessing/
        sig_kmer_file="$gwas_dir/$trait/results/kmers/pass_threshold_5per"
        tail -n +2 $sig_kmer_file > $gwas_dir/$trait/results/kmers/pass_threshold_5per_headerless
        #Print the p-value in the header and the kmer as the fasta sequence while removing the underscore and the numbers from the end of the sequence
        awk '{ gsub(/_.*/, "", $2); print ">seq"NR-0"_"$10"\n"$2; }' < $gwas_dir/$trait/results/kmers/pass_threshold_5per_headerless > $gwas_dir/$trait/postprocessing/sig_kmer.fasta
    done

    }

function prepare_manhattan () {
#############################################################################################################################
###Use BWA to map the kmers/reads against the reference genome and extract their position for creation of a manhattan plot###
#############################################################################################################################
    for trait in ${traits[@]}; do
        seqs="$gwas_dir/$trait/postprocessing/sig_${1}.fasta"
        mkdir -p $gwas_dir/$trait/postprocessing/mapping_out/
        mapping_out="$gwas_dir/$trait/postprocessing/mapping_out"
        reference_genome="/home/sarah/sarah2/LettuceKnow/RefGenome_Salinas/GCF_002870075.3_Lsat_Salinas_v8_genomic.fna"
        #map the kmers with a mapping length of at least 31 bp and no mismatches and convert to bam to save some space
        bwa mem -k 31 -T 31 -vv -a -c 5000 $reference_genome $seqs | samtools view -Sb > $mapping_out/mapping_${1}.bam

        #process the mapping file by splitting the columns and writing the result to a csv file
        samtools view -F 4 $mapping_out/mapping_${1}.bam | awk 'BEGIN{print "SEQ,CHR,BP,P"}; { split($3,a,"|"); split($1,b,"_"); print b[1]","a[4]","$4","b[2] }' > $mapping_out/data_for_manhattan_${1}.csv



    done
    }

function manhattan_plot () {
######################################################################
###Create a manhattan plot (PDF) using Rscript                     ###
###This function is both used for the kmer and reads manhattan plot###
######################################################################
    for trait in ${traits[@]}; do
        manhattan_out="$gwas_dir/$trait/postprocessing/manhattan_${trait}_${1}.png"
        mapping_out_dir="$gwas_dir/$trait/postprocessing/mapping_out"
        threshold="$gwas_dir/$trait/results/kmers/threshold_5per"
        Rscript $justin/Manhattan_Plot.R $mapping_out_dir/data_for_manhattan_"${1}".csv $manhattan_out $threshold "$trait"
    done
    }

function sig_kmers_presence_absence () {
################################################################################################################################################
###Extract the significant kmers from the pass_threshold_5per file and create a list for getting the presence/absence patterns in each sample###
################################################################################################################################################
    for trait in ${traits[@]}; do
        mkdir -p $gwas_dir/$trait/postprocessing/read_extraction/sig_kmers_per_sample
        sig_kmer_file_headerless=$gwas_dir/$trait/results/kmers/pass_threshold_5per_headerless
        #Create list of significant kmers for extracting out of the kmer table for a significant kmer absence/presence table
        awk '{ gsub(/_.*/, "", $2); print $2 }' < $sig_kmer_file_headerless > $gwas_dir/$trait/postprocessing/read_extraction/sig_kmers.txt
        awk '{ gsub(/_.*/, "", $2); printf $2"\t"$10"\n" }' $sig_kmer_file_headerless > $gwas_dir/$trait/postprocessing/read_extraction/sig_kmers_and_pvalue.txt
        outfile_filter_kmers="$gwas_dir/$trait/postprocessing/read_extraction/presence_absence_pattern_sig_kmers.txt"

        if [[ ! -e $outfile_filter_kmers ]]; then
            $voichek/bin/filter_kmers -t $heap_dir/kmer_table -k $gwas_dir/$trait/postprocessing/read_extraction/sig_kmers.txt -o $outfile_filter_kmers
        else
            echo "~/"$outfile_filter_kmers "already exists, delete this file if you want to rerun extracting the kmers from the kmer-table again"
        fi

        #Go through the presence absence table and write the kmer sequence to a fasta file with the name of the sample if it is present in the sampl
        out_dir_extract="$gwas_dir/$trait/postprocessing/read_extraction/sig_kmers_per_sample"
        mkdir -p $out_dir_extract
        presence_absence_table="$gwas_dir/$trait/postprocessing/read_extraction/presence_absence_pattern_sig_kmers.txt"
        sig_kmers_and_pvalue="$gwas_dir/$trait/postprocessing/read_extraction/sig_kmers_and_pvalue.txt"
        awk -v o=$out_dir_extract 'BEGIN{FS = OFS = "\t";} FNR==NR{pval[$1]=$2; next;} FNR==1 { for(i=1;i<=NF;i++){ heading[i]=$i; } } FNR!=1 { for(i=2;i<=NF;i++){ if ($i=="1"){ output[i]=output[i]">kmer"FNR-1"|"pval[$1]"\n"$1"\n"; }}} END{for (i in output){print output[i] >> o"/"heading[i]".fasta"; }}' $sig_kmers_and_pvalue $presence_absence_table

    done
    }



function BBduk () {
###########################################################################
###Extract the reads associated with the significant kmers through BBduk###
###########################################################################
    for trait in ${traits[@]}; do
        mkdir -p $gwas_dir/$trait/postprocessing/read_extraction/extracted_reads
        mkdir -p $gwas_dir/$trait/postprocessing/read_extraction/BBduk_stats
        mkdir -p $gwas_dir/$trait/logs/BBduk/
        samples_with_phenotypes=($(awk 'NR>1 { print $1 }' < $gwas_dir/$trait/logs/phenotypes.tsv))
        base_read_processing_dir="$gwas_dir/$trait/postprocessing/read_extraction"
        for basename in "${samples_with_phenotypes[@]}"; do
            bbduk.sh -Xmx10g in1=$reads_dir/${basename}_1.fq.gz in2=$reads_dir/${basename}_2.fq.gz maskmiddle=f rename=t \
            outm=$base_read_processing_dir/extracted_reads/${basename}_assoc_reads.fq \
            ref=$base_read_processing_dir/sig_kmers_per_sample/${basename}.fasta k=31 hdist=0 \
            stats=$gwas_dir/$trait//postprocessing/read_extraction/BBduk_stats/${basename}_stats.txt 2> $gwas_dir/$trait/logs/BBduk/${basename}_log.txt
        done
    done
    }

# function parse_reads_into_fasta () {
# ###############################################################################
# ###Parse the reads into a fasta format that can be used by prepare_manhattan###
# ###############################################################################
#     for trait in ${traits[@]}; do
#         extracted_read_dir="$gwas_dir/$trait/postprocessing/read_extraction/extracted_reads"
#         mkdir -p $gwas_dir/$trait/postprocessing/read_extraction/reads_fasta_format/
#         samples_with_phenotypes=($(awk 'NR>1 { print $1 }' < $gwas_dir/$trait/logs/phenotypes.tsv))
#         for basename in "${samples_with_phenotypes[@]}"; do
#             #This awk function checks if the line is a header of a read and has a hit with a kmer (to circumvent the paired read without significant value from also being mapped) \
#             #By having 3 seperate field seperators (-F parameter), the p-values are isolated and the lowest is attached to the header in the same way as with the kmer sequenes \
#             #which allows for the same functions to be used
#             awk -F "[\t|=]" '{if(NR%4==1 && $0 ~ /kmer/){m=$0;for(i=1;i<=NF;i++)if($i<m)m=$i; printf(">"substr($1,2)"_"m"\n");getline;print}}' $extracted_read_dir/${basename}_assoc_reads.fq \
#             > $gwas_dir/$trait/postprocessing/read_extraction/reads_fasta_format/${basename}.fasta
#             cat $gwas_dir/$trait/postprocessing/read_extraction/reads_fasta_format/*.fasta > $gwas_dir/$trait/postprocessing/sig_reads.fasta
#             sed -i 's/^\([^ _]*\)_/\1/' $gwas_dir/$trait/postprocessing/sig_reads.fasta
#             sed -i 's/ //g' $gwas_dir/$trait/postprocessing/sig_reads.fasta
#         done
#     done
#     }

function assembly_assoc_reads () {
##############################################################################
###Assemble the extracted reads using SPAdes, all associated reads together###
##############################################################################
    for trait in ${traits[@]}; do
        mkdir -p $gwas_dir/$trait/postprocessing/assembly_combined/
        base_read_processing_dir="$gwas_dir/$trait/postprocessing/read_extraction"
        cat $base_read_processing_dir/extracted_reads/* > $base_read_processing_dir/assoc_reads.fq
        spades.py --12 $base_read_processing_dir/extracted_reads/assoc_reads.fq -o $gwas_dir/$trait/postprocessing/assembly_combined/
    done
    }


function assembly_per_sample () {
########################################################################
###Assemble the extracted reads using SPAdes seperate for each sample###
########################################################################
    for trait in ${traits[@]}; do
        mkdir -p $gwas_dir/$trait/postprocessing/assembly_separate/
        samples_with_phenotypes=($(awk 'NR>1 { print $1 }' < $gwas_dir/$trait/logs/phenotypes.tsv))
        base_read_processing_dir="$gwas_dir/$trait/postprocessing/read_extraction"
        for basename in "${samples_with_phenotypes[@]}"; do
            spades.py --12 $base_read_processing_dir/extracted_reads/${basename}_assoc_reads.fq -o $gwas_dir/$trait/postprocessing/assembly_separate/$basename/
        done
    done

    }

function mapping_pangenome () {
#######################################################################
###Map the created contigs and the associated reads to the pangenome###
#######################################################################
    for trait in ${traits[@]}; do
        mkdir -p $gwas_dir/$trait/postprocessing/mapping_out_pangenome/
        minimap2 -a $pangenome $gwas_dir/$trait/postprocessing/assembly_combined/contigs.fasta | samtools sort -o $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_pangenome_contigs_all_reads.bam
        bwa mem -k 31 -T 31 -vv -a -c 5000 -t 10 $pangenome $gwas_dir/$trait/postprocessing/read_extraction/assoc_reads.fq | samtools view -Sb > $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_assoc_reads_pangenome.bam
        #sort and index reads
        samtools sort $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_assoc_reads_pangenome.bam > $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_assoc_reads_pangenome_sorted.bam
        samtools index $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_assoc_reads_pangenome_sorted.bam
        #sort and index contigs
        samtools sort $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_pangenome_contigs_all_reads.bam > $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_pangenome_contigs_all_reads_sorted.bam
        samtools index $gwas_dir/$trait/postprocessing/mapping_out_pangenome/mapping_pangenome_contigs_all_reads_sorted.bam
    done

    }

function BLASTn () {
##################################################
###Use the contigs created by SPAdes for BLASTn###
##################################################
    for trait in ${traits[@]}; do

    mkdir -p $gwas_dir/$trait/postprocessing/BLASTn/
    seqkit seq -m 500 $gwas_dir/$trait/postprocessing/assembly_combined/contigs.fasta > $gwas_dir/$trait/postprocessing/assembly_combined/contigs_500.fasta
    blastn -db $blast_db -out $gwas_dir/$trait/postprocessing/BLASTn/assembly_ombined.tsv -outfmt 6 -query $gwas_dir/$trait/postprocessing/assembly_combined/contigs_500.fasta

        done
    done
    }







#######################################################################################################################################################################################################

#main
    # Step 1: Parse arguments and set defaults if not set or null
    echo "Parsing arguments"
    parse_args "$@"
    check_empty_values
    parse_traits
    echo "Setting defaults"
    set_values
    echo "Sourcing helpers.sh"
    source $justin/helpers.sh
#    echo "  Sourcing cut-tsv.sh"
#    source $justin/cut-tsv.sh
    cut_pheno_file $pheno_file $pheno_id
    # Step 2: Perform kmerGWAS steps
    # Step 2.1: Perform kmer sampling
    if [[ "$steps" =~ 1 ]]; then
        echo "Performing kmer sampling"
        sample_kmers
    fi
    # Step 2.2: Perform kmer heaping
    if [[ "$steps" =~ 2 ]]; then
        echo "Performing kmer heaping"
        heap-kmers
    fi
    # Step 2.3: Perform kmerGWAS
    if [[ "$steps" =~ 3 ]]; then
        echo "Performing kmerGWAS"
        kmer_gwas
    fi
    # Step 3 Post processing
    if [[ "$steps" =~ 4 ]]; then
        echo "Performing post processing"
        create_fasta_sig_kmers
        echo "Creating manhattan plot using kmer sequences"
        prepare_manhattan kmer
        manhattan_plot kmer
        echo "Extracting significant kmers for read extraction"
        sig_kmers_presence_absence
        echo "Extracting reads"
        BBduk
#         parse_reads_into_fasta
        echo "Assembling significant reads"
        assembly_assoc_reads
        assembly_per_sample
        echo "Mapping contigs and reads against the pangenome"
        mapping_pangenome
        echo "Perfoming BLASTn"
        BLASTn
    fi
