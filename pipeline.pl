#!/usr/local/bin/perl
#Author: Tomislav Ilicic
#Oragnisation: Wellcome Trust Sanger Institute
##Group: Teichmann Group
##Date: 26/01/2015
#Description PAPERLINE. Pipeline to work with RNA-sequencing data


use strict;
use warnings;
use Switch;
use Carp;
use English qw(-no_match_vars);
use Env;
use File::Path qw(make_path remove_tree);
use File::Copy;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use lib '/homes/ti1/tools/perl_libs/lib/perl5';
use Config::Simple;
use List::MoreUtils qw(uniq);

#DIRECTORIES
my $ROOT_DIR = "";
my $TEMP_DIR="";
my $REF_DIR = "";
my $TOOLS_MAPPING = "";
my $TOOLS_QUANTIFICATION = "";

#EXTENSIONS
my $EXT_BAM = "";
my $EXT_FASTQ= "";
my $EXT_FASTQ_GZ = "";
my $EXT_SAM = "";
my $EXT_SORTED = "";
my $EXT_COUNTS = "";
my $EXT_METRICS = "";
my $EXT_MERGED = "";
my $EXT_MARKED_DUPL= "";
my $EXT_FASTA = "";
my $EXT_LOG = "";
my $EXT_GTF= "";
my $EXT_GSNAP_SPLICE = "";
my $EXT_GSNAP_IIT = "";
my $EXT_DICT = "";
my $EXT_SUMMARY = "";
my $FORWARD_STRAND = "";
my $REVERSE_STRAND = "";


#BINARY SOFTWARES
#MAPPING
my $GSNAP="";
my $GSNAP_BUILD="";
my $GTF_SPLICE_TOOL = "";
my $GTF_IIT_TOOL = "";
my $HTSEQTOOL="";
my $BOWTIE1_BUILD="";
my $BOWTIE1="";
my $BOWTIE2_BUILD="";
my $BOWTIE2="";
my $BWA_BUILD="";
my $BWA="";
my $STAR_BUILD="";
my $STAR = "";
my $SALMON_BUILD = "";
#QUANTIFICATION
my $TOPHAT="";
my $CUFFLINKS="";

#HELPER TOOLS
my $PICARD_TOOL = "";
my $BAM_2_FASTQ="";
my $SAMTOOLS="";

my $CONFIG_FILE="";

#HELPER
my %available_genomes;
my %available_mapper;
my %available_counter;

#PIPELINE HASH LABELS
my $BULK_NAME_LABEL = "BULK";
my $DATA_TYPE_LABEL = "data_type";
my $SINGLE_END_LABEL = "single-end";
my $PAIRED_END_LABEL = "paired-end";
my $BUILD_DB_LABEL = "reference_genome";
my $COUNT_TOOL_LABEL ="count_tool";
my $MAPPING_ARGS_LABEL="mapping_arguments";
my $COUNT_ARGS_LABEL="counting_arguments";
my $USER_IN_DIR_LABEL = "user_in_dir";
my $USER_OUT_DIR_LABEL = "user_out_dir";
my $JOB_INDEX_LABEL = "job_index";
my $MAPPING_TOOL_LABEL = "mapping_tool";
my $GENOME_LABEL = "genome";
my $GTF_FILE_LABEL = "gtf_file";
my $OVERWRITE_LABEL = "overwrite_files";
my $data_type = "null";
my $CONFIG_FILE_LABEL = "config_file";

#MAIN FUNCTION IS RUN.
#THERE HAPPENS EVERYTHING.
#FIRST IT CHECKS WHAT USER HAS SPECIFIED AND IF HE SPECIFIED IT CORRECTLY
#SECOND IT TESTS IF INPUT FILES ARE AVAILABLE
#THEN IT EXECUTES THE TOOL
#ON EACH EXECUTION, OUTPUTFILES ARE CHECKED IF THEY EXIST AND ARE NOT 0
#IF OUTPUFIL EXIST AT THE BEGINNING ALREADY, EXECUTION WILL NOT HAPPEN

#REGARDING MAPPING:

#TOPHAT IS A SPECIAL CASE AND BWA AS WELL.
#IF TOPHAT IS SELECTED, IN THE BACKGROUND IT SUGGESTS BOWTIE REFERENCE DATABSAES
#TO THE USER. THEREFORE CODE IS A BIT UNCLEAR
#BWA NEEDS AN ADDITIONAL STEP TO GENERATE .sam FILE AFTER MAPPING

#REGARDING QUANTIFICATION:
#HTSEQ AND CUFFLINKS DIFFER IN THEIR OUTPUT.

#MAIN FUNCTION
sub run {
    #ALL USER INFORMATION STORED
    my %parameters = &parse_command_line(\@ARGV);
    #CHECKS IF PARAMETERS ARE CORRECT
    %parameters = &check_parameters(\%parameters);
    &write_log_file(\%parameters);
    
    #CREATE REFERENCE DATABSE
    if (exists $parameters{$BUILD_DB_LABEL}) {
        my $ref_dir = $REF_DIR;
        &create_database(\%parameters);
        
    #READ PREPROCESSING,MAPPING,QUANTIFICATION
    } else {
	my $user_in=$ROOT_DIR.$parameters{$USER_IN_DIR_LABEL};
	my $sample_prefix=$parameters{$USER_IN_DIR_LABEL};
	my $user_out=$TEMP_DIR.$parameters{$USER_OUT_DIR_LABEL};
	my $index = $parameters{$JOB_INDEX_LABEL};       	 
        #INITIAL INPUT FILES
        my $input_dir;
        my @input_files;
        
        #PIPELINE MAIN OUTPUT DIRECTORIES
        my $preprocessing       = $user_out."/preprocessed";
        my $mapping_root        = $user_out."/mapped";
        my $counts_dir          = $user_out."/counts";
        
        #MAPPINGOUTPUT DIRECTORIES
        my $mapping_dir         = $mapping_root."/sam";
        my $sorted_dir          = $mapping_root."/sorted_bam";
        my $mapping_stats       = $mapping_root."/stats";
        my $deduplication_dir   = $mapping_root."/deduplicated_bam";

        #QUANTIFICATION OUTPUT DIRECTORIES
        my $standard_counts = $counts_dir."/standard";
        my $de_counts_dir = $counts_dir."/deduplicated";

        #MAP DATA IF USER SPECIFIED
        if (exists $parameters{$MAPPING_TOOL_LABEL}) {
            
            #INPUT DIRECTORY WHERE FILES MUST BE STORED
            $input_dir = $user_in."/raw";
            
            #LOOK FOR CELL WITH NUMBER $INDEX AND ONE OF THE EXTENSIONS
            #EACH ELEMENT IN ARRAY SHOULD BE ONE CELL WITH ONE OF THE EXTENSIONS. SUPPORTS  MULTIPLE LANES PER CELL AS LONG THE SAME SUFFIX
            my @supported_mapping_extensions = ($EXT_BAM, $EXT_FASTQ, $EXT_FASTQ_GZ);
            @input_files = &grep_input_file($input_dir, $index, \@supported_mapping_extensions);
            
            #IF NO FILES FOUND EXIT SCRIPT
            if (&is_empty_array(\@input_files)) {
                &print_error("No raw input files found in: $input_dir");
            }
            
            #CHECK IF SINGLE OR PAIRED END BASED ON FILE NAME
            &detect_if_paired_end(\@input_files, $index);
            #PREPROCESS INPUT DATA SO THAT IT IS READY FOR MAPPING
            my %ouput_hash = &preprocess(\@input_files, $preprocessing, $index, $parameters{$OVERWRITE_LABEL}, $sample_prefix);
            
            #PREPROCESSED INPUT FILES
            @input_files = @{$ouput_hash{"input_files"}};
            #PREPROCESSED BAM FILES FOR SUBMISSION
            my @concat_bam = @{$ouput_hash{"concat_bam"}};
            
            #ACTUAL MAPPING
            @input_files = &map_data(\@input_files, $mapping_dir, $index, $parameters{$MAPPING_TOOL_LABEL},  $parameters{$MAPPING_ARGS_LABEL} ,$parameters{$GENOME_LABEL}, $parameters{$OVERWRITE_LABEL}, $sample_prefix);
            
            #CONVERSION TO BAM AND SORTING
            @input_files = &sort_bam(\@input_files, $sorted_dir, $index, $parameters{$OVERWRITE_LABEL},  $sample_prefix);
            
            #DEDUPLICATION OF READS
            if (&is_paired_end()) {
                @input_files = (@concat_bam, @input_files);
                @input_files = &remove_duplicates(\@input_files, $deduplication_dir, $parameters{$GENOME_LABEL}, $index, $parameters{$OVERWRITE_LABEL}, $sample_prefix);
            }
            #COLLECT SUMMARY INFORMATION
            #&collect_alignment_summary(\@input_files, $mapping_stats, $parameters{$GENOME_LABEL}, $index, $parameters{$OVERWRITE_LABEL});
        }
        
        #QUANTIFICATION OF GENES/ISOFORMS
        #QUANTIFIES DEDUPLICATED AND NON-DUPLICATED READS
        if (exists $parameters{$COUNT_TOOL_LABEL}) {
    
            
            #INPUT DIRECTORY MUST BE WHERE MAPPED AND SORTED BAM FILES ARE STORED
            $input_dir = $sorted_dir;
            my @supported_mapping_extensions = ($EXT_SORTED);
            #MAPPED BAM FILES MUST EXIST
            @input_files = &grep_input_file($input_dir, $index, \@supported_mapping_extensions);
            if (&is_empty_array(\@input_files)) {
                &print_error("No raw input files found in: $input_dir");
            }
            
            #QUANTIFY INPUT FILES OF NON-DEDUPLCIATED READ COUNBTS
            #GTF FILE MUST BE AVAILABLE IN GENOME DIRECTORY
            @input_files = &quantify_expression(\@input_files, $standard_counts, $index, $parameters{$GENOME_LABEL}, $parameters{$COUNT_TOOL_LABEL}, $parameters{$COUNT_ARGS_LABEL}, $parameters{$OVERWRITE_LABEL}, $sample_prefix);
            
            #IF READS WERE DEDUPLICATED QUANTIFY THE EXPRESSION AGAIN
            if (-d $deduplication_dir) {
                $input_dir = $deduplication_dir;
                my @supported_mapping_extensions = ($EXT_BAM);
                #DEDUPLICATED MAPPED BAM FILES MUST EXIST
                @input_files = &grep_input_file($input_dir, $index, \@supported_mapping_extensions);
                if (&is_empty_array(\@input_files)) {
                    &print_error("No raw input files found in: $input_dir");
                }
                
                #OUTPUT DIR
                #QUANTIFY GENE EXPRESSION USING DEDUPLICATED READS
                @input_files = &quantify_expression(\@input_files, $de_counts_dir, $index, $parameters{$GENOME_LABEL}, $parameters{$COUNT_TOOL_LABEL}, $parameters{$COUNT_ARGS_LABEL}, $parameters{$OVERWRITE_LABEL}, $sample_prefix);
                
            }
            
        }
        &print_info_message("Pipeline finished.");
    }
}

#PREPROCESS INPUT DATA
sub preprocess {
    
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $overwrite = $_[3];
    my $sample_prefix = $_[4];
    
    my $bam_to_fastq_dir = $output_dir."/bam2fastq";
    my $fastq_concat_dir = $output_dir."/concatinated_fastq";
    my $concatinated_bam_dir = $output_dir."/concatinated_bam";
    my $gzip_to_fastq_dir = $output_dir."/unzipped_fastq";
    my $ext = &get_file_ext($input_files[0]);
    &print_info_message("Preprocessign data");
    
    switch ($ext) {
        #BAM TO FASTQ CONVERSION - HAPPENS AUTOMATICALLY AND ONLY IF NECCESSARY
        case ($EXT_BAM)   { &print_info_message("BAM detected"); @input_files = &bam_to_fastq(\@input_files, $bam_to_fastq_dir, $index, $overwrite, $sample_prefix) }
        case ($EXT_FASTQ)   { &print_info_message("FASTQ detected");}
        case ($EXT_FASTQ_GZ)   { &print_info_message("FASTQ.GZ detected");  @input_files = &gzip_to_fastq(\@input_files, $gzip_to_fastq_dir, $index, $overwrite)}
        else            { &print_error("Invalid Extension"); }
    }
    
    #SORT FILES TO KNOW WHICH TO CONCATINATE
    @input_files = sort @input_files;
    #CONCATINATION OF FASTQ FILES - HAPPENS AUTOMATICALLY AND ONLY IF NECCESSARY
    @input_files = &concat_fastq(\@input_files, $fastq_concat_dir, $index, $overwrite, $sample_prefix);
    
    #CONVERT CONCATINATED FASTQ TO BAM FILES
    my @concat_bam = &fastq_to_bam(\@input_files, $concatinated_bam_dir, $index, $overwrite, $sample_prefix);
    
    my %out_hash;
    $out_hash{"input_files"} = [ @input_files ];
    $out_hash{"concat_bam"} = [ @concat_bam ];
    
    return (%out_hash);
}

#COLLECT ALIGNMENT SUMMARY
sub collect_alignment_summary {
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $genome_name = $_[2];
    my $index = $_[3];
    my $overwrite = $_[4];
    my $sample_name = $_[5];
    my @output_files;
    
    my $genome_file = $REF_DIR.$genome_name."/".$genome_name.$EXT_FASTA;
    my $output_file = $output_dir."/".$sample_name."_".$index.$EXT_SUMMARY;
    &print_info_message("Mapping statistics");
    push(@output_files, $output_file);
    
    
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    
    #LOG FILE
    my $log_file = $output_dir."/mapping_statistics_$index"."$EXT_LOG";
    unlink $log_file;
    #ASSUMES THAT FILE IS SORTED
    my @command;
    push(@command, "java -jar");
    push(@command, "$PICARD_TOOL/CollectAlignmentSummaryMetrics.jar");
    push(@command, "INPUT=$input_files[0]");
    push(@command, "OUTPUT=$output_file");
    #push(@command, "REFERENCE_SEQUENCE=$genome_file");
    
    
    #EXECUTE COMMAND
    &print_info_message("Input: ".join(",",@input_files), "\n");
    &exec_cmd($log_file, \@command);
    
    #POST CHECKS
    &post_checks(\@output_files);
    &print_info_message("Done");
    
    
}

#REMOVE READ DUPLICATES FROM ALIGNMENTS
#ACCORDING TO PICARDTOOLS
sub remove_duplicates {
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $genome_name = $_[2];
    my $index = $_[3];
    my $overwrite = $_[4];
    my $sample_name = $_[5];
    my @output_files;
    
    #INPUT FILES
    my $unmapped_bam = $input_files[0];
    my $mapped_sorted_bam = $input_files[1];
    #REFERENCE GENOME FASTA FILE NEEDED
    my $genome_file = $REF_DIR.$genome_name."/".$genome_name.$EXT_FASTA;
    
    
    #OUTPUT FILES
    my $output_file = $output_dir."/".$sample_name."_".$index.$EXT_BAM;
    my $output_metrics_file = $output_dir."/".$sample_name."_".$index.$EXT_METRICS;
    my $output_merged_file = $output_dir."/".$sample_name."_".$index.$EXT_MERGED;
    my $output_marked_duplicates_file = $output_dir."/".$sample_name."_".$index.$EXT_MARKED_DUPL;
    
    my $output_name = &get_file_name($output_file);
    &print_info_message("Duplicate removal");
    push(@output_files, $output_file);
    
    
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    
    #LOG FILE
    my $log_file = $output_dir."/duplicate_removal_$index"."$EXT_LOG";
    unlink $log_file;
    
    my @command;
    #MERGE UNALIGNED AND ALIGNED BAM FILES
    push(@command, "java -jar");
    push(@command, "$PICARD_TOOL/MergeBamAlignment.jar");
    push(@command, "UNMAPPED_BAM=$unmapped_bam");
    push(@command, "ALIGNED_BAM=$mapped_sorted_bam");
    push(@command, "OUTPUT=$output_merged_file");
    push(@command, "REFERENCE_SEQUENCE=$genome_file");
    push(@command, "VALIDATION_STRINGENCY=LENIENT");
    
    if (&is_paired_end()) {
        push(@command, "PE=true");
    } else {
        push(@command, "PE=false");
    }
    
    #EXECUTE COMMAND
    &print_info_message("Input: ".join(",",@input_files), "\n");
    &exec_cmd($log_file, \@command);
    
    #MARK DUPLICATE READS
    @command = ();
    push(@command, "java -jar");
    push(@command, "$PICARD_TOOL/MarkDuplicates.jar");
    push(@command, "INPUT=$output_merged_file");
    push(@command, "OUTPUT=$output_marked_duplicates_file");
    push(@command, "METRICS_FILE=$output_metrics_file");
    push(@command, "VALIDATION_STRINGENCY=LENIENT");
    
    #print @command;
    &exec_cmd($log_file, \@command);
    
    #USE SAMTOOLS TO REMOVE DUPLICATE READS AND CREATE BAM FILES
    @command = ();
    push(@command, "$SAMTOOLS view -F 0x400");
    push(@command, "-h");
    push(@command, "-b $output_marked_duplicates_file");
    &exec_cmd_redirect($log_file, \@command, $output_file);
    
    #POST CHECKS
    &post_checks(\@output_files);
    &print_info_message("Done");
    return (@output_files);
}


#CONVERT FASTQ TO BAM FILE USING PICARD TOOLS
sub fastq_to_bam {
    
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $overwrite = $_[3];
    my @output_files;
    
    #INPUT FILES
    my $input_forward = $input_files[0];
    my $input_reverse = $input_files[1];
    
    #SET OUTPUT FILE
    my $output_file = $output_dir."/raw_".$index.$EXT_BAM;
    my $output_name = &get_file_name($output_file);
    &print_info_message("Fastq2Bam");
    push(@output_files, $output_file);
    
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    
    my @command;
    
    push(@command, "java -jar");
    push(@command, "$PICARD_TOOL/FastqToSam.jar");
    push(@command, "FASTQ=$input_forward");
    #IN CASE OF PAIRED END SEQUENING
    if (&is_paired_end()) {
        push(@command, "FASTQ2=$input_reverse");
    }
    push(@command, "OUTPUT=$output_file");
    push(@command, "SAMPLE_NAME=$output_name");
    
    #LOG FILE
    my $log_file = $output_dir."/fastq_bam_$index"."$EXT_LOG";
    unlink $log_file;
    
    #EXECUTE COMMAND
    &print_info_message("Input: ".join(",",@input_files), "\n");
    &exec_cmd_redirect($log_file, \@command, $output_file);
    
    #POST CHECKS
    &post_checks(\@output_files);
    &print_info_message("Done");
    
    return(@output_files);
    
}


#QUANTIFY EXPRESSION OF GENES/ISOFORMS USING MAPPED SAM FILES
sub quantify_expression {
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $genome = $_[3];
    my $quant_tool = $_[4];
    my $args = $_[5];
    my $overwrite = $_[6];
    my $sample_name = $_[7];
    my @output_files;
    my @args = &get_args($args);
    
    #SET OUTPUT FILE
    my $output_file = $output_dir."/".$sample_name."_".$index.$EXT_COUNTS;
    &print_info_message("Quantification");
    push(@output_files, $output_file);
    
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    
    #ERROR IF MULTIPLE INPUT FILES BUT BULK
    &check_multiple_files_bulk(\@input_files, $index);
    &print_info_message("Input: ".join(",",@input_files), "\n");
    
    #GTF FILE USED FOR COUNTING
    my $gtf_file = $REF_DIR.$genome."/$genome".$EXT_GTF;
    my @command;
    
    #CHECK WHICH QUANTIFICATION TOOL WAS SELECTED
    
    #HTSEQ
    #OUTPUTS gene_id, counts
    if(&string_contains($quant_tool, "htseq")) {
        $HTSEQTOOL = $TOOLS_QUANTIFICATION.$quant_tool.$HTSEQTOOL;
        push(@command, $HTSEQTOOL);
        push(@command, "-f bam");
        push(@command, "-r name");
        push(@command, "-a 0");
        push(@command, "-s no");
        push(@command, $input_files[0]);
        push(@command, $gtf_file);
    #CUFFLINKS
    #OUTPUTS COMPLETE GTF FILES + gene_id, counts
    } elsif (&string_contains($quant_tool, "cufflinks")) {
        $CUFFLINKS = $TOOLS_QUANTIFICATION.$quant_tool.$CUFFLINKS;
        push(@command, $CUFFLINKS);
        push(@command, "-g $gtf_file");
        push(@command, $input_files[0]);
    }
    
    #ADDITIONAL ARGUMENTS PASSED BY USER
    @command = (@command, @args);
    
    #LOG FILE
    my $log_file = $output_dir."/counting_$index"."$EXT_LOG";
    unlink $log_file;
    
    #EXECUTE COMMAND
    my $exit_code = &exec_cmd_redirect($log_file, \@command, $output_file);
    
    #POST CHECKS
    &post_checks(\@output_files);
    if ($exit_code != 0) {
        &print_error("Exit code not 0. Check $log_file for problems.");
    }
    &print_info_message("DONE");
    return(@output_files);
    
}

#SORT BAM FILES ACCORDING TO ITS POSITION USING PICARD TOOLS
sub sort_bam {
    
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $overwrite = $_[3];
    my $sample_name = $_[4];
    my @output_files;
    
    #SET OUTPUT FILE
    my $output_file = $output_dir."/".$sample_name."_".$index.$EXT_SORTED;
    push(@output_files, $output_file);
    &print_info_message("Sam2Bam conversion and positional ordering");
    
    #PRECHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    
    #ERROR IF MULTIPLE INPUT FILES BUT BULK
    &check_multiple_files_bulk(\@input_files, $index);
    my @command;
    push(@command, "java -jar");
    push(@command, "$PICARD_TOOL/SortSam.jar");
    push(@command, "INPUT=$input_files[0]");
    push(@command, "OUTPUT=$output_file");
    push(@command, "SORT_ORDER=queryname");
    push(@command, "CREATE_INDEX=true");
    push(@command, "VALIDATION_STRINGENCY=LENIENT");
    
    #LOG FILE
    my $log_file = $output_dir."/sorting_$index"."$EXT_LOG";
    unlink $log_file;
    
    #EXECUTE COMMAND
    &print_info_message("Sorting of: ".join(",",@input_files), "\n");
    &exec_cmd($log_file, \@command);
    
    #POST CHECKS
    &post_checks(\@output_files);
    &print_info_message("DONE");
    return(@output_files);

}

#MAP PROCESSED RAW READS ONTO GENOME OR TRANSCRIPTOME
#TOPHAT AND BWA ARE SPECIAL CASES AS THEY TREAT DATA DIFFERENTLY
#THEREFORE CODE NOT SUPER NICE
sub map_data {
    
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $mapper = $_[3];
    my @args = &get_args($_[4]);
    my $used_genome = $_[5];
    my $overwrite = $_[6];
    my $sample_name = $_[7];
    my @output_files;
    my $output_file = $output_dir."/$sample_name"."_".$index.$EXT_SAM;
    
    push(@output_files, $output_file);
    &print_info_message("Read mapping");
    
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    my $genome_name = &get_file_name_no_ext($used_genome);
    my $genome_dir = $REF_DIR."/".$genome_name."/".$mapper;
    
    my $file_f;
    my $file_r;
    my $log_file = $output_dir."/mappping_".$index.$EXT_LOG;
    unlink $log_file;
    my @command;
    
    #READ FIRST TWO LINES OF INPUT FILE TO GET THE READ LENGTH
    #IMPORTAT TO KNOW FOR SOME MAPPING ALGORITHMS
    my $read_length = 0;
    open my $info, $input_files[0] or die "Could not open $input_files[0]: $!";
    while( my $line = <$info>)  {
        $read_length = length($line)-1;
        last if $. == 2;
    }
    close $info;
    print "Read length:$read_length\n";
    
    #MAP DATA THROUGH USER SPECIFIED MAPPING TOOL
    
    #GSNAP MAPPING
    if (&string_contains($mapper, "gmap")) {
        $GSNAP = $TOOLS_MAPPING.$mapper.$GSNAP;
        push(@command, $GSNAP);
        push(@command, "-A sam -B 5 -t 10");
        push(@command, "-s $genome_dir/$genome_name");
        
        if (-e $genome_dir) {
            push(@command, "-d $genome_name");
            push(@command, "-D $genome_dir");
        } else {
            &print_error("Database does not exist '$genome_dir'");
        }
        #PAIRED END DATA
        if (&is_paired_end()) {
            $file_f = $input_files[0];
            $file_r = $input_files[1];
            push(@command,$file_f);
            push(@command,$file_r);
            #SINGLE END DATA
        } else {
            $file_f = $input_files[0];
            push(@command,$file_f);
        }
        
    #BOWTIE1, BOTWIE2 MAPPING
    } elsif (&string_contains($mapper, "bowtie-1")) {
            my $BOWTIE = $TOOLS_MAPPING.$mapper.$BOWTIE1;
        #EXPORT LOCATION OF GENOME DATABASE. OTHERWISE IT DOESN"T KNOW WHERE TO LOOK FORE
        push(@command, "export BOWTIE_INDEXES=$genome_dir && ");
        push(@command, $BOWTIE);
        push(@command, "-S");
        push(@command, $genome_dir."/".$genome_name);
        
        #PAIRED END DATA
        if (&is_paired_end()) {
            $file_f = $input_files[0];
            $file_r = $input_files[1];
            push(@command, "-1");
            push(@command,$file_f);
            push(@command, "-2");
            push(@command,$file_r);
            #SINGLE END DATA
        } else {
            $file_f = $input_files[0];
            push(@command,$file_f);
        }
    } elsif (&string_contains($mapper, "bowtie2")) {
        my $BOWTIE = $TOOLS_MAPPING.$mapper.$BOWTIE2;
        push(@command, "export BOWTIE2_INDEXES=$genome_dir && ");
        push(@command, $BOWTIE);
        push(@command, "-x");
        push(@command, $genome_dir."/".$genome_name);
        #PAIRED END DATA
        if (&is_paired_end()) {
            $file_f = $input_files[0];
            $file_r = $input_files[1];
            push(@command, "-1");
            push(@command,$file_f);
            push(@command, "-2");
            push(@command,$file_r);
            #SINGLE END DATA
        } else {
            $file_f = $input_files[0];
            push(@command,$file_f);
        }   
    #BWA1, BW2 MAPPING
    } elsif (&string_contains($mapper, "bwa")) {
        $BWA = $TOOLS_MAPPING.$mapper.$BWA;
        my $output_sai_file = $output_dir."/$sample_name"."_".$index.".sai";
        my @output_sai_files;
        my $exit_code = 0;
        
        #IF READ LENGTH IS SMALLER THAN 100bp, RUN "aln" ALGORITHM FROM BWA
        if ($read_length  < 70) {
            #ALIGN EACH FILE SEPERATLY (output .sai files)
            for (my $i = 0; $i <= $#input_files; $i++) {
                @command = ();
                push(@command, $BWA);
                push(@command, "aln");
                $output_sai_file = $output_dir."/$sample_name"."_".$index."_".($i+1).".sai";
                push(@command, $genome_dir."/".$genome_name);
                
                $file_f = $input_files[$i];
                push(@command,$file_f);
                push(@command, $output_sai_file);
                
                #ADDITIONAL ARGUMENTS PASSED BY USER
                @command = (@command, @args);
                #EXECUTE COMMAND
                $exit_code = &exec_cmd_redirect($log_file, \@command, $output_sai_file);
                push(@output_sai_files, $output_sai_file);
            }
            
            #CHECK IF SAI FILES EXIST
            &post_checks(\@output_sai_files);
            if ($exit_code != 0) {
                &print_error("Exit code not 0. Check $log_file for problems.");
            }
            @command = ();
            push(@command, $BWA);
            #DEPENDING IF PAIRED OR SINGLE, RUN SAMPE, SAMSE TO GENERATE
            #.sam FILE
            if (&is_paired_end()) {
                push(@command, "sampe");
                push(@command, $genome_dir."/".$genome_name);
                push(@command, $output_sai_files[0]);
                push(@command, $output_sai_files[1]);
                push(@command, $input_files[0]);
                push(@command, $input_files[1]);
                #SINGLE END DATA
            } else {
                push(@command, "samse");
                push(@command, $genome_dir."/".$genome_name);
                push(@command, $output_sai_files[0]);
                push(@command, $input_files[0]);
            }
            #SET USER SPECIFIED ARGUMENTS TO NULL, SINCE THEY HAVE BEEN USED ALREADY
            @args = ();
        } else {
            push(@command, $BWA);
            push(@command, "mem");
            push(@command, $genome_dir."/".$genome_name);
            push(@command, $input_files[0]);
            if (&is_paired_end()) {
                push(@command, $input_files[1]);
            }
        }
    #TOPHAT MAPPING
    } elsif (&string_contains($mapper, "tophat")) {
        $TOPHAT = $TOOLS_MAPPING.$mapper.$TOPHAT;
        #TOPHAT USES BOWTIE INDEX AND IS NOT ABLE TO GENERATE INDEX BY ITSELF
        #THEREFORE LOOK UP MOST RECENT DATABASE FOR MAPPING USING BOWTIE2
        $mapper = &get_latest_mapper(\%available_mapper, "bowtie2", "tophat", 1);
        $genome_dir = $REF_DIR."/".$genome_name."/".$mapper;
        push(@command, "export BOWTIE2_INDEXES=$genome_dir/ && ");
        push(@command, "export PATH=\$PATH:$TOOLS_MAPPING$mapper && ");
        push(@command, $TOPHAT);
        @command = (@command, @args);
        @args = ();
        push(@command, $genome_name);
        push(@command, @input_files);

    #STAR MAPPER
    } elsif (&string_contains($mapper, "star")) {
        $STAR = $TOOLS_MAPPING.$mapper.$STAR;
        push(@command, $STAR);
        push(@command, "--genomeDir");
        push(@command, $genome_dir);
        push(@command, "--readFilesIn");
        push(@command,  join(' ', @input_files));
        push(@command, "--genomeLoad");
        push(@command, "LoadAndKeep");
        push(@command, "--outStd");
        push(@command, "SAM");
        my $TEMP = $output_dir."/$sample_name"."_".$index."/";
        remove_tree($TEMP);
        push(@command, "--outTmpDir");
        push(@command, $TEMP);
        my @arg_new;
        foreach $a (@args) {
            push(@arg_new , "-".$a);
        } 
        @args = @arg_new;   
    }
    print @command;
    #ADD ADDITIONAL ARGUMENTS
    @command = (@command, @args);
    #EXECUTE COMMAND
    &print_info_message("Input: ".join(",",@input_files), "\n");
    my $exit_code = &exec_cmd_redirect($log_file, \@command, $output_file);
    
    #POST CHECKS
    &post_checks(\@output_files);
    if ($exit_code != 0) {
        &print_error("Exit code not 0. Check $log_file for problems.");
    }
    &print_info_message("DONE");
    
    return(@output_files);
}

#CREATE DATABASES FROM FASTA FILE
sub create_database {
    
    my %parameters = %{$_[0]};
    #SPECIFIED GENOME AND MAPPER
    my $mapper = $parameters{$MAPPING_TOOL_LABEL};
    my $genome_file = $parameters{$BUILD_DB_LABEL};
    my $genome_dir = "";
    my $genome_name = &get_file_name($genome_file);
    $genome_name =~ s/\./_/g;
    
    #ROOT/ORGANISM_NAME/MAPPER_AND_VERSION/FILES
    my @ar = split("/", $mapper);
    my $version = $ar[$#ar];
    my $organism_dir = $REF_DIR."/".$genome_name;
    $genome_dir = $organism_dir."/".$version;
    
    
    #HERE YOU MAYBE WANT TO DO MORE CHECKS TO VERIFY WHETHER DB EXISTS OR NOT
    if (-d $genome_dir) {
        &print_error("Database already exists in: $genome_dir. Please manually delete or rename input_file");
    } else {
        make_path($genome_dir);
    }
    print $genome_dir;
    
    &print_info_message("Database creation");
    my @command;
    my $log_file = $genome_dir."/".$genome_name."_genome_build".$EXT_LOG;
    unlink $log_file;
    my $gtf_f = $parameters{$GTF_FILE_LABEL};
    
    #CHOOSE PROGRAM FOR DATABASE CREATION
    
    #GSNAP. GSNAP GENERATES SPLICEING FILE AND OTHER FILES THAT IT NEEDS FOR MAPPING. THEREFORE MULTIPLE COMMANDS NEED TO BE EXECUTED
    if (&string_contains($mapper, "gmap")) {
        #FULL PATH TO BINARIES
        $GSNAP_BUILD = $TOOLS_MAPPING.$mapper.$GSNAP_BUILD;
        $GTF_IIT_TOOL = $TOOLS_MAPPING.$mapper.$GTF_IIT_TOOL;
        $GTF_SPLICE_TOOL = $TOOLS_MAPPING.$mapper.$GTF_SPLICE_TOOL;
        
        my $gtf_out = $genome_dir."/".$genome_name."$EXT_GSNAP_SPLICE";
        push(@command, "cat");
        push(@command, $gtf_f);
        push(@command, "| $GTF_SPLICE_TOOL -");
        &print_info_message("Creating splice_sites file $gtf_out");
        &exec_cmd_redirect($log_file, \@command, $gtf_out);
        
        @command = ();
        my $splice_file = $genome_dir."/".$genome_name."$EXT_GSNAP_IIT";
        push(@command, "cat  $gtf_out | ");
        push(@command, $GTF_IIT_TOOL);
        push(@command, "-o $splice_file");
        push(@command, $gtf_out);
        &print_info_message("Creating splice_sites $EXT_GSNAP_IIT file $splice_file");
        &exec_cmd($log_file, \@command);
        
        @command = ();
        push(@command, $GSNAP_BUILD);
        push(@command, "-d $genome_name");
        push(@command, $genome_file);
        push(@command, "-D $genome_dir");
        
    #BOWTIE
    #BOTWIE1 AND BOWTIE2 ARE EXECUTABLE IN THE SAME WAY
    #GENERATED FILES NEED TO BE TRANSFERED MANUALLY
    } elsif (&string_contains($mapper, "bowtie")) {
        @command = ();
        my $BOWTIE_BUILD ="";
        if (&string_contains($mapper, "bowtie-1")) {
            $BOWTIE_BUILD = $TOOLS_MAPPING.$mapper.$BOWTIE1_BUILD;
        } elsif (&string_contains($mapper, "bowtie2")) {
            $BOWTIE_BUILD = $TOOLS_MAPPING.$mapper.$BOWTIE2_BUILD;
        }
        push(@command, "export BOWTIE_INDEXES=$genome_dir && (");
        push(@command, $BOWTIE_BUILD);
        push(@command, $genome_file);
        push(@command, $genome_name);
        push(@command, "; mv $genome_name*.ebwt $genome_name*.bt2 $genome_dir)");
    #BWA
    #GENERATED FILES NEED TO BE TRANSFERED MANUALLY
    } elsif (&string_contains($mapper, "bwa")) {
        $BWA_BUILD = $TOOLS_MAPPING.$mapper.$BWA_BUILD;
        push(@command, "(");
        push(@command, "$BWA_BUILD");
        push(@command, "index");
        push(@command, "-a bwtsw");
        push(@command, "-p $genome_name");
        push(@command, $genome_file);
        push(@command, "; mv $genome_name*.fai $genome_name*.rpac $genome_name*.amb $genome_name*.ann $genome_name*.pac $genome_name*.bwt $genome_name*.rbwt $genome_name*.rsa $genome_name*.sa $genome_dir )");
    } elsif (&string_contains($mapper, "STAR")) {

        $STAR_BUILD = $TOOLS_MAPPING.$mapper.$STAR_BUILD;
        push(@command, "$STAR_BUILD");
        push(@command, "--runMode genomeGenerate");
        push(@command, "--genomeDir");
        push(@command, $genome_dir);
        push(@command, "--genomeFastaFiles");
        push(@command, $genome_file);
        push(@command, "--limitGenomeGenerateRAM 53379893120 ");

    } elsif (&string_contains($mapper, "salmon")) {
	$SALMON_BUILD = $TOOLS_MAPPING.$mapper.$SALMON_BUILD;
        push(@command, "DYLD_LIBRARY_PATH=$TOOLS_MAPPING$mapper/lib && ");
	push(@command, "$SALMON_BUILD");
	push(@command, "index");
	push(@command, "-t");
	push(@command, $genome_file);
	push(@command, "-i");
	push(@command, $genome_dir);
    }

    &print_info_message("Creating index files");
    &exec_cmd($log_file, \@command);
    
    #COPY GTF FILE TO ORGANISM DIR
    my $gtf_copy = $organism_dir."/".$genome_name.$EXT_GTF;
    &print_info_message("Copying gtf file to $gtf_copy");
    
    unless (-e $gtf_copy) {
        copy($gtf_f, $gtf_copy) or die "Copy failed: $!";
        
    }
    #COPY REFERNCE FASTA FILE TO ORGANISM DIR
    my $ref_genome_f_out = $organism_dir."/".$genome_name.$EXT_FASTA;
    &print_info_message("Copying FASTA file to $ref_genome_f_out");
    
    unless (-e $ref_genome_f_out) {
        copy($genome_file, $ref_genome_f_out);
    }
    
    #CREATE SEQEUENCE DICTIONARY
    @command = ();
    push(@command, "java -jar");
    push(@command, $PICARD_TOOL."/CreateSequenceDictionary.jar");
    push(@command, "REFERENCE=$ref_genome_f_out");
    push(@command, "OUTPUT=$ref_genome_f_out".$EXT_DICT);
    &print_info_message("Create Sequencing directory $ref_genome_f_out");
    &exec_cmd($log_file, \@command);
    
}

#BAM TO FASTQ. CONVERTS BAM FILES TO FASTQ FILES.
#WORKS WITH SINGLE CELL AND BULK DATA
#BOTH SINGLE AND PAIRED END
#TEST SINGLE END DATA MISSING
sub bam_to_fastq {
    
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $overwrite = $_[3];
    my $sample_name = $_[4];
    
    #DEFINE OUTPTU_FILES
    my @output_files;
    my $output_file = "";
    my @program_output_files;
    foreach my $file (@input_files) {
        #PAIRED END DATA
        if (index($file, '#') != -1) {
            my @split = split('\#', &get_file_name_no_ext($file));
            #IMPORTANT TO USE ORIGNAL NAME NOT SAMPLE NAME AS OTHERWISE INFORMATION ABOUT LANES GET LOST
            $output_file = $output_dir."/".$split[0];
            my $num = $split[$#split];
            push(@output_files,  $output_file."_".$num."_1".$EXT_FASTQ);
            push(@output_files,  $output_file."_".$num."_2".$EXT_FASTQ);
            $output_file = $output_file."_".$num."#".$EXT_FASTQ;
            push(@program_output_files, $output_file);
            #DEFINING DATA TYPE
            $data_type = $PAIRED_END_LABEL;
            #SINGLE END DATA
        } else {
            $output_file = $output_dir."/".&get_file_name_no_ext($file).$EXT_FASTQ;
            push(@output_files,  $output_file);
            push(@output_files,  $output_file);
            push(@program_output_files, $output_file);
            #DEFINING DATA TPYE
            $data_type = $SINGLE_END_LABEL;
        }
    }
    
    &print_info_message("Bam2Fastq conversion");
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    #ERROR IF MULTIPLE INPUT FILES BUT BULK
    &check_multiple_files_bulk(\@input_files, $index);
    #CONVERT EACH INPUT FILE
    my $log_file = $output_dir."/bam2Fastq_$index"."$EXT_LOG";
    print $log_file;
    unlink $log_file;
    &print_info_message("Input: ".join(",",@input_files), "\n");
    
    for (my $i = 0; $i <= $#program_output_files; $i++) {
        my $file = $input_files[$i];
        my $output_file = $program_output_files[$i];
        my @command;
        push(@command, $BAM_2_FASTQ);
        push(@command,  "-o $output_file");
        push(@command,  $file);
        push(@command, "--force");
        &exec_cmd($log_file, \@command);
    }
    
    #POST CHECKS
    &post_checks(\@output_files);
    &print_info_message("DONE");
    return(@output_files);
}

#UNZIPS GZIP DATA AND GIVES FASTQ FILES BACK
sub gzip_to_fastq {
    
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $overwrite = $_[3];
    my $sample_name = $_[4];
    
    #DEFINE OUTPTU_FILES
    my @output_files;
    my $output_file = "";
    foreach my $file (@input_files) {
        $output_file = $output_dir."/".&get_file_name_no_ext($file).$EXT_FASTQ;
        push (@output_files, $output_file);
    }
    
    &print_info_message("Unzipping files");
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    
    #ERROR IF MULTIPLE INPUT FILES BUT BULK
    &check_multiple_files_bulk(\@input_files, $index);
    
    #UNZIP FILES
    for (my $i=0; $i <= $#input_files; $i++) {
        my $status = gunzip $input_files[$i] => $output_files[$i], AutoClose=>1 or die "gunzip failed: $GunzipError\n";
    }
    #POST CHECKS
    &post_checks(\@output_files);
    &print_info_message("DONE");
    return(@output_files);
}

#CONCATINATE FASTQ FILES
sub concat_fastq {
    my @input_files = @{$_[0]};
    my $output_dir = $_[1];
    my $index = $_[2];
    my $overwrite = $_[3];
    my $sample_name = $_[4];
    my @output_files;
    
    #IF TWO FILES (BULK OR ONLY ONE LANE FOR SINGLE CELL DATA) (_1.fastq, _2.fastq)
    if ($#input_files == 1) {
        return(@input_files);
    }
    
    &check_multiple_files_bulk(\@input_files, $index);
    &print_info_message("Fastq concat");
    #IF MORE THEN TWO FILES (MULTIPLE LANES FROM SINGLE CELL DATA)
    #DEFINE OUTPTU_FILES
    my @in_forward;
    my @in_reverse;
    my $output_file_f = "";
    my $output_file_r = "";
    my $file;
    #PAIRED END DATA
    if (&is_paired_end()) {
        for (my $i= 0; $i <= $#input_files; $i++) {
            
            $file = $input_files[$i];
            if ($file =~ /.*$FORWARD_STRAND$EXT_FASTQ/) {
                push(@in_forward, $file);
            } else {
                push(@in_reverse, $file);
            }
        }
        $output_file_f = $output_dir."/".$sample_name."_".$index.$FORWARD_STRAND.$EXT_FASTQ;
        $output_file_r = $output_dir."/".$sample_name."_".$index.$REVERSE_STRAND.$EXT_FASTQ;
        push(@output_files, $output_file_f);
        push(@output_files, $output_file_r);
        #SINGLE END DATA
    } else {
        @in_forward = @input_files;
        $output_file_f = $output_dir."/".$sample_name."_".$index."_single_end$EXT_FASTQ";
        push(@output_files, $output_file_f);
    }
    
    #PRE CHECKS
    if (!&pre_checks(\@output_files, $output_dir, \@input_files, $index, $overwrite)) {
        return @output_files;
    }
    #ERROR IF MULTIPLE INPUT FILES BUT BULK
    &check_multiple_files_bulk(\@input_files, $index);
    
    #CONCAT FORWARD FILES
    my $log_file = "concat_$index"."$EXT_LOG";
    unlink $log_file;
    &print_info_message("Input: ".join(",",@input_files), "\n");
    &concat_files(\@in_forward, $output_file_f);
    &concat_files(\@in_reverse, $output_file_r);
    
    #CHECK IF OUTPUT EXISTS
    #MAYBE YOU WANT TO DO ADDITIONAL CHECKS HERE OF THE LOG FILES
    &post_checks(\@output_files);
    &print_info_message("DONE");
    return(@output_files);
    
}


########################################################
##############HELPER FUNCTIONS##########################
########################################################

#DETECTS IF PAIRED END DATA
sub detect_if_paired_end {
    my @input_files = @{$_[0]};
    my $ext = "";
    my $index = $_[1];
    $data_type = $SINGLE_END_LABEL;
    my @suffixes;
    foreach my $file (@input_files) {
        $ext = &get_file_ext($file);
        $file = &get_file_name_no_ext($file);
        switch ($ext) {
            
            #IF BAM FILE: FILE NAME NEEDS TO CONTAIN A HASH '#'
            case ($EXT_BAM) {
                #PAIRED END DATA
                if (index($file, '#') != -1) {
                    #DEFINING DATA TYPE
                    $data_type = $PAIRED_END_LABEL;
                }
                
            }
            #IF FASTQ FILE: EITHER _1 OR _2 NEED TO MAP BEFORE THE EXTENSION
            case ($EXT_FASTQ) {
                if ($file =~ /(.*)(_[1-2]$)/g) {
                    push(@suffixes, $2);
                }
            }
            case ($EXT_FASTQ_GZ) {
                if ($file =~ /(.*)(_[1-2]$)/g) {
                    push(@suffixes, $2);
                }
            }
            
        }
    }
    my @unique_suffix = uniq @suffixes;
    if ($#unique_suffix > 0) { 
        my %params = map { $_ => 1 } @unique_suffix;
        if(exists($params{"_1"}) & exists($params{"_2"})) {
            $data_type = $PAIRED_END_LABEL;
        }

    }
    &print_info_message("Data type: " . $data_type);
}

#RETURN IF PAIRED OR NOT BASED ON GLOBAL VARIABLE
sub is_paired_end {
    if ($data_type eq $PAIRED_END_LABEL) {
        return 1;
    } else {
        return 0;
    }
}

#GET USER ARGUMENTS SEPERATE BY A COMMA ','
sub get_args {
    my $args = $_[0];
    
    if (!defined $args) {
        $args = "";
    }
    $args =~ s/=/ /g;
    my @a = split(",", $args);
    
    if (!@a) {
        @a = ();
    } else {
        foreach my $arg (@a) {
            $arg = "-".$arg;
        }
    }
    return(@a);
}

#PRE CHECKS TO SEE IF INPUT FILE EXISTS, OUTPUT EXISTS ETC
sub pre_checks {
    my @output_files = @{$_[0]};
    my $output_dir = $_[1];
    my @input_files = @{$_[2]};
    my $index = $_[3];
    my $overwrite = $_[4];
    
    #CHECK IF OUTPUT FILES EXIST
    my $exists = &files_exist(\@output_files);
    my $ok = 1;
    #IF DON"T EXIST CREATE OUTPUT FOLDER
    if ($exists == 1) {
        #Maybe you want to add additonal checks on log files before you return the output files
        if (!defined $overwrite) {
            &print_info_message("Output already exists: ".join(",", @output_files));
            $ok =0
        }
    } else {
        unless (-e $output_dir) {
            &make_path($output_dir);
        }
    }
    
    if ($#input_files == -1) {
        &print_error("No input files found.\n");
    }
    
    return $ok;
}


#CHECK IF INPUT DIRECTORIES CONTAIN MULITIPLE FILES BUT ARE ONE BULK SAMPLE
#IF THIS IS THE CASE THROW ERRORß
sub check_multiple_files_bulk {
    my @input_files = @{$_[0]};
    my $index = $_[1];
    if ($#input_files > 0 && $index eq $BULK_NAME_LABEL) {  &print_error("Mulitple input files found (only one needed): ".join("\n", @input_files)."\n"); }
}

#CHECK IF FILES HAVE BEEN GENERATED AND EXIST
sub post_checks {
    my @output_files = @{$_[0]};
    #CHECK IF OUTPUT EXISTS
    #MAYBE YOU WANT TO DO ADDITIONAL CHECKS HERE OF THE LOG FILES
    if (!&files_exist(\@output_files) == 1) {
        &print_error("No output created..\n");
    } else {
        &print_info_message("Output: ".join(", ", @output_files));
    }
    
}

#GET EXTENSION OF FILE NAME
#ALL EXTENSIONS .ext1.ext2.ext3
sub get_file_ext {
    my $input_f = $_[0];
    my ($ext) = $input_f =~ /((\.[^.\s]+)+)$/;
    return $ext;
}

#GET FILE NAME WITH EXTENSIONS
sub get_file_name {
    my $input_f = $_[0];
    my @dirs = split('/', $input_f);
    my $file = $dirs[$#dirs];
    return $file;
}

#GET FILE NAME WITHOUT EXTENSIONS
sub get_file_name_no_ext {
    my $file = &get_file_name($_[0]);
    my @split = split('\.', $file);
    my $name = $split[0];
    return $name;
}

#CHECK IF ARRAY IS EMPTY
sub is_empty_array {
    my @files = @{$_[0]};
    if ($#files > -1) {
        return 0;
    } else {
        return 1;
    }
}

#PRINT INFORMATION MESSAGE
sub print_info_message {
    print("[INFO:] $_[0]\n");
}


#WRITE STRING TO FILE
sub write_string_to_file {
    
    my $string = $_[0];
    my $out_file = $_[1];
    open (MYFILE, $_[2], $out_file) or die "Error opening $out_file: $!";
    print MYFILE $string;
    close(MYFILE);
    
}

#EXECUTE COMMAND VERSION 1
sub exec_cmd {
    
    my @cmd = @{$_[1]};
    my $log =  $_[0];
    
    my $cmd_string = join(' ', @cmd);
    my $output = "[Executed command]: $cmd_string\n";
    &write_string_to_file($output, $log, ">>");
    
    $output = `@cmd 2>&1`;
    &write_string_to_file($output, $log, ">>");
}

#EXECUTE COMMAND VERSION 2
sub exec_cmd_redirect {
    my @cmd = @{$_[1]};
    my $log =  $_[0];
    my $output = $_[2];
    my $com = "[Executed command]:".join(" ", @cmd)."\n";
    &write_string_to_file($com, $log, ">>");
    
    push(@cmd, "1> $output");
    push(@cmd, "2>> $log");
    
    my $code = system("@cmd");
    &write_string_to_file("[Exit code]: $code\n", $log, ">>");
    return ($code);
}



#RETURN 0 iF FILE DOESNT EXIST
sub files_exist {
    my @files = @{$_[0]};
    my $exists = 1;
    foreach my $item (@files) {
        unless (-e $item) {
            $exists = 0;
            last;
        }
    }
    return ($exists);
}

#GO INTO DIR, LIST ALL FILES, AND PICK ALL THAT END WITH .$extension.
sub grep_input_file {
    my $directory = $_[0];
    my $index = $_[1];
    my @extensions = @{$_[2]};
   
    opendir (DIR, $directory) or die $!;
    my @files;
    
    while (my $file = readdir(DIR)) {
        
        for my $ext (@extensions){
            if(&check_naming_convention($ext, $file, $index)) {
                push(@files, $directory."/".$file);
                last;
            }
        }
    }
    return @files;
}

#CHECK IF NAMING CONVENTION CORRECT OF FILES
sub check_naming_convention {
    my $ext = $_[0];
    my $file = $_[1];
    my $index = $_[2];
    my $found = 0;
    
    #NAMING CONVENTION FOR FASTQ AND FASTQ.GZ FILES
    #BASICALLY THEY CONTAIN _1 or _2 before extension
    if ($ext eq $EXT_FASTQ | $ext eq $EXT_FASTQ_GZ) {
        if ($file =~ /(.*#)(\d+)(.*)/g) {
            $found = $index eq $2;
        }        
    #ALL OTHERS
    }else {
        if($file =~ /(.*[^1-9]$index|^$index)$ext$/) {
            $found = 1;
        }
    }
    return ($found);
}

#CONCATINATE FILES TO OUTPUTFILE
sub concat_files {
    my @files = @{$_[0]};
    my $output_file = $_[1];
    
    if ($#files != -1) {
        
        open (MYFILE, '>', $output_file) or die "[ERROR]: Can't write to file '$output_file' $!";
        foreach (@files) {
            open(FILE, $_) || ((warn "Can't open file $_\n"));
            while (<FILE>) {
                print MYFILE;
            }
            close(FILE);
        }
        close(MYFILE);
    }
}

#CHECK IF FOLDER IS EMPTY
sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}


sub assign_variables {
    my %Config = %{$_[0]};

#while (my ($k,$v)=each %Config){print "$k $v\n"}
#Aexit 2;
$ROOT_DIR = $Config{"DIRECTORIES.ROOT_DIR"};
$TEMP_DIR=$Config{"DIRECTORIES.TEMP_DIR"};;
$REF_DIR = $Config{"DIRECTORIES.REF_DIR"};
$TOOLS_MAPPING = $Config{"DIRECTORIES.TOOLS_MAPPING"};
$TOOLS_QUANTIFICATION = $Config{"DIRECTORIES.TOOLS_QUANTIFICATION"};

#EXTENSIONS
$EXT_BAM = $Config{"EXTENSIONS.EXT_BAM"};
$EXT_FASTQ= $Config{"EXTENSIONS.EXT_FASTQ"};
$EXT_FASTQ_GZ = $Config{"EXTENSIONS.EXT_FASTQ_GZ"};
$EXT_SAM = $Config{"EXTENSIONS.EXT_SAM"};
$EXT_SORTED = $Config{"EXTENSIONS.EXT_SORTED"};
$EXT_COUNTS = $Config{"EXTENSIONS.EXT_COUNTS"};
$EXT_METRICS = $Config{"EXTENSIONS.EXT_METRICS"};
 $EXT_MERGED = $Config{"EXTENSIONS.EXT_MERGED"};
 $EXT_MARKED_DUPL= $Config{"EXTENSIONS.EXT_MARKED_DUPL"};
 $EXT_FASTA = $Config{"EXTENSIONS.EXT_FASTA"};
 $EXT_LOG = $Config{"EXTENSIONS.EXT_LOG"};
 $EXT_GTF= $Config{"EXTENSIONS.EXT_GTF"};
 $EXT_GSNAP_SPLICE = $Config{"EXTENSIONS.EXT_GSNAP_SPLICE"};
 $EXT_GSNAP_IIT = $Config{"EXTENSIONS.EXT_GSNAP_IIT"};
 $EXT_DICT = $Config{"EXTENSIONS.EXT_DICT"};
 $EXT_SUMMARY = $Config{"EXTENSIONS.EXT_SUMMARY"};
 $FORWARD_STRAND = $Config{"EXTENSIONS.FORWARD_STRAND"};
 $REVERSE_STRAND = $Config{"EXTENSIONS.REVERSE_STRAND"};


#BINARY SOFTWARES
#MAPPING
 $GSNAP=$Config{"SOFTWARE.GSNAP"};
 $GSNAP_BUILD=$Config{"SOFTWARE.GSNAP_BUILD"};
 $GTF_SPLICE_TOOL = $Config{"SOFTWARE.GTF_SPLICE_TOOL"};
 $GTF_IIT_TOOL = $Config{"SOFTWARE.GTF_IIT_TOOL"};
 $HTSEQTOOL=$Config{"SOFTWARE.HTSEQTOOL"};
 $BOWTIE1_BUILD=$Config{"SOFTWARE.BOWTIE1_BUILD"};
 $BOWTIE1=$Config{"SOFTWARE.BOWTIE1"};
 $BOWTIE2_BUILD=$Config{"SOFTWARE.BOWTIE2_BUILD"};
 $BOWTIE2=$Config{"SOFTWARE.BOWTIE2"}; 
 $STAR=$Config{"SOFTWARE.STAR"};
 $BWA_BUILD=$Config{"SOFTWARE.BWA_BUILD"};
 $STAR_BUILD=$Config{"SOFTWARE.STAR_BUILD"};
 $STAR=$Config{"SOFTWARE.STAR"};
 $SALMON_BUILD=$Config{"SOFTWARE.SALMON"};
 $BWA=$Config{"SOFTWARE.BWA"};
#QUANTIFICATION
 $TOPHAT=$Config{"SOFTWARE.TOPHAT"};
 $CUFFLINKS=$Config{"SOFTWARE.CUFFLINKS"};
 
#HELPER TOOLS
$PICARD_TOOL = $Config{"SOFTWARE.PICARD_TOOL"};
$BAM_2_FASTQ=$Config{"SOFTWARE.BAM_2_FASTQ"};
$SAMTOOLS=$Config{"SOFTWARE.SAMTOOLS"};

}
#CHECK ALL PARAMETERS
#THERE IS SOME WORKAROUND FOR TOPHAT HERE
sub check_parameters {
    my %parameters = %{$_[0]};
    
    #LOAD CONFIG FILE   
    $CONFIG_FILE = $parameters{$CONFIG_FILE_LABEL}; 
    if (!defined $CONFIG_FILE || $CONFIG_FILE_LABEL eq "") {
        &print_error("Config file not specified");
    }
    my %Config;
    Config::Simple->import_from($CONFIG_FILE, \%Config);
    &assign_variables(\%Config);
    
    #AVAILABLE SOFTWARE AND REFERENCE GENOMES
    %available_genomes = &list_dirs($REF_DIR);
    %available_counter = &list_dirs($TOOLS_QUANTIFICATION);
    %available_mapper = &list_dirs($TOOLS_MAPPING);
    
    #TWO PROGRAMS. EITHER CREATING REFERENCE GENOME DATABASES OR WORKING ON RNA_SEQ DATA
    #DATABASE CREATION
    
    #IF REF DB CREATION
    if (exists $parameters{$BUILD_DB_LABEL}) {
        
        #CHECK IF GTF AND GENOME FASTA FILE HAVE BEEN SET
        my $genome_file =  $parameters{$BUILD_DB_LABEL};
        unless (-e $genome_file) {
            if (!defined $genome_file) { $genome_file="";}
            &print_error("Reference genome file does not exist '$genome_file'");
        }
        my $gtf;
        $gtf = $parameters{$GTF_FILE_LABEL};
        if (!(defined $gtf)) {
            &print_error("Gtf file (-gtf) not specified.");
        }
        unless (-e $gtf) {
            &print_error("Gtf file '$gtf' does not exist.");
        }
        
        #CHECK IF MAPPING TOOL WAS SPECIFIED
        my $mapper = $parameters{$MAPPING_TOOL_LABEL};
        &check_object_specified("Mapping tool (-m)", $mapper, \%available_mapper);
        
    #MAPPING, QUANTIFICATION
    } else {
        #check input defined
        my $input = $parameters{$USER_IN_DIR_LABEL};
        $input =~ s/\///g;
        $input = "" unless defined $input;
        &dir_path_ok($input);
        &dir_exists($ROOT_DIR.$input);

        if (!(defined $input)) {
            &print_error("Input file (-i) not specified.");
        }
        
        #Check if running bulk or single cell data
        $parameters{$JOB_INDEX_LABEL} = &check_job_id($parameters{$JOB_INDEX_LABEL});
        #check output  defined
        my $output=$parameters{$USER_OUT_DIR_LABEL};
        if (!(defined $output)) {
            $output = "/output";
        }
        
        $output = $input."/".$output;
        $parameters{$USER_OUT_DIR_LABEL} = $output;
        
        my $mapper = "";
        
        #IF MAPPING WAS SELECTED
        if (exists $parameters{$MAPPING_TOOL_LABEL}) {
            
            #Check if raw folder existis
            my $raw_folder = $ROOT_DIR.$input."/raw";
            unless (-d $raw_folder) {
                &print_error("Raw folder does not exist. Please create and put raw data inside.: $raw_folder");
            }
            
            #THE PROBLEM WITH TOPHAT IS THAT IT RELIES ON BOWTIE
            #HENCE IT CAN'T BE TREATED INDIVIDUALLY.
            #THEREFORE IF TOPHAT WAS SELECTED, IT LOOKS UP WHICH
            #BOWTIE VERSION IS THE NEWEST AND GETS ALL AVAILABLE
            #REFERENCE GENOME FOR BOWTIE
            $mapper =  $parameters{$MAPPING_TOOL_LABEL};
            if (&string_contains($mapper, "tophat")) {
                $mapper = &get_latest_mapper(\%available_mapper, "bowtie2", "tophat", 0);
            }
            &check_object_specified("Mapping tool (-m)", $mapper, \%available_mapper);
            #FILTER OUT GENOMES THAT ARE NOT AVAIALBLE FOR A SPECIFIC MAPPING TOOL
            my %mapper_available_genomes;
            while ((my $key, my $value) = each %available_genomes) {
                my %g = &list_dirs($REF_DIR."/".$key);
                if (exists $g{$mapper}) {
                    $mapper_available_genomes{$key}=1;
                }
            }
            
            #THROW ERROR IF GENOME NOT AVAILABLE
            my $used_genome =  $parameters{$GENOME_LABEL};
            &check_object_specified("Reference genome (-g)", $used_genome, \%mapper_available_genomes);
            if(exists $parameters{$GENOME_LABEL} && exists $parameters{$BUILD_DB_LABEL}) {
                &print_error("Can't use -g and -build_db option at the same time.");
            }
        }
        
        
        #CHECK IF QUANTIFICATION TOOL HAS BEEN SELECTED
        my $counter;
        if (exists $parameters{$COUNT_TOOL_LABEL}) {
            $counter =  $parameters{$COUNT_TOOL_LABEL};
            &check_object_specified("Quantification tool (-h)", $counter, \%available_counter);
            my $used_genome =  $parameters{$GENOME_LABEL};
            &check_object_specified("Reference genome (-g)", $used_genome, \%available_genomes);
        }
    }
    return (%parameters);
}

#CHECK IF OBJECT HAS BEEN SPECIFIED
sub check_object_specified {
    my $object = $_[0];
    #OBJECTS AVAILABLE
    my %objects_available = %{$_[2]};
    
    #SET TO EMPTY IF NOTHING HAS BEEN SPECIFIED
    my $object_specified = $_[1];
    if (!defined $object_specified) {
        $object_specified ="";
    }
    
    if (keys %objects_available == 0) {
        &print_error("No $object are available. Please create one.");
    } else {
        #IF OBJECT IS NOT AVAILABLE
        if (!exists($objects_available{$object_specified})) {
            my $counter = 0;
            my @objects_available_num;
            #NUMERATE ALL AVAILABLE OBJECTS
            foreach my $key (sort {lc $a cmp lc $b} keys %objects_available) {
                push(@objects_available_num, ++$counter.". ".$key)
            }
            #AND SUGGEST TO USER
            &print_error("$object not available or not specified. Create $object or choose from:\n".join("\n", @objects_available_num));
        }
    }
}

#GET LATEST MAPPER VERSION
sub get_latest_mapper {
    my %available_mapper = %{$_[0]};
    my $prefix = $_[1];
    my $software = $_[2];
    my $quite = $_[3];
    my @available = keys %available_mapper ;
    @available = sort @available;
    for (my $i = 0; $i <= $#available; $i++) {
        if(&string_contains($available[$i], $prefix)) {
            if ($quite == 1) {
                &print_info_message("$software requires $prefix. Chose latest installed version $available[$i]");
            }
            return ($available[$i]);
        }
    }
    &print_error("$prefix not available to run $software. Please install a version of it");
}

#LIST ALL DIRECTORIES IN SELECTED DIRECTORY. NOT RECURSIVE
sub list_dirs {
    my $dirname = $_[0];
    opendir my($dh), $dirname or die "Couldn't open dir '$dirname': $!";
    
    my @dirs = grep {-d "$dirname/$_" && ! /^\.{1,2}$/} readdir($dh);
    closedir $dh;
    
    my %dirs = map { $_ => 1 } @dirs;
    return %dirs;
    
}


#CHECK AND SET PROPER JOB ID
sub check_job_id {
    
    my $job = $_[0];
    if (!(defined $job && $job ne "")) {
        $job = $ENV{'LSB_JOBINDEX'};
        if (!(defined $job && $job ne "")) {
	    $job = $ENV{'SGE_TASK_ID'};
        }
	if (!(defined $job && $job ne "")) {
		$job = $BULK_NAME_LABEL;
	}
    }
    return ($job);
}
sub write_log_file {
    my %parameters = %{$_[0]};
    
    foreach my $param (sort keys %parameters) {
        print "$param: $parameters{$param}\n";
    }
}


sub write_to_file {
    
    my $file = $_[0];
    my $output = $_[1];
    # check if the file exists
    if (-f $file) {
        unlink $file
        or croak "Cannot delete $file: $!";
    }
    
    
    # use a variable for the file handle
    my $OUTFILE;
    
    # use the three arguments version of open
    # and check for errors
    open $OUTFILE, '>>', $file
    or croak "Cannot open $file: $OS_ERROR";
    
    
    # you can check for errors (e.g., if after opening the disk gets full)
    print { $OUTFILE } "Something\n"
    or croak "Cannot write to $file: $OS_ERROR";
    
    # check for errors
    close $OUTFILE
    or croak "Cannot close $file: $OS_ERROR";
}

sub print_error {
    
    my $message = $_[0];
    print "[ERROR]: $message\n";
    exit 1;
}


sub dir_path_ok{
    
    my $dir = $_[0];
    if ($dir !~ /\w+/) {
        &print_error("Invalid directory path: '$dir'");
    }
}



sub dir_exists {
    my $dir = $_[0];
    unless (-d $dir) {
        &print_error("Directory does not exist: '$dir'");
    }
}

sub string_contains {
    my $origin = lc($_[0]);
    my $sub_string = lc($_[1]);
    if (defined $origin && defined $sub_string && index($origin, $sub_string) != -1) {
        return 1;
    }
    return 0;
    
}

#PARSE COMMAND LINE INPUT
sub parse_command_line {
    
    my @arguments = @{$_[0]};
    my %parameters = ();
    
    if ($#arguments == -1) {
        &print_program_info();
    }
    for (my $i = 0; $i<= $#arguments; $i++) {
        my $arg = $arguments[$i];
        switch ($arg) {
            case "-i"
            {
                $parameters{$USER_IN_DIR_LABEL} =  $arguments[++$i];
            }
            case "-o"
            {
                $parameters{$USER_OUT_DIR_LABEL} = $arguments[++$i];
                
            }
            case "-r"
            {
                $parameters{$OVERWRITE_LABEL} =  "yes";
                
            }
            case "-index"
            {
                $parameters{$JOB_INDEX_LABEL} =  $arguments[++$i];
                my $index = $parameters{$JOB_INDEX_LABEL};
                $index = "" unless defined $index;
		$parameters{$JOB_INDEX_LABEL} = $index;

                if ($index !~ /[0-9]+/) {
                    &print_error("Invalid job index: '$index'");
                }
            }
            case "-m"
            {
                $parameters{$MAPPING_TOOL_LABEL} =  $arguments[++$i];
                
            }
            
            case "-g"
            {
                $parameters{$GENOME_LABEL} =  $arguments[++$i];
            }
            
            case "-gtf"
            {
                $parameters{$GTF_FILE_LABEL} =  $arguments[++$i];
            }
            case "-build_db"
            {
                $parameters{$BUILD_DB_LABEL} =  $arguments[++$i];
            }
            case "-q"
            {
                $parameters{$COUNT_TOOL_LABEL} =  $arguments[++$i];
            }
            
            case "-margs"
            {
                $parameters{$MAPPING_ARGS_LABEL} =  $arguments[++$i];
            }
            case "-c"
            {
                
                $parameters{$CONFIG_FILE_LABEL} =  $arguments[++$i];
            }
            case "-qargs"
            {
                $parameters{$COUNT_ARGS_LABEL} = $arguments[++$i];
            }
            else { &print_error("Invalid option: '$arg'");
            }
        }
    }
    return %parameters;
}
sub print_program_info
{
    print "SINGLE CELL RNA SEQUENCING PIPELINE 1.0\n";
    print "Usage: perl pipeline.pl [options]\n\n";
    print "This pipeline is designed for RNA-seq data. Paired and single end, bulk and single cell data. For questions contact ti1\@sanger.ac.uk or look at the code\n";
    print "Options:\n";
    print "RNA-seq data\n";
    print "\t-i INPUT_DIRECTORY\tInput directory from $ROOT_DIR\n";
    print "\t-o OUTPUT_DIRECTORY\tOutput directory in $TEMP_DIR\n";
    print "\t-r\tOverwrite results\n";
    print "\t-index 0-9*\tIndex to use. Generally for single cell data only. Each index will be taken as file name for cells. Only specify if you want to run a particular cell.\n";
    print "\t-m MAPPER\tMapping program to use for mapping (Run with empty selection to see available programs). Requires -i, -g option\n";
    print "\t-margs MAPPER_ARGUMENTS\tAdditional arguments that will be directly passed to the executable. Please only used if you know what you are doing.\n";
    print "\t-g GENOME\tGenome to use (Only in combination with -m or -c)\n";
    print "\t-c COUNT_TOOL\tTool to use for gene count table creation (Run with empty selection to see available programs). Requires -i,-g option\n";
    print "\t-cargs COUNT_ARGUMENTS\tAdditional arguments that will be directly passed to the executable. Please only used if you know what you are doing.\n";
    print "\n\n";
    
    print "Reference genome creation\n";
    print "\t-build_db MAPPER\tMapping program to use to create index files (Run with empty selection to see available programs). Requires -g, -gtf\n";
    print "\t-g REFERENCE_FILE\tFile used to create reference genome. Cannot be empty.\n";
    print "\t-gtf GTF_FILE\t File used to create reference genome. Cannot be empty.\n";
    exit;
}


&run();
