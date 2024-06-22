/**
 * @file    script/Cscript/scanning_signalsnp.c
 * @brief   script description
 * 
 * @author  hanfc
 * @date    2024/06/21
 * @license MIT License
 * 
 * @version 1.0
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>


#define MAX_LINE_LEN 1024


#define INFO    "INFO"
#define ERROR   "ERROR"
#define WARNING "WARNING"

void log_print(const char *level, const char *fmt, ...) {
    va_list args;
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    char time_buffer[80];
    strftime(time_buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);

    fprintf(stderr, "[%s] %s: ", time_buffer, level);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
}


// Function prototypes
void print_usage(const char *prog_name) {
    fprintf(stdout, "Usage:%s -g <gemma_output_file> -s <snpEff_annotation_file> -n <sample_number> [-t <threshold>] [-pre <prefix>] [-o <output_path>] [-h]\n", prog_name);
    fprintf(stdout, "Required options:\n");
    fprintf(stdout, "   -g, --gemma  Intput the result of gemma analysis\n");
    fprintf(stdout, "   -s, --snpAnn    Intput the result of snpEff annotation\n");
    fprintf(stdout, "   -n, --number  The sample number in gemma model\n");

    fprintf(stdout, "Optional options:\n");

    fprintf(stdout, "   -t, --threshold  The threshold for p-value. default: 0.05/total_snps\n");
    fprintf(stdout, "   -pre, --prefix  Prefix of the output\n");
    fprintf(stdout, "   -o, --output The output path\n");
    fprintf(stdout, "   -h, --help      Display this help message\n");
}



void parse_arguments(int argc, char *argv[], char *gemma_file, char *prefix, char *snpAnn_file, char *output_file, float *threshold, int *number) {
    int i;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--gemma") == 0) {
            if (i + 1 < argc) {
                // genbank_file = argv[++i];
                char *gemma_path = argv[++i];
                strcpy(gemma_file, gemma_path);
            } else {
                log_print(ERROR, "Missing gemma output file argument");
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-pre") == 0 || strcmp(argv[i], "--prefix") == 0) {
            if (i + 1 < argc) {
                // *prefix = argv[++i];
                char *prefix_path = argv[++i];
                strcpy(prefix, prefix_path);
            }
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--snpAnn") == 0) {
            if (i + 1 < argc) {
                char *snpAnn_path = argv[++i];
                strcpy(snpAnn_file, snpAnn_path);
            }
        } else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--theshold") == 0) {
            if (i + 1 < argc) {
                *threshold = atof(argv[++i]);
            }
        } else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--number") == 0) {
            if (i + 1 < argc) {
                *number = atoi(argv[++i]);
            }                    
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                char *output_path = argv[++i];
                strcpy(output_file, output_path);
                // output_file = argv[++i];
            }
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(EXIT_SUCCESS);
        } else {
            log_print(ERROR, "Invalid option '%s'", argv[i]);
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (strlen(gemma_file) == 0 || strlen(snpAnn_file) == 0 || *number == 0) {
        log_print(ERROR, "Please provide all required arguments");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    } else {
        char *ext = strrchr(gemma_file, '.'); // get the extension of the file name
        if (strcmp(ext, ".txt") != 0) { // check if the extension is "assoc.txt"
            log_print(ERROR, "gemma output file must have a extension (assoc.txt)");
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        } else {
            size_t bnmlen = ext - gemma_file; // get the length of the base name
            char *bnm = (char *)malloc(bnmlen +  1);

            if (bnm == NULL) 
            {
                log_print(ERROR, "Failed to allocate memory for base name");
                exit(EXIT_FAILURE);
            }
            strncpy(bnm, gemma_file, bnmlen);
            bnm[bnmlen] = '\0';

            if (strlen(prefix) == 0) {
                // prefix = strrchr(bnm, '/') + 1;
                if (strstr(bnm, "/") == NULL) {
                    strcpy(prefix, bnm);
                } else {
                    strcpy(prefix, strrchr(bnm, '/') + 1);
                }
            }

            if (strlen(output_file) == 0)
            {
                if (strstr(bnm, "/") == NULL) {
                    strcpy(output_file, "./");
                } else {
                    char *last_prefix = strrchr(bnm, '/') + 1;
                    size_t outlen = last_prefix - bnm;
                    strncpy(output_file, bnm, outlen);
                    output_file[outlen] = '\0';
                }
            } 
            free(bnm);
        }
    }
}


typedef struct {
    char *chr;
    char *snp;
    float p_wald;
    float pve;
    float log_p_value;
} snp_info;


void read_gemma_file(char *gemma_file, snp_info **snp_list, int *snp_num, float threshold, int number) {

    FILE *fp = fopen(gemma_file, "r");
    if (fp == NULL) {
        log_print(ERROR, "Failed to open file %s", gemma_file);
        exit(EXIT_FAILURE);
    }


    char line[MAX_LINE_LEN];
    int line_num = 0;
    float Pthreshold = 0.0;

    while (fgets(line, MAX_LINE_LEN, fp)!= NULL) {
        line_num++;
    }
    log_print(INFO, "Total snps: %d", line_num);

    if (threshold == 0.0) {
        threshold = -log10(0.05 / (float)line_num);
        log_print(INFO, "The threshold (0.05 / total snps) for p-value (-log10) is: %f", threshold);
        Pthreshold = pow(10, -threshold);
        
    } else {
        log_print(INFO, "The threshold for p-value (-log10) is: %f", threshold);
        Pthreshold = pow(10, -threshold);
    }

    // rewind the file
    rewind(fp);

    *snp_list = malloc(100 * sizeof(snp_info));
    int snp_list_size = 100;

    char chr[100];
    char rs[100];
    float beta = 0.0;
    float af = 0.0;
    float se = 0.0;
    int n_miss = 0.0;
    float p_wald = 0.0;

    float pve = 0.0;
    while (fgets(line, MAX_LINE_LEN, fp)!= NULL) {
        sscanf(
                line, "%s %s %*d %d %*s %*s %f %e %e %*e %*e %e", 
                chr,
                rs, 
                &n_miss, 
                &af, 
                &beta,
                &se,
                &p_wald
            );
        // skip the header line
        if (strstr(line, "p_wald")) {
            continue;
        } else if (p_wald <= Pthreshold) {
            // calculate the Pve
            pve = (2 * (pow(beta, 2) * af * (1 - af))) / (2 * pow(beta, 2) * af * (1 - af) + pow(se, 2) * 2 * (number - n_miss) * af * (1 - af));
            // log_print(INFO, "SNP: %s, P-value: %e, Pve: %f", rs, p_wald, pve);
            if (*snp_num >= snp_list_size) {
                snp_list_size += 100;
                *snp_list = realloc(*snp_list, snp_list_size * sizeof(snp_info));
                if (*snp_list == NULL) {
                    log_print(ERROR, "Failed to reallocate memory for snp_list");
                    exit(EXIT_FAILURE);
                }
            }

            (*snp_list)[*snp_num].snp = malloc(strlen(rs) + 1);
            (*snp_list)[*snp_num].chr = malloc(strlen(chr) + 1);
            if ((*snp_list)[*snp_num].snp == NULL || (*snp_list)[*snp_num].chr == NULL) {
                log_print(ERROR, "Failed to allocate memory for snp");
                exit(EXIT_FAILURE);
            }
            
            strcpy((*snp_list)[*snp_num].snp, rs);
            strcpy((*snp_list)[*snp_num].chr, chr);
            (*snp_list)[*snp_num].p_wald = p_wald;
            (*snp_list)[*snp_num].pve = pve;
            (*snp_list)[*snp_num].log_p_value = -log10(p_wald);
            (*snp_num)++;
        }
    }
    fclose(fp);
}


void read_snpAnn_file(char *snpAnn_file, snp_info *snp_list, int *snp_num, char *output_file) {
    FILE *fw = fopen(output_file, "w");
    if (fw == NULL) {
        log_print(ERROR, "Failed to open file %s", output_file);
        exit(EXIT_FAILURE);
    }
    FILE *fr = fopen(snpAnn_file, "r");
    if (fr == NULL) {
        log_print(ERROR, "Failed to open file %s", snpAnn_file);
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LEN];
    int line_num = 0;

    int pos;
    char id[100];
    char ref[10];
    char alt[10];
    char ann[1024];

    fprintf(fw, "Chr\tPos\tID\tRef\tAlt\tPve\tP_wald\t-log(10)\tAnn\n");
    while (fgets(line, MAX_LINE_LEN, fr)!= NULL) {
        if (strstr(line, "##")) {
            continue;
        } else {
            sscanf(line, "%*s %d %s %s %s %*s %*s %s", &pos, id, ref, alt, ann);
            for (int i = 0; i < *snp_num; i++) {
                if (strcmp(snp_list[i].snp, id) == 0) {
                    // log_print(INFO, "SNP: %s, Annotation: %s", id, ann);
                    // write to file
                    fprintf(fw, "%s\t%d\t%s\t%s\t%s\t%f\t%e\t%f\t%s\n", snp_list[i].chr, pos, id, ref, alt, snp_list[i].pve, snp_list[i].p_wald, snp_list[i].log_p_value,ann);
                }
            }
        }
    }
    fclose(fr);
    fclose(fw);
}

    


int main(int argc, char *argv[]) {
    char *gemma_file = malloc(1024);
    char *snpAnn_file = malloc(1024);
    char *prefix = malloc(1024);
    char *output_file = malloc(1024);
    float threshold = 0.0;
    int number = 0;
    parse_arguments(argc, argv, gemma_file, prefix, snpAnn_file, output_file, &threshold, &number);

    if (access(gemma_file, F_OK) == -1) {
        log_print(ERROR, "%s does not exist!", gemma_file);
        print_usage(argv[0]);
        free(gemma_file);
        free(prefix);
        free(output_file);
        free(snpAnn_file);
        exit(1);
    }    

    if (access(snpAnn_file, F_OK) == -1) {
        log_print(ERROR, "%s does not exist!", snpAnn_file);
        print_usage(argv[0]);
        free(gemma_file);
        free(prefix);
        free(output_file);
        free(snpAnn_file);
        exit(1);
    }

    if (access(output_file, F_OK) == -1) {
        log_print(ERROR, "Output Path does not exist!");
        print_usage(argv[0]);
        free(gemma_file);
        free(prefix);
        free(output_file);
        free(snpAnn_file);
        exit(1);
    }

    log_print(INFO, "Gemma file: %s", gemma_file);
    log_print(INFO, "SnpEff annotation file: %s", snpAnn_file);
    log_print(INFO, "Prefix: %s", prefix);
    log_print(INFO, "Output file: %s", output_file);
    log_print(INFO, "Sample number: %d", number);

    // TODO: implement the scanning of signalsnp
    snp_info *snp_list = NULL;
    int snp_num = 0;
    read_gemma_file(gemma_file, &snp_list, &snp_num, threshold, number);


    if (snp_num > 0) {
        // printf("______ ______ ______\n");
        log_print(INFO, "The number of signal snps found is %d", snp_num);
        printf(" -------------  ------------  ------------  ---------\n");
        printf(" SNP            P-value       Pve           -Log\n");
        printf(" -------------  ------------  ------------  ---------\n");
        for (int i = 0; i < snp_num; i++) {

            printf(" %s | %e | %f | %f\n", snp_list[i].snp, snp_list[i].p_wald, snp_list[i].pve, snp_list[i].log_p_value);
            // log_print(INFO, "SNP: %s, P-value: %e, Pve: %f", snp_list[i].snp, snp_list[i].p_wald, snp_list[i].pve);
        }
        printf(" -------------  ------------  ------------  ---------\n");
    } else {
        log_print(ERROR, "No signal snps found");
        exit(EXIT_FAILURE);
    }

    // 
    size_t len = strlen(output_file);
    if (len > 0 && output_file[len - 1] != '/') {
        strcat(output_file, "/");
    }
    strcat(output_file, prefix);
    strcat(output_file, ".scanning_signalsnp.txt");
    log_print(INFO, "getting the annotation of signal snps from SnpEff annotation file ...");
    read_snpAnn_file(snpAnn_file, snp_list, &snp_num, output_file);
    log_print(INFO, "The result file is %s", output_file);
    log_print(INFO, "Done");

    // free memory
    if (snp_list ) free(snp_list); 
    free(gemma_file);
    free(snpAnn_file);
    free(prefix);
    free(output_file);

    return 0;

}