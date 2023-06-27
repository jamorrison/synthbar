/*
 * The MIT License
 *
 * Copyright (c) 2022-2023 Jacob Morrison <jacob.morrison@vai.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <sys/time.h>
#include <zlib.h>

#include "kstring.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define SB_VERSION "1.0.0" /* synthbar version */
#define N_EXTRA_CHARS 7    /* 4 newlines + 1 space + 1 separator + 1 null-terminator */

// What the function name says!
// Returns time in seconds
static inline double get_current_time() {
    struct timeval  tp;
    struct timezone tzp;

    gettimeofday(&tp, &tzp);

    return tp.tv_sec + 1e-6*tp.tv_usec;
}

// Configuration variables
typedef struct {
    char     *outfn;         /* name of output file */
    char     *barcode;       /* barcode to add to each read */
    uint8_t   umi_first;     /* print the UMI before the barcode in each read */
    uint8_t   remove_linker; /* remove linker (1) or not (0) */
    int32_t   linker_length; /* number of bases in linker */
    int32_t   umi_length;    /* number of bases in UMI */
} sb_conf_t;

// Initialize config variables
sb_conf_t init_sb_conf() {
    sb_conf_t conf = {0};

    conf.outfn         = (char *)"-";
    conf.barcode       = (char *)"CATATAC";
    conf.umi_first     = 0;
    conf.remove_linker = 0;
    conf.linker_length = 6;
    conf.umi_length    = 8;

    return conf;
}

// Print version of code
static int print_version() {
    fprintf(stderr, "Program: synthbar\n");
    fprintf(stderr, "Version: %s\n", SB_VERSION);
    fprintf(stderr, "Contact: Jacob Morrison <jacob.morrison@vai.org>\n");

    return 0;
}

// Print usage information for help
static int usage(sb_conf_t *conf) {
    fprintf(stderr, "\n");
    print_version();
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: synthbar [options] <FASTQ with UMIs>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o, --output STR           name of output file [stdout]\n");
    fprintf(stderr, "Processing Options:\n");
    fprintf(stderr, "    -b, --barcode STR          barcode to prepend to each read [%s]\n", conf->barcode);
    fprintf(stderr, "    -U, --umi-first            add barcode to read after the UMI [off]\n");
    fprintf(stderr, "    -r, --remove-linker        remove linker from read [not removed]\n");
    fprintf(stderr, "    -l, --linker-length INT    length of linker to remove [%i]\n", conf->linker_length);
    fprintf(stderr, "    -u, --umi-length INT       length of UMI before linker [%i]\n", conf->umi_length);
    fprintf(stderr, "    -h, --help                 print usage and exit\n");
    fprintf(stderr, "        --version              print version and exit\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note 1: Input FASTQ can be gzip compressed or uncompressed\n");
    fprintf(stderr, "\n");

    return 0;
}

int main(int argc, char *argv[]) {
    // Init variables
    sb_conf_t conf = init_sb_conf();
    int c;

    // Command line arguments
    static const struct option loptions[] = {
        {"output"       , required_argument, NULL, 'o'},
        {"barcode"      , required_argument, NULL, 'b'},
        {"umi-first"    , no_argument      , NULL, 'U'},
        {"remove-linker", no_argument      , NULL, 'r'},
        {"linker-length", required_argument, NULL, 'l'},
        {"umi-length"   , required_argument, NULL, 'u'},
        {"help"         , no_argument      , NULL, 'h'},
        {"version"      , no_argument      , NULL,  1 },
        {NULL, 0, NULL, 0}
    };

    // Parse CLI
    if (argc < 2) {
        usage(&conf);
        return 0;
    }

    while ((c = getopt_long(argc, argv, "b:l:o:u:Uhrz", loptions, NULL)) >= 0) {
        switch (c) {
            case 'o':
                conf.outfn = optarg;
                break;
            case 'b':
                conf.barcode = optarg;
                break;
            case 'U':
                conf.umi_first = 1;
                break;
            case 'r':
                conf.remove_linker = 1;
                break;
            case 'l':
                conf.linker_length = (int32_t)atoi(optarg);
                break;
            case 'u':
                conf.umi_length = (int32_t)atoi(optarg);
                break;
            case 'h':
                usage(&conf);
                return 0;
            case 1:
                print_version();
                return 0;
            default:
                usage(&conf);
                return 0;
        }
    }

    // Check for input file
    char *infn = optind < argc ? argv[optind++] : NULL;
    if (!infn) {
        usage(&conf);
        fprintf(stderr, "Please provide an input FASTQ\n");
        return 1;
    }

    // Check linker and UMI lengths
    if (conf.umi_length < 0 || conf.linker_length < 0) {
        fprintf(stderr, "Linker (%i) and UMI (%i) lengths must both be >= 0\n", conf.linker_length, conf.umi_length);
        return 1;
    }

    // Init files and handle errors
    gzFile fh1 = gzopen(infn, "r");
    if (!fh1) {
        fprintf(stderr, "Could not open input file: %s\n", infn);
        return 1;
    }

    FILE *oh1 = strcmp(conf.outfn, "-") == 0 ? stdout : fopen(conf.outfn, "w");
    if (strcmp(conf.outfn, "-") != 0 && !oh1) {
        fprintf(stderr, "Could not open output file: %s\n", conf.outfn);
        gzclose(fh1);
        return 1;
    }

    // Create qual string to add
    size_t  bc_len   = strlen(conf.barcode);
    char   *pre_qual = malloc(bc_len + 1);
    memset(pre_qual, 'I', bc_len);
    pre_qual[bc_len] = '\0';

    // Variable initialization
    int        ret_code   = 0;
    kseq_t    *ks1        = kseq_init(fh1);
    int32_t    u_plus_l   = conf.umi_length + conf.linker_length;
    int32_t    link_start = conf.remove_linker ? u_plus_l : conf.umi_length;
    uint32_t   read_count = 0;
    kstring_t *str        = (kstring_t *)calloc(1, sizeof(kstring_t));

    // Process reads
    double t1 = get_current_time();
    while (kseq_read(ks1) >= 0) {
        read_count++;

        // Handle error case of too short read, seq and qual should be same length, so only check seq
        if (conf.remove_linker && ks1->seq.l < u_plus_l) {
            fprintf(stderr, "Read shorter than UMI and linker lengths provided (%li < %i)\n", ks1->seq.l, u_plus_l);
            ret_code = 1;
            goto end;
        }

        // Pre-allocate space, or expand space ahead of time to reduce the number of allocations needed
        size_t str_len = ks1->name.l + ks1->comment.l + ks1->seq.l + ks1->qual.l + (size_t)N_EXTRA_CHARS;

        if (str_len > str->m) {
            int err = ks_resize(str, str_len);
            if (err < 0) {
                ret_code = 1;
                fprintf(stderr, "Unable to reallocate sufficient space\n");
                goto end;
            }
        }

        // Read name
        ksprintf(str, "@%s", ks1->name.s);

        // Read comment (if applicable)
        if (ks1->comment.l > 0) {
            ksprintf(str, " %s", ks1->comment.s);
        }

        // UMI and barcode (seq)
        if (!conf.umi_first) {
            ksprintf(str, "\n%s%.*s", conf.barcode, conf.umi_length, ks1->seq.s);
        } else {
            ksprintf(str, "\n%.*s%s", conf.umi_length, ks1->seq.s, conf.barcode);
        }

        // Linker (seq), sequence, and separator
        ksprintf(str, "%s\n+\n", ks1->seq.s + link_start);

        // UMI and barcode (qual)
        if (!conf.umi_first) {
            ksprintf(str, "%s%.*s", pre_qual, conf.umi_length, ks1->qual.s);
        } else {
            ksprintf(str, "%.*s%s", conf.umi_length, ks1->qual.s, pre_qual);
        }

        // Linker (qual) and quality
        ksprintf(str, "%s\n", ks1->qual.s + link_start);

        // Print out read
        fprintf(oh1, "%s", str->s);

        // Reset kstring for next read
        str->s[0] = '\0';
        str->l = 0;
    }

    goto end;

end:
    double t2 = get_current_time();

    // Clean up
    free(str->s);
    free(str);
    kseq_destroy(ks1);
    free(pre_qual);
    if (strcmp(conf.outfn, "-") != 0) { fclose(oh1); }
    gzclose(fh1);

    fprintf(stderr, "[synthbar:%s] %u reads processed in %.3f seconds (wall time)\n", __func__, read_count, t2-t1);

    return ret_code;
}
