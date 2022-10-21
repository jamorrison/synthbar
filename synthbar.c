/*
 * The MIT License
 *
 * Copyright (c) 2022 Jacob Morrison <jacob.morrison@vai.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <getopt.h>
#include <sys/time.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// What the function name says!
static inline double get_current_time() {
    struct timeval  tp;
    struct timezone tzp;

    gettimeofday(&tp, &tzp);

    return tp.tv_sec + 1e-6*tp.tv_usec;
}

// Configuration variables
typedef struct {
    char     *outfn;         /* name of output file */
    uint8_t   remove_linker; /* remove linker (1) or not (0) */
    uint16_t  linker_length; /* number of bases in linker */
    uint16_t  umi_length;    /* number of bases in UMI */
} sb_conf_t;

// Initialize config variables
sb_conf_t init_sb_conf() {
    sb_conf_t conf = {0};

    conf.outfn = (char *)"-";
    conf.remove_linker = 0;
    conf.linker_length = 6;
    conf.umi_length    = 8;

    return conf;
}

// Print version of code
static int print_version() {
    fprintf(stderr, "Program: synthbar\n");
    fprintf(stderr, "Version: 1.0.0\n");
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
    fprintf(stderr, "    -r, --remove-linker        remove linker from read [not removed]\n");
    fprintf(stderr, "    -l, --linker-length INT    length of linker to remove [%u]\n", conf->linker_length);
    fprintf(stderr, "    -u, --umi-length INT       length of UMI before linker [%u]\n", conf->umi_length);
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

    static const struct option loptions[] = {
        {"output"       , required_argument, NULL, 'o'},
        {"remove-linker", no_argument      , NULL, 'r'},
        {"linker-length", required_argument, NULL, 'l'},
        {"umi-length"   , required_argument, NULL, 'u'},
        {"help"         , no_argument      , NULL, 'h'},
        {"version"      , no_argument      , NULL,  1 },
        {NULL, 0, NULL, 0}
    };

    // Parse CLI
    if (argc < 2) { usage(&conf); return 0; }
    while ((c = getopt_long(argc, argv, "l:o:u:hrz", loptions, NULL)) >= 0) {
        switch (c) {
            case 'o':
                conf.outfn = optarg;
                break;
            case 'r':
                conf.remove_linker = 1;
                break;
            case 'l':
                conf.linker_length = (uint16_t)atoi(optarg);
                break;
            case 'u':
                conf.umi_length = (uint16_t)atoi(optarg);
                break;
            case 'h':
                usage(&conf);
                return 1;
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

    // Read file
    gzFile  fh1 = gzopen(infn, "r");
    FILE   *oh1 = strcmp(conf.outfn, "-")  == 0 ? stdout : fopen(conf.outfn, "w");
    kseq_t *ks1 = kseq_init(fh1);

    uint16_t total_length = conf.umi_length + conf.linker_length;
    double t1 = get_current_time();
    while (kseq_read(ks1) >= 0) {
        fprintf(oh1, "%s %s\nCATATAC", ks1->name.s, ks1->comment.s);

        if (conf.remove_linker) {
            // Handle error case, seq and qual should be same length, so only check seq
            if (ks1->seq.l < conf.umi_length+conf.linker_length) {
                fprintf(stderr, "Read not long enough for UMI and linker lengths provided (%li<%u)\n",
                        ks1->seq.l, conf.umi_length+conf.linker_length);
                kseq_destroy(ks1);
                if (strcmp(conf.outfn, "-") != 0) { fclose(oh1); }
                gzclose(fh1);

                return 1;
            }

            fprintf(oh1, "%.*s", conf.umi_length, ks1->seq.s);
            fprintf(oh1, "%s\n+\nIIIIIII", ks1->seq.s + total_length);

            fprintf(oh1, "%.*s", conf.umi_length, ks1->qual.s);
            fprintf(oh1, "%s\n", ks1->qual.s + total_length);
        } else {
            fprintf(oh1, "%s\n+\nIIIIIII%s\n", ks1->seq.s, ks1->qual.s);
        }
        //fprintf(oh1, "%s %s\nCATATAC%s\n+\nIIIIIII%s\n", ks1->name.s, ks1->comment.s, ks1->seq.s, ks1->qual.s);
    }
    //while (kseq_read(ks1) >= 0) {
    //    ksprintf(&str, "%s %s\nCATATAC", ks1->name.s, ks1->comment.s);

    //    // Create sequence and quality strings without linker if desired
    //    if (conf.remove_linker) {
    //        // Handle error case, seq and qual should be same length, so only check seq
    //        if (ks1->seq.l < conf.umi_length+conf.linker_length) {
    //            fprintf(stderr, "Read not long enough for UMI and linker lengths provided (%li<%u)\n",
    //                    ks1->seq.l, conf.umi_length+conf.linker_length);
    //            kseq_destroy(ks1);
    //            if (strcmp(conf.outfn, "-") != 0) { fclose(oh1); }
    //            gzclose(fh1);

    //            return 1;
    //        }

    //        // Set up seq string
    //        kputsn(ks1->seq.s, conf.umi_length, &str);
    //        kputsn(ks1->seq.s + total_length, ks1->seq.l - total_length, &str);

    //        ksprintf(&str, "\n+\nIIIIIII");

    //        // Set up qual string
    //        kputsn(ks1->qual.s, conf.umi_length, &str);
    //        kputsn(ks1->qual.s + total_length, ks1->qual.l - total_length, &str);
    //    } else {
    //        ksprintf(&str, "%s\n+\nIIIIIII%s", ks1->seq.s, ks1->qual.s);
    //    }

    //    //fprintf(oh1, "%s %s\nCATATAC%s\n+\nIIIIIII%s\n", ks1->name.s, ks1->comment.s, ks1->seq.s, ks1->qual.s);
    //    fprintf(oh1, "%s\n", str.s);

    //    ks_release(&str);
    //}
    double t2 = get_current_time();

    // Clean up
    kseq_destroy(ks1);
    if (strcmp(conf.outfn, "-") != 0) { fclose(oh1); }
    gzclose(fh1);

    fprintf(stderr, "[synthbar:%s] Wall time: %.3f seconds\n", __func__, t2-t1);

    return 0;
}
