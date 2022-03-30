/*
  Copyright (C) 2021  The Blosc Developers <blosc@blosc.org>
  https://blosc.org
  License: BSD 3-Clause (see LICENSE.txt)

  Example program demonstrating use of the Blosc entropy probe codec from C code.

  To compile this program:

  $ gcc entropy_probe.c -o entropy_probe -lblosc2

  To run:

  $ ./entropy_probe

  To run in entropy mode:
  
  $ ./entropy_probe -e output.csv
*/

#include <stdio.h>
#include "blosc2.h"

#define ENTROPY_PROBE_ID 244

#define KB  1024.
#define MB  (1024*KB)
#define GB  (1024*MB)

#define MAX_COPY 32U
#define MAX_DISTANCE 8191
#define MAX_FARDISTANCE (65535 + MAX_DISTANCE - 1)

// The hash length (1 << HASH_LOG2) can be tuned for performance (12 -> 15)
#define HASH_LOG2 (12U)
//#define HASH_LOG2 (13U)
//#define HASH_LOG2 (14U)

#define HASH_FUNCTION(v, s, h) {      \
  v = (s * 2654435761U) >> (32U - h); \
}

#define BLOSCLZ_READU16(p) *((const uint16_t*)(p))
#define BLOSCLZ_READU32(p) *((const uint32_t*)(p))

#define LITERAL2(ip, anchor, copy) {                    \
  oc++; anchor++;                                       \
  ip = anchor;                                          \
  copy++;                                               \
  if (copy == MAX_COPY) {                               \
    copy = 0;                                           \
    oc++;                                               \
  }                                                     \
}


static uint8_t *get_run(uint8_t *ip, const uint8_t *ip_bound, const uint8_t *ref) {
  uint8_t x = ip[-1];
  int64_t value, value2;
  /* Broadcast the value for every byte in a 64-bit register */
  memset(&value, x, 8);
  /* safe because the outer check against ip limit */
  while (ip < (ip_bound - sizeof(int64_t))) {
    value2 = ((int64_t *) ref)[0];
    if (value != value2) {
      /* Return the byte that starts to differ */
      while (*ref++ == x) ip++;
      return ip;
    } else {
      ip += 8;
      ref += 8;
    }
  }
  /* Look into the remainder */
  while ((ip < ip_bound) && (*ref++ == x)) ip++;
  return ip;
}


static uint8_t *get_match(uint8_t *ip, const uint8_t *ip_bound, const uint8_t *ref) {
  while (ip < (ip_bound - sizeof(int64_t))) {
    if (*(int64_t *) ref != *(int64_t *) ip) {
      /* Return the byte that starts to differ */
      while (*ref++ == *ip++) {}
      return ip;
    } else {
      ip += sizeof(int64_t);
      ref += sizeof(int64_t);
    }
  }
  /* Look into the remainder */
  while ((ip < ip_bound) && (*ref++ == *ip++)) {}
  return ip;
}


static uint8_t *get_run_or_match(uint8_t *ip, uint8_t *ip_bound, const uint8_t *ref, bool run) {
  if (run) {
    ip = get_run(ip, ip_bound, ref);
  } else {
    ip = get_match(ip, ip_bound, ref);
  }

  return ip;
}


// Get a guess for the compressed size of a buffer
static float get_cratio(const uint8_t *ibase, int maxlen, int minlen, int ipshift) {
  const uint8_t *ip = ibase;
  int32_t oc = 0;
  const uint16_t hashlen = (1U << (uint8_t) HASH_LOG2);
  uint16_t htab[1U << (uint8_t) HASH_LOG2];
  uint32_t hval;
  uint32_t seq;
  uint8_t copy;
  // Make a tradeoff between testing too much and too little
  uint16_t limit = (maxlen > hashlen) ? hashlen : maxlen;
  const uint8_t *ip_bound = ibase + limit - 1;
  const uint8_t *ip_limit = ibase + limit - 12;

  // Initialize the hash table to distances of 0
  memset(htab, 0, hashlen * sizeof(uint16_t));

  /* we start with literal copy */
  copy = 4;
  oc += 5;

  /* main loop */
  while (ip < ip_limit) {
    const uint8_t *ref;
    unsigned distance;
    const uint8_t *anchor = ip;    /* comparison starting-point */

    /* find potential match */
    seq = BLOSCLZ_READU32(ip);
    HASH_FUNCTION(hval, seq, HASH_LOG2)
    ref = ibase + htab[hval];

    /* calculate distance to the match */
    distance = (unsigned int) (anchor - ref);

    /* update hash table */
    htab[hval] = (uint16_t) (anchor - ibase);

    if (distance == 0 || (distance >= MAX_FARDISTANCE)) {
      LITERAL2(ip, anchor, copy)
      continue;
    }

    /* is this a match? check the first 4 bytes */
    if (BLOSCLZ_READU32(ref) == BLOSCLZ_READU32(ip)) {
      ref += 4;
    } else {
      /* no luck, copy as a literal */
      LITERAL2(ip, anchor, copy)
      continue;
    }

    /* last matched byte */
    ip = anchor + 4;

    /* distance is biased */
    distance--;

    /* get runs or matches; zero distance means a run */
    ip = get_run_or_match((uint8_t*)ip, (uint8_t*)ip_bound, ref, !distance);

    ip -= ipshift;
    int32_t len = (int32_t) (ip - anchor);
    if (len < minlen) {
      LITERAL2(ip, anchor, copy)
      continue;
    }

    /* if we haven't copied anything, adjust the output counter */
    if (!copy)
      oc--;
    /* reset literal counter */
    copy = 0;

    /* encode the match */
    if (distance < MAX_DISTANCE) {
      if (len >= 7) {
        oc += ((len - 7) / 255) + 1;
      }
      oc += 2;
    } else {
      /* far away, but not yet in the another galaxy... */
      if (len >= 7) {
        oc += ((len - 7) / 255) + 1;
      }
      oc += 4;
    }

    /* update the hash at match boundary */
    seq = BLOSCLZ_READU32(ip);
    HASH_FUNCTION(hval, seq, HASH_LOG2)
    htab[hval] = (uint16_t) (ip++ - ibase);
    ip++;
    /* assuming literal copy */
    oc++;
  }

  float ic = (float) (ip - ibase);
  return ic / (float) oc;
}


int entropy_probe(const uint8_t *input, int32_t input_len,
                  uint8_t *output, int32_t output_len,
                  uint8_t meta,
                  blosc2_cparams *cparams, const void *chunk) {
  if (cparams->splitmode != BLOSC_ALWAYS_SPLIT) {
    BLOSC_TRACE_ERROR("Entropy probe can only be used in SPLIT mode. Aborting!");
    return BLOSC2_ERROR_CODEC_PARAM;
  }
  // Get the cratio.  minlen and ipshift are decent defaults, but one can try with (4, 4) or (3, 4) or (4, 3)
  float cratio = get_cratio(input, input_len, 3, 3);
  int cbytes = (int) ((float)input_len / cratio);
  if (cbytes > input_len) {
    cbytes = input_len;
  }
  return cbytes;
}


// This function receives an instrumented chunk having nstreams
int extr_data(FILE *csv_file, int nchunk, uint8_t *chunk, int nstreams, int category) {
  blosc2_instr *instr_data = (blosc2_instr *) chunk;

  for (int nstream = 0; nstream < nstreams; nstream++) {
    float cratio = instr_data->cratio;
    float cspeed = instr_data->cspeed;
    bool special_val = instr_data->flags[0];
    if (!special_val) {
      // Fill csv file
      int special_vals = 0;
      fprintf(csv_file, "%.3g, %.3g, %d, %d, %d\n", cratio, cspeed, special_vals, nchunk, category);
      printf("Chunk %d, block %d: cratio %.3g, speed %.3g\n", nchunk, nstream, cratio, cspeed);
    } else {
        // Fill csv file
        int special_vals = 1;
        fprintf(csv_file, "%.3g, %.3g, %d, %d, %d\n", cratio, cspeed, special_vals, nchunk, category);
        printf("Chunk %d, block %d: cratio %.3g, speed %.3g\n", nchunk, nstream, cratio, cspeed);
    }
    instr_data++;
  }

  return 0;
}

//Categorize the data according to the combination of filter-codec-splitmode
int categorize_data(char codec, char filter, char splitmode) {
    printf(" AAAAAAAAAAAAAAAAAAAAAAAAAAAAA el nombre del codec es: %c\n", codec);

    int category;
    if (codec == BLOSC_BLOSCLZ && filter == BLOSC_NOFILTER && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 0;
    } else if (codec == BLOSC_BLOSCLZ && filter == BLOSC_SHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 1;
    } else if (codec == BLOSC_BLOSCLZ && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 2;
    } else if (codec == BLOSC_LZ4 && filter == BLOSC_NOFILTER && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 3;
    } else if (codec == BLOSC_LZ4 && filter == BLOSC_SHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 4;
    } else if (codec == BLOSC_LZ4 && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 5;
    } else if (codec == BLOSC_LZ4HC && filter == BLOSC_NOFILTER && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 6;
    } else if (codec == BLOSC_LZ4HC && filter == BLOSC_SHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 7;
    } else if (codec == BLOSC_LZ4HC && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 8;
    } else if (codec == BLOSC_ZLIB && filter == BLOSC_NOFILTER && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 9;
    } else if (codec == BLOSC_ZLIB && filter == BLOSC_SHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 10;
    } else if (codec == BLOSC_ZLIB && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 11;
    } else if (codec == BLOSC_ZSTD && filter == BLOSC_NOFILTER && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 12;
    } else if (codec == BLOSC_ZSTD && filter == BLOSC_SHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 13;
    } else if (codec == BLOSC_ZSTD && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_ALWAYS_SPLIT) {
        category = 14;

    } else if (codec == BLOSC_BLOSCLZ && filter == BLOSC_NOFILTER && splitmode == BLOSC_NEVER_SPLIT) {
        category = 15;
    } else if (codec == BLOSC_BLOSCLZ && filter == BLOSC_SHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 16;
    } else if (codec == BLOSC_BLOSCLZ && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 17;
    } else if (codec == BLOSC_LZ4 && filter == BLOSC_NOFILTER && splitmode == BLOSC_NEVER_SPLIT) {
        category = 18;
    } else if (codec == BLOSC_LZ4 && filter == BLOSC_SHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 19;
    } else if (codec == BLOSC_LZ4 && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 20;
    } else if (codec == BLOSC_LZ4HC && filter == BLOSC_NOFILTER && splitmode == BLOSC_NEVER_SPLIT) {
        category = 21;
    } else if (codec == BLOSC_LZ4HC && filter == BLOSC_SHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 22;
    } else if (codec == BLOSC_LZ4HC && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 23;
    } else if (codec == BLOSC_ZLIB && filter == BLOSC_NOFILTER && splitmode == BLOSC_NEVER_SPLIT) {
        category = 24;
    } else if (codec == BLOSC_ZLIB && filter == BLOSC_SHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 25;
    } else if (codec == BLOSC_ZLIB && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 26;
    } else if (codec == BLOSC_ZSTD && filter == BLOSC_NOFILTER && splitmode == BLOSC_NEVER_SPLIT) {
        category = 27;
    } else if (codec == BLOSC_ZSTD && filter == BLOSC_SHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 28;
    } else if (codec == BLOSC_ZSTD && filter == BLOSC_BITSHUFFLE && splitmode == BLOSC_NEVER_SPLIT) {
        category = 29;
    }

    return category;
}

void print_compress_info(void) {
  char *name = NULL, *version = NULL;
  int ret;

  printf("Blosc version: %s (%s)\n", BLOSC_VERSION_STRING, BLOSC_VERSION_DATE);

  printf("List of supported compressors in this build: %s\n",
         blosc_list_compressors());

  printf("Supported compression libraries:\n");
  ret = blosc_get_complib_info("blosclz", &name, &version);
  if (ret >= 0) printf("  %s: %s\n", name, version);
  free(name);
  free(version);
  ret = blosc_get_complib_info("lz4", &name, &version);
  if (ret >= 0) printf("  %s: %s\n", name, version);
  free(name);
  free(version);
  ret = blosc_get_complib_info("zlib", &name, &version);
  if (ret >= 0) printf("  %s: %s\n", name, version);
  free(name);
  free(version);
  ret = blosc_get_complib_info("zstd", &name, &version);
  if (ret >= 0) printf("  %s: %s\n", name, version);
  free(name);
  free(version);
}


int main(int argc, char *argv[]) {
  char usage[256];
  char csv_filename[256];
  char data_filename[256];
  bool entropy_probe_mode = false;

  print_compress_info();

  strcpy(usage, "Usage: entropy_probe [-e] data_filename");

  if (argc < 2) {
    printf("%s\n", usage);
    exit(1);
  }

  if (argc >= 3) {
    if (strcmp("-e", argv[1]) != 0) {
      printf("%s\n", usage);
      exit(1);
    }
    strcpy(data_filename, argv[2]);
    entropy_probe_mode = true;
  }
  else {
    strcpy(data_filename, argv[1]);
  }
  printf("fitxer %s\n", data_filename);
  blosc2_cparams cparams = BLOSC2_CPARAMS_DEFAULTS;
  cparams.instr_codec = true;
  cparams.blocksize = 256 * (int)KB;
  cparams.splitmode = BLOSC_ALWAYS_SPLIT;

  if (entropy_probe_mode) {
    // The entropy probe detector is meant to always be used in SPLIT mode
    cparams.splitmode = BLOSC_ALWAYS_SPLIT;
  }
  cparams.typesize = sizeof(float);

  if (entropy_probe_mode) {
    blosc2_codec udcodec;
    udcodec.compcode = ENTROPY_PROBE_ID;
    udcodec.compver = 1;
    udcodec.complib = 1;
    udcodec.compname = "entropy_probe";
    udcodec.encoder = entropy_probe;
    udcodec.decoder = NULL;
    blosc2_register_codec(&udcodec);
    cparams.compcode = ENTROPY_PROBE_ID;
  } //else {
  //  cparams.compcode = BLOSC_BLOSCLZ;
  //}

  blosc2_dparams dparams = BLOSC2_DPARAMS_DEFAULTS;

  blosc2_schunk *schunk = blosc2_schunk_open(data_filename);
  if (schunk == NULL) {
    printf("Cannot open the data file\n");
    exit(1);
  }
  uint8_t *chunk = malloc(schunk->chunksize);
  uint8_t *chunk2 = malloc(schunk->chunksize);
  printf("nchunks in dataset: %d\n", schunk->nchunks);

  if (!entropy_probe_mode) {
      // Loop over different filters and codecs
      int v_codecs[] = {0, 1, 2, 4, 5};
      int v_splits[] = {1, 2};

      for (int ncodec = 0; ncodec < sizeof(v_codecs); ++ncodec) {
          for (int nfilter = 0; nfilter <= BLOSC_BITSHUFFLE; ++nfilter) {
              for (int nsplit = 0; nsplit < BLOSC_NEVER_SPLIT; ++nsplit) {

                  cparams.splitmode = v_splits[nsplit];
                  cparams.compcode = v_codecs[ncodec];
                  cparams.filters[BLOSC2_MAX_FILTERS - 1] = nfilter;
                  blosc2_context *cctx = blosc2_create_cctx(cparams);
                  blosc2_context *dctx = blosc2_create_dctx(dparams);

                  const char *compname;
                  switch (cparams.compcode) {
                      case BLOSC_BLOSCLZ:
                          compname = BLOSC_BLOSCLZ_COMPNAME;
                          break;
                      case BLOSC_LZ4:
                          compname = BLOSC_LZ4_COMPNAME;
                          break;
                      case BLOSC_LZ4HC:
                          compname = BLOSC_LZ4HC_COMPNAME;
                          break;
                      case BLOSC_ZLIB:
                          compname = BLOSC_ZLIB_COMPNAME;
                          break;
                      case BLOSC_ZSTD:
                          compname = BLOSC_ZSTD_COMPNAME;
                          break;
                      case ENTROPY_PROBE_ID:
                          compname = "entropy";
                          break;
                      default:
                          printf("Unsupported codec!");
                          exit(1);
                  }

                  char *sfilter;
                  switch (nfilter) {
                      case BLOSC_NOFILTER:
                          sfilter = "nofilter";
                          break;
                      case BLOSC_SHUFFLE:
                          sfilter = "shuffle";
                          break;
                      case BLOSC_BITSHUFFLE:
                          sfilter = "bitshuffle";
                          break;
                      default:
                          printf("Unsupported filter!");
                          exit(1);
                  }

                  char *ssplit;
                  switch (cparams.splitmode) {
                      case BLOSC_ALWAYS_SPLIT:
                          ssplit = "split";
                          break;
                      case BLOSC_NEVER_SPLIT:
                          ssplit = "nosplit";
                          break;

                      default:
                          printf("Unsupported splitmode!");
                          exit(1);
                  }

                  // Create csv file
                  sprintf(csv_filename, "%s-%s-%s.csv", compname, sfilter, ssplit);
                  printf("EEEEEEEEEEEEEEEEEE el nombre del codec es: %s\n", compname);
                  printf("EEEEEEEEEEEEEEEEEE el nombre del filtre es: %s\n", sfilter);
                  printf("EEEEEEEEEEEEEEEEEE el nombre del split es: %s\n", ssplit);

                  FILE *csv_file = fopen(csv_filename, "w");
                  if (csv_file == NULL) {
                      printf("Error creating the file\n");
                      return -1;
                  }

                  int category;
                  if (strcmp(compname, BLOSC_BLOSCLZ_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 0;
                  } else if (strcmp(compname, BLOSC_BLOSCLZ_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 1;
                  } else if (strcmp(compname, BLOSC_BLOSCLZ_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 2;
                  } else if (strcmp(compname, BLOSC_LZ4_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 3;
                  } else if (strcmp(compname, BLOSC_LZ4_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 4;
                  } else if (strcmp(compname, BLOSC_LZ4_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 5;
                  } else if (strcmp(compname, BLOSC_LZ4HC_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 6;
                  } else if (strcmp(compname, BLOSC_LZ4HC_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 7;
                  } else if (strcmp(compname, BLOSC_LZ4HC_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 8;
                  } else if (strcmp(compname, BLOSC_ZLIB_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 9;
                  } else if (strcmp(compname, BLOSC_ZLIB_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "split") == 00) {
                      category = 10;
                  } else if (strcmp(compname, BLOSC_ZLIB_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 11;
                  } else if (strcmp(compname, BLOSC_ZSTD_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 12;
                  } else if (strcmp(compname, BLOSC_ZSTD_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 13;
                  } else if (strcmp(compname, BLOSC_ZSTD_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "split") == 0) {
                      category = 14;

                  } else if (strcmp(compname, BLOSC_BLOSCLZ_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 15;
                  } else if (strcmp(compname, BLOSC_BLOSCLZ_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 16;
                  } else if (strcmp(compname, BLOSC_BLOSCLZ_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 17;
                  } else if (strcmp(compname, BLOSC_LZ4_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 18;
                  } else if (strcmp(compname, BLOSC_LZ4_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 19;
                  } else if (strcmp(compname, BLOSC_LZ4_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 20;
                  } else if (strcmp(compname, BLOSC_LZ4HC_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 21;
                  } else if (strcmp(compname, BLOSC_LZ4HC_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 22;
                  } else if (strcmp(compname, BLOSC_LZ4HC_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 23;
                  } else if (strcmp(compname, BLOSC_ZLIB_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 24;
                  } else if (strcmp(compname, BLOSC_ZLIB_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 25;
                  } else if (strcmp(compname, BLOSC_ZLIB_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 26;
                  } else if (strcmp(compname, BLOSC_ZSTD_COMPNAME) == 0 && strcmp(sfilter, "nofilter") == 0 && strcmp(ssplit, "nosplit") == 00) {
                      category = 27;
                  } else if (strcmp(compname, BLOSC_ZSTD_COMPNAME) == 0 && strcmp(sfilter, "shuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 28;
                  } else if (strcmp(compname, BLOSC_ZSTD_COMPNAME) == 0 && strcmp(sfilter, "bitshuffle") == 0 && strcmp(ssplit, "nosplit") == 0) {
                      category = 29;
                  }

                  fprintf(csv_file, "cratio, speed, special_vals, nchunk, category\n");
                  for (int nchunk = 0; nchunk < schunk->nchunks; nchunk++) {
                      printf("decompressing chunk # %d (out of %d)\n", nchunk, schunk->nchunks);
                      int dsize = blosc2_schunk_decompress_chunk(schunk, nchunk, chunk, schunk->chunksize);
                      if (dsize < 0) {
                          printf("Error decompressing chunk in schunk.  Error code: %d\n", dsize);
                          return dsize;
                      }
                      int csize = blosc2_compress_ctx(cctx, chunk, dsize, chunk2, schunk->chunksize);
                      if (csize < 0) {
                          printf("Error compressing chunk.  Error code: %d\n", csize);
                          return csize;
                      }
                      int dsize2 = blosc2_decompress_ctx(dctx, chunk2, csize, chunk, dsize);
                      if (dsize2 < 0) {
                          printf("Error decompressing chunk.  Error code: %d\n", dsize2);
                          return dsize2;
                      }

                      int nstreams = dsize2 / (int) sizeof(blosc2_instr);
                      printf("Chunk %d data with %d streams:\n", nchunk, nstreams);

                      int categoria = categorize_data(compname, *sfilter, *ssplit);
                      extr_data(csv_file, nchunk, chunk, nstreams, category);
                  }

                  fclose(csv_file);
              }
          }
      }

  } else {
      // Loop over different filters
      for (int nfilter = 0; nfilter <= BLOSC_BITSHUFFLE; ++nfilter) {
          cparams.filters[BLOSC2_MAX_FILTERS - 1] = nfilter;
          blosc2_context *cctx = blosc2_create_cctx(cparams);
          blosc2_context *dctx = blosc2_create_dctx(dparams);

          const char *compname;
          switch (cparams.compcode) {
              case BLOSC_BLOSCLZ:
                  compname = BLOSC_BLOSCLZ_COMPNAME;
                  break;
              case BLOSC_LZ4:
                  compname = BLOSC_LZ4_COMPNAME;
                  break;
              case BLOSC_LZ4HC:
                  compname = BLOSC_LZ4HC_COMPNAME;
                  break;
              case BLOSC_ZLIB:
                  compname = BLOSC_ZLIB_COMPNAME;
                  break;
              case BLOSC_ZSTD:
                  compname = BLOSC_ZSTD_COMPNAME;
                  break;
              case ENTROPY_PROBE_ID:
                  compname = "entropy";
                  break;
              default:
                  printf("Unsupported codec!");
                  exit(1);
          }

          char *sfilter;
          switch (nfilter) {
              case BLOSC_NOFILTER:
                  sfilter = "nofilter";
                  break;
              case BLOSC_SHUFFLE:
                  sfilter = "shuffle";
                  break;
              case BLOSC_BITSHUFFLE:
                  sfilter = "bitshuffle";
                  break;
              default:
                  printf("Unsupported filter!");
                  exit(1);
          }

          // Create csv file
          sprintf(csv_filename, "%s-%s.csv", compname, sfilter);
          FILE *csv_file = fopen(csv_filename, "w");
          if (csv_file == NULL) {
              printf("Error creating the file\n");
              return -1;
          }

          fprintf(csv_file, "cratio, speed, special_vals, nchunk, category\n");
          for (int nchunk = 0; nchunk < schunk->nchunks; nchunk++) {
              printf("decompressing chunk # %d (out of %d)\n", nchunk, schunk->nchunks);
              int dsize = blosc2_schunk_decompress_chunk(schunk, nchunk, chunk, schunk->chunksize);
              if (dsize < 0) {
                  printf("Error decompressing chunk in schunk.  Error code: %d\n", dsize);
                  return dsize;
              }
              int csize = blosc2_compress_ctx(cctx, chunk, dsize, chunk2, schunk->chunksize);
              if (csize < 0) {
                  printf("Error compressing chunk.  Error code: %d\n", csize);
                  return csize;
              }
              int dsize2 = blosc2_decompress_ctx(dctx, chunk2, csize, chunk, dsize);
              if (dsize2 < 0) {
                  printf("Error decompressing chunk.  Error code: %d\n", dsize2);
                  return dsize2;
              }

              int nstreams = dsize2 / (int) sizeof(blosc2_instr);
              printf("Chunk %d data with %d streams:\n", nchunk, nstreams);

              int category = -1;
              extr_data(csv_file, nchunk, chunk, nstreams, category);
          }

          fclose(csv_file);

      }
  }

  /* Free resources */
  blosc2_schunk_free(schunk);
  free(chunk);
  free(chunk2);
  printf("Success!\n");

  return 0;
}
