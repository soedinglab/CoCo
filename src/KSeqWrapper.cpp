/*
 * Copyright (C) 2017 Martin Steinegger
 *
 * This file is part of MMSeqs2
 * https://github.com/soedinglab/MMseqs2
 */

#include "KSeqWrapper.h"
#include "kseq.h"
#include "filehandling.h"
#include "Info.h"
#include "util.h"
#include <unistd.h>

namespace KSEQFILE {
    KSEQ_INIT(int, read)
}

KSeqFile::KSeqFile(const char* fileName) {
    file = openFileOrDie(fileName, "r");
    seq = (void*) KSEQFILE::kseq_init(fileno(file));
    type = KSEQ_FILE;
}

bool KSeqFile::ReadEntry() {
    KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*) seq;
    int result = KSEQFILE::kseq_read(s);
    if (result < 0)
        return false;
    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqFile::~KSeqFile() {
    kseq_destroy((KSEQFILE::kseq_t*)seq);
    if (fclose(file) != 0) {
        Info(Info::ERROR) << "Cannot close KSeq input file\n";
        EXIT(EXIT_FAILURE);
    }
}


namespace KSEQSTREAM {
    KSEQ_INIT(int, read)
}

KSeqStream::KSeqStream() {
    seq = (void*) KSEQSTREAM::kseq_init(STDIN_FILENO);
    type = KSEQ_STREAM;
}

bool KSeqStream::ReadEntry() {
    KSEQSTREAM::kseq_t* s = (KSEQSTREAM::kseq_t*) seq;
    int result = KSEQSTREAM::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqStream::~KSeqStream() {
    kseq_destroy((KSEQSTREAM::kseq_t*)seq);
}

#ifdef HAVE_ZLIB
namespace KSEQGZIP {
    KSEQ_INIT(gzFile, gzread)
}

KSeqGzip::KSeqGzip(const char* fileName) {
    if(fileExists(fileName) == false) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }

    file = gzopen(fileName, "r");
    if(file == NULL) {
        perror(fileName); EXIT(EXIT_FAILURE);
    }

    seq = (void*) KSEQGZIP::kseq_init(file);
    type = KSEQ_GZIP;
}

bool KSeqGzip::ReadEntry() {
    KSEQGZIP::kseq_t* s = (KSEQGZIP::kseq_t*) seq;
    int result = KSEQGZIP::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;
    return true;
}

KSeqGzip::~KSeqGzip() {
    kseq_destroy((KSEQGZIP::kseq_t*)seq);
    gzclose(file);
}
#endif


#ifdef HAVE_BZLIB
namespace KSEQBZIP {
    KSEQ_INIT(BZFILE *, BZ2_bzread)
}

KSeqBzip::KSeqBzip(const char* fileName) {
    if(fileExists(fileName) == false) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }
    FILE *fp = openFileOrDie(fileName, "r+b");
    int bzError;
    file = BZ2_bzReadOpen(&bzError, fp, 0, 0, NULL, 0);
    if(bzError != 0){
        perror(fileName); EXIT(EXIT_FAILURE);
    }
    seq = (void*) KSEQBZIP::kseq_init(file);
    type = KSEQ_BZIP;
}

bool KSeqBzip::ReadEntry() {
    KSEQBZIP::kseq_t* s = (KSEQBZIP::kseq_t*) seq;
    int result = KSEQBZIP::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;
    return true;
}

KSeqBzip::~KSeqBzip() {
    kseq_destroy((KSEQBZIP::kseq_t*)seq);
    int bzError;
    BZ2_bzReadClose(&bzError, file);
}
#endif

KSeqWrapper* KSeqFactory(const char* file) {
    KSeqWrapper* kseq = NULL;
    if( strcmp(file, "stdin") == 0 ){
        kseq = new KSeqStream();
        return kseq;
    }

    if(endsWith(".gz", file) == false && endsWith(".bz2", file) == false ) {
        kseq = new KSeqFile(file);
        return kseq;
    }
#ifdef HAVE_ZLIB
    else if(endsWith(".gz", file) == true) {
        kseq = new KSeqGzip(file);
        return kseq;
    }
#else
    else if(Util::endsWith(".gz", file) == true) {
        Info(Info::ERROR) << "CoCo was not compiled with zlib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

#ifdef HAVE_BZLIB
    else if(endsWith(".bz2", file) == true) {
        kseq = new KSeqBzip(file);
        return kseq;
    }
#else
    else if(endsWith(".bz2", file) == true) {
        Info(Info::ERROR) << "CoCo was not compiled with bz2lib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

    return kseq;
}

