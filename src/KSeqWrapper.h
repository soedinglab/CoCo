/*
 * Copyright (C) 2017 Martin Steinegger
 *
 * This file is part of MMSeqs2
 * https://github.com/soedinglab/MMseqs2
 */

#ifndef MMSEQS_KSEQWRAPPER_H
#define MMSEQS_KSEQWRAPPER_H

#include "kseq.h"

#include <string>

class KSeqWrapper {
public:
    struct KSeqEntry {
        kstring_t name;
        kstring_t sequence;
        kstring_t comment;
        kstring_t qual;
    } entry;

    enum kseq_type {
        KSEQ_FILE,
        KSEQ_STREAM,
        KSEQ_GZIP,
        KSEQ_BZIP,
        KSEQ_BUFFER
    };
    kseq_type type;

    virtual bool ReadEntry() = 0;
    virtual ~KSeqWrapper() {};

protected:
    void* seq;
};

class KSeqFile : public KSeqWrapper {
public:
    KSeqFile(const char* file);
    bool ReadEntry();
    ~KSeqFile();
private:
    FILE* file;
};


class KSeqStream : public KSeqWrapper {
public:
    KSeqStream();
    bool ReadEntry();
    ~KSeqStream();
};

#ifdef HAVE_ZLIB
#include <zlib.h>

class KSeqGzip : public KSeqWrapper {
public:
    KSeqGzip(const char* file);
    bool ReadEntry();
    ~KSeqGzip();
private:
    gzFile file;
};
#endif

#ifdef HAVE_BZLIB
#include <bzlib.h>

class KSeqBzip : public KSeqWrapper {
public:
    KSeqBzip(const char* file);
    bool ReadEntry();
    ~KSeqBzip();
private:
    BZFILE *file;
};
#endif


KSeqWrapper* KSeqFactory(const char* file);

#endif //MMSEQS_KSEQWRAPPER_H
