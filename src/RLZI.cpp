/*
 *  This file is part of BWTIL.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   BWTIL is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   BWTIL is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

//============================================================================
// Name        : sFM-index.cpp
// Author      : Nicola Prezza
// Version     : 1.0
// Copyright   : GNU General Public License (http://www.gnu.org/copyleft/gpl.html)
// Description : test class for the succinct FM index
//============================================================================

#include <iostream>
#include <queue>
#include <string>
#include <boost/heap/pairing_heap.hpp>
#include <bits/stdc++.h>

#include "BWTIL/data_structures/succinctFMIndex.h"
// #include "BWTIL/data_structures/RMaxQBlockDecomp.h"

#include <sdsl/rmq_support.hpp>
#include <sdsl/rmq_support_sparse_table.hpp>
#include <sdsl/bit_vectors.hpp>
#include <vector>

#include <chrono>
#include <assert.h>
#include <limits>

using namespace bwtil;
using namespace std;
using namespace sdsl;

#include <string>
#include <chrono>
#include <algorithm>
#include <fstream>

#include <thread>
#include "optparse.h"

// Heng Li - KSEQ LIB for parsing FASTA files
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "kstring.h"
KSEQ_INIT(gzFile, gzread)

struct ProfileResult
{
    std::string Name;
    long long Start, End;
};

struct InstrumentationSession
{
    std::string Name;
};

class Instrumentor
{
private:
    InstrumentationSession* m_CurrentSession;
    std::ofstream m_OutputStream;
    int m_ProfileCount;
public:
    Instrumentor()
        : m_CurrentSession(nullptr), m_ProfileCount(0)
    {
    }

    void BeginSession(const std::string& name, const std::string& filepath = "results.json")
    {
        m_OutputStream.open(filepath);
        WriteHeader();
        m_CurrentSession = new InstrumentationSession{ name };
    }

    void EndSession()
    {
        WriteFooter();
        m_OutputStream.close();
        delete m_CurrentSession;
        m_CurrentSession = nullptr;
        m_ProfileCount = 0;
    }

    void WriteProfile(const ProfileResult& result)
    {
        if (m_ProfileCount++ > 0)
            m_OutputStream << ",";

        std::string name = result.Name;
        std::replace(name.begin(), name.end(), '"', '\'');

        m_OutputStream << "{";
        m_OutputStream << "\"cat\":\"function\",";
        m_OutputStream << "\"dur\":" << (result.End - result.Start) << ',';
        m_OutputStream << "\"name\":\"" << name << "\",";
        m_OutputStream << "\"ph\":\"X\",";
        m_OutputStream << "\"pid\":0,";
        m_OutputStream << "\"tid\":0,";
        m_OutputStream << "\"ts\":" << result.Start;
        m_OutputStream << "}";

        m_OutputStream.flush();
    }

    void WriteHeader()
    {
        // m_OutputStream << "{\"otherData\": {},\"traceEvents\":[";
        m_OutputStream << "[";
        m_OutputStream.flush();
    }

    void WriteFooter()
    {
        // m_OutputStream << "]}";
        m_OutputStream << "]";
        m_OutputStream.flush();
    }

    static Instrumentor& Get()
    {
        static Instrumentor instance;
        return instance;
    }
};

class InstrumentationTimer
{
public:
    InstrumentationTimer(const char* name)
        : m_Name(name), m_Stopped(false)
    {
        m_StartTimepoint = std::chrono::high_resolution_clock::now();
    }

    ~InstrumentationTimer()
    {
        if (!m_Stopped)
            Stop();
    }

    void Stop()
    {
        auto endTimepoint = std::chrono::high_resolution_clock::now();

        long long start = std::chrono::time_point_cast<std::chrono::nanoseconds>(m_StartTimepoint).time_since_epoch().count();
        long long end = std::chrono::time_point_cast<std::chrono::nanoseconds>(endTimepoint).time_since_epoch().count();


        if ((end - start) > 0){
          Instrumentor::Get().WriteProfile({ m_Name, start, end});
        }
        m_Stopped = true;
    }
private:
    const char* m_Name;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
    bool m_Stopped;
};

class Timer
{
public:
  Timer()
  {
    m_StartTimepoint = std::chrono::high_resolution_clock::now();
  }

  ~Timer()
  {
    Stop();
  }

  void Stop()
  {
    auto endTimepoint = std::chrono::high_resolution_clock::now();

    auto start = std::chrono::time_point_cast<std::chrono::nanoseconds>(m_StartTimepoint).time_since_epoch().count();

    auto end = std::chrono::time_point_cast<std::chrono::nanoseconds>(endTimepoint).time_since_epoch().count();

    auto duration = end - start;
    double ms = duration * 0.001;

    std::cout << duration << "us (" << ms << ")\n";

  }
private:
  std::chrono::time_point< std::chrono::high_resolution_clock> m_StartTimepoint;
};


char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static inline int aln_score(const int m, const int o) {
	// return m*p->mm_score + o*p->gapo_score + e*p->gape_score;
  return m*3 - o*1;
}

struct query_bases_t {
  uchar b;
  ulint L;
  ulint U;
  ulint maxAr;

  query_bases_t(){}
  query_bases_t(uchar b, ulint L, ulint U, ulint maxAr){}
};

struct aln_entry_t {
	// bwtint_t L;
	// bwtint_t U; // (L,U): SA interval of [i,n-1]
  ulint i;
  ulint L;
  ulint U;
  ulint next_p;
	uint32_t num_mm, num_gapo, num_snps;
	uint32_t score, state, aln_length;
  uchar c;
  uchar b;
  std::vector<int> mis_pos;
  std::vector<char> mis_bases;
	//uint8_t score; // aln score so far
	//uint8_t i; // position in the read
	//uint8_t num_mm;
	//uint8_t num_gapo;
	//uint8_t num_gape;
	//uint8_t state;
	//int n_seed_mm;
	//uint8_t last_diff_pos;
	//uint8_t padding;

	// edit transcript
	//uint8_t aln_length;
	// char aln_path[ALN_PATH_ALLOC];

  aln_entry_t(){} // DEFAULT CONSTRUCTOR

  aln_entry_t(ulint i, ulint L, ulint U, uint32_t num_mm, uint32_t num_gapo, uint32_t next_p, uint32_t state, uint32_t num_snps, uint32_t aln_length, uchar& uu, uchar& uuu): i(i), L(L), U(U), num_mm(num_mm), num_gapo(num_gapo), next_p(next_p), state(state), num_snps(num_snps), aln_length(aln_length), c(uu), b(uuu)
  {
    mis_pos = {};
    mis_bases = {};
    // c = *uu;
    // c = strcpy(c, uu);
  }

  aln_entry_t(ulint i, ulint L, ulint U, uint32_t num_mm, uint32_t num_gapo, uint32_t next_p, uint32_t state, uint32_t num_snps, uint32_t aln_length, uchar& uu, uchar& uuu, std::vector<int>& mis_pos, std::vector<char>& mis_bases): i(i), L(L), U(U), num_mm(num_mm), num_gapo(num_gapo), next_p(next_p), state(state), num_snps(num_snps), aln_length(aln_length), c(uu), b(uuu), mis_pos(mis_pos), mis_bases(mis_bases)
  {
    // c = *uu;
    // c = strcpy(c, uu);
  }

};


struct CompareAlnScores {
  bool operator()(aln_entry_t const& p1, aln_entry_t const& p2){
    // return true if alignment p1 has lower score
    // than alignment p2
    return aln_score(p1.num_mm, 0) < aln_score(p2.num_mm, 0);
  }
};

// Function which return string by concatenating it.
string repeat(string s, int n)
{
    // Copying given string to temparory string.
    string s1 = s;

    for (int i=1; i<n;i++)
        s += s1; // Concatinating strings

    return s;
}

#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)

int main(int argc,char** argv) {

  optparse::OptionParser parser =         optparse::OptionParser().description(R"(
  _____  _      ___________
 |  __ \| |    |___  /_   _|
 | |__) | |       / /  | |
 |  _  /| |      / /   | |
 | | \ \| |____ / /__ _| |_
 |_|  \_\______/_____|_____|


RLZI - Relative Lempel-Ziv with Inexact matchings)").usage("usage: %prog [options] -R <ref.fasta> -S <samp.fasta>").version("%prog 1.0");

  parser.add_option("-k", "--seed").dest("seedlen")
          .type("int")
          .set_default("10")
          .help("Minimum seed length.");

  parser.add_option("-m", "--num-mismatches").dest("nummis")
          .type("int")
          .set_default("2")
          .help("Report phrases with at most <int> mismatches.");

  parser.add_option("-R", "--referenceseq").dest("refseq")
          .help("Fasta or multi-fasta with the reference sequence").metavar("FILE");

  parser.add_option("-S", "--sampleseq").dest("sampseq")
          .help("Fasta with the single sequence that we will parse").metavar("FILE");

  parser.add_option("-o", "--outprefix").dest("outprefix")
          .set_default("out")
          .help("Prefix to be used in outputs").metavar("FILE");

  parser.add_option("-v", "--verbose")
          .action("store_true")
          .dest("verbose")
          .set_default("0")
          .help("don't print status messages to stdout");

  parser.add_option("-b", "--readable")
          .action("store_true")
          .dest("readable")
          .set_default("0")
          .help("output auxiliary bitvectors such as B, I and Ch in txt format (it requires a lot of space if parsing genomes greater than 1 GB)");

  parser.add_option("-c", "--reverse-complement").dest("rc")
          .action("store_true")
          .set_default("0")
          .help("Include the reverse-complement of the reference sequence");

    const optparse::Values options = parser.parse_args(argc, argv);
    const std::vector<std::string> args = parser.args();


    // if (options["referenceseq"].size() == 0 | options["sampleseq"].size() == 0){
        // parser.error("You didn't specify a reference or sample sequence.");
    // }

    std::cout << "VERBOSE OPTION= " << options["verbose"] << std::endl;
    if (options.get("verbose"))
    {
        std::cout << options["filename"] << "\n";
    }

	// if(argc != 4 and argc != 3 and argc != 5){
	// 	cout << "*** succinct FM-index data structure : a wavelet-tree based uncompressed FM index ***\n";
	// 	cout << "Usage: sFM-index option file [pattern]\n";
	// 	cout << "where:\n";
	// 	cout <<	"- option = build|search|lz. \n";
	// 	cout << "- file = path of the text file (if build mode) or .sfm sFM-index file (if search mode). \n";
	// 	cout << "- pattern = must be specified in search mode. It is the pattern to be searched in the index.\n";
	// 	exit(0);
	// }

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

	int build=0,search=1,lz=2;
  int readablefiles=options.get("readable");
	int mode;

	// if(string(argv[1]).compare("build")==0)
	// 	mode=build;
	// else if(string(argv[1]).compare("search")==0)
	// 	mode=search;
	// else if(string(argv[1]).compare("lz")==0)
	// 	mode=lz;
	// else{
	// 	cout << "Unrecognized option "<<argv[1]<<endl;
	// 	exit(0);
	// }
  mode = 2;

  string pref(options["outprefix"]);

	string in(options["refseq"]);
	string inS(options["sampseq"]);

  string out(in);
	out.append(".sfm");

  string outS(inS);
	outS.append(".sfm");

	string pattern;
  //
	// if(mode==search){
	// 	pattern = string(argv[3]);
	// }

  int k_seed;
  int debug = options.get("verbose");


  // if(mode==lz){
  //   k_seed = atoi(options.get("num-mismatches"));
  //   inS = argv[4];
  // }

    auto t1 = high_resolution_clock::now();

	succinctFMIndex SFMI;
  // RL
  succinctFMIndex SFMI_S;

	// if(mode==build){

  int rc = atoi(options.get("rc"));
  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(in.c_str(), "r");
  seq = kseq_init(fp);
  char *s;
  string ref_string;
  ulint refl;
  refl = 0;
  vector<ulint> chrom_lens;
  while ((l = kseq_read(seq)) >= 0){
    s = seq->seq.s;
    cout << "seq: \n" <<  seq->name.s << endl;
    cout << "len: " << seq->seq.l << endl;
    ref_string = ref_string + static_cast<string>(s);
    ref_string.append("$");
    chrom_lens.push_back(refl);
    refl = refl + seq->seq.l + 1;
    if(rc){
      int c0, c1;
			for (int i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
				c0 = comp_tab[(int)seq->seq.s[i]];
				c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
				seq->seq.s[i] = c1;
				seq->seq.s[seq->seq.l - 1 - i] = c0;
			}
			if (seq->seq.l & 1){ // complement the remaining base
				seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];
      }
      ref_string = ref_string + static_cast<string>(seq->seq.s);
      ref_string.append("$");
      chrom_lens.push_back(refl);
      refl = refl + seq->seq.l + 1;
    }
  }
  // cout << ref_string << endl;
  // cout << ref_string.size() << endl;
  // cout << "RC=" << rc << endl;
  // ref_string = static_cast<string>(s);

      kseq_destroy(seq);
      gzclose(fp);
  cout << "len total: " << refl << endl;

  bit_vector Chref(refl-1,0);

  for(int u=0; u < chrom_lens.size(); u++){
    cout << "LEN chrom[" << u << "]=" << chrom_lens << endl;
    Chref[chrom_lens[u]] = 1;
  }

  chrom_lens.push_back(refl);

  cout << "Storing Chref bitvector... " << endl;
  sd_vector<> sdbCh(Chref);
  store_to_file(sdbCh, pref+"_sdbChref.sdsl");

  if (readablefiles){
    ofstream outfileCh;
    outfileCh.open("outChref.txt");
    outfileCh << Chref;
    outfileCh.close();
  }


  sd_vector<>::rank_1_type Chref_rank(&sdbCh);

  // del Chref;

  // sd_vector<>::rank_1_type Ch_rank(&Chref);

    std::reverse(ref_string.begin(), ref_string.end());


    ifstream frefbwt(out.c_str());

    if (frefbwt.good() == 0){
      cout << "Building succinct FM-index of REF file "<< in << endl;
		    SFMI = succinctFMIndex(ref_string,refl,false);

		    cout << "\nStoring succinct FM-index of REF in "<< out << endl;
		    SFMI.saveToFile(out);
    }
    gzFile fpS;
    kseq_t *seqS;
    int ls;
    fpS = gzopen(inS.c_str(), "r");
    seqS = kseq_init(fpS);
    char *sS;
    string samp_string;
    ulint sampl;
    sampl = 0;
    while ((ls = kseq_read(seqS)) >= 0){
      sS = seqS->seq.s;
      cout << "seq: \n" <<  seqS->name.s << endl;
      cout << "len: " << seqS->seq.l << endl;
      samp_string = static_cast<string>(sS);
      samp_string.append("$");
      sampl = seqS->seq.l + 1;

    }
    // ref_string = static_cast<string>(s);
    cout << "len total: " << sampl << endl;
    kseq_destroy(seqS);
    gzclose(fpS);


    std::reverse(samp_string.begin(), samp_string.end());
    cout << "Building succinct FM-index of SAMPLE file "<< inS << endl;
		SFMI = succinctFMIndex(samp_string,sampl,false);

		cout << "\nStoring succinct FM-index of SAMPLE in "<< outS << endl;
		SFMI.saveToFile(outS);

		cout << "Done.\n";

	// }

	if(mode==search){

		cout << "Loading succinct FM-index from file "<< in <<endl;
		SFMI = succinctFMIndex::loadFromFile(in);
		cout << "Done." << endl;

		cout << "\nSearching pattern \""<< pattern << "\""<<endl;

		vector<ulint> occ = SFMI.getOccurrencies( pattern );

		cout << "The pattern occurs " << occ.size() << " times in the text at the following positions : \n";

		for(uint i=0;i<occ.size();i++)
			cout << occ.at(i) << " ";

		cout << "\n\nDone.\n";

	}

	if (mode==lz){


    k_seed = atoi(options.get("seedlen")) - 1;

    Instrumentor::Get().BeginSession("Profile");

    ofstream outfile;
    string outfilename;
    outfilename = pref + "_sdsl.txt";
    outfile.open(outfilename.c_str());

		cout << "Loading succinct FM-index from file " << in << endl;
		SFMI = succinctFMIndex::loadFromFile(out);
		cout << "Done." << endl;
		IndexedBWT* idxBWT;
    cout << "Time SFMI pointer" << endl;
    // {
    //     InstrumentationTimer timer("SFMI pointer");
        idxBWT = SFMI.get_idxBWTPtr();
    // }

    // RMaxQBlockDecomp RMQ(idxBWT);
    // RLZ - start
    SFMI_S = succinctFMIndex::loadFromFile(outS);
    IndexedBWT* idxBWTS = SFMI_S.get_idxBWTPtr();

		// vector<factor> enc;

		ulint i = 0;
		uint j = 0;
		ulint p = 1;
		uchar c;
		ulint cidx;
		pair<ulint, ulint> interval;
		ulint max_Ar_idx, max_Ar;
		uint iprime;
    ulint pos_factor;
    int which_chr;
    // LZ-parser - take the Text length
		// ulint n = SFMI.textLength();
    // RLZ-parser - take the Sample text length
    ulint n = SFMI_S.textLength();
    ulint nR = SFMI.textLength();

    cout << "Size of bitvectors = " << nR << endl;
    bit_vector B(n,0);
    bit_vector I(n,0);
    // sd_vector<> B(bv);
    std::vector<char> M;


    // string inFileName = "inpChref.txt";
    // ifstream inFile;
    // inFile.open(inFileName.c_str());
    // int countch = 0;
    // while (!inFile.eof()){
    //   inFile >> Chref[countch];
    //   // std::cout << M[count] << " ";
    //   countch = countch + 1;
    // }
    //
    // cout << "Inicio dos cromossomos em R = " << Chref << endl;

    int_vector<> vec(nR);
    for (size_t i=0; i < nR; i++)
      vec[i] = idxBWT->convertToTextCoordinate(i);
    rmq_succinct_sct<false> RMQsdsl(&vec);
    // rmq_succinct_sada<false> RMQsdsl(&vec);
    // rmq_support_sparse_table<int_vector<>,false> RMQsdsl(&vec);
    util::clear(vec);
    // RLZ-parser - END snippet

    cout << "TEXT LENGTH = " << nR << endl;
		ulint sp, ep;
    ulint next_p;
		int perc, last_perc= -1;
		// cout << "Letter " << c << endl;
		priority_queue<aln_entry_t, vector<aln_entry_t>, CompareAlnScores> heap;

    int best_score, num_best, max_entries, total_entries, last_num_entries, alns_get_num_entries;

    int diff_left, best_diff;
    uint32_t matches, max_diff;

    if (debug){
      cout << "###################################################" << endl;
      cout << "INITIALIZING ALGORITHM" << endl;
      cout << "###################################################\n\n" << endl;
    }
    string s = "##";
    int nhashtags;
    int allow_diff;

    int ii, old_p;
    p = 1;
    // for(ii=0; ii <= n; ii++){
    //   old_p = p;
    //   c = idxBWTS->at( p );
    //   p = idxBWTS->LF( p );
    //   cout << "CHAR " << c << " int(CHAR) = " << int(c) << " WITH P = " << old_p <<", NEXT P = " << p << endl;
    // }

    ulint track_len = 0;
    I[0] = 1;

    // int ref_len = 13;
    // int k_seed = 3;
    j = i + 1;
		while (i < (n)){

      if (debug){
        cout << "###################################################" << endl;
        cout << " FOR i = " << i << endl;
      }
      nhashtags = 1;
			// cout << (int)(n/10.) << endl;
			// perc = (100*j)/n;
			// if (perc > last_perc and (perc%5)) {
			// 	cout << " " << perc << "% done." << endl;
			// 	last_perc=perc;
			// }

			// j = i + 1;


      // {
      //       InstrumentationTimer timer("idxBWTS_");
			         c = idxBWTS->at( p );
      // }
			// p = idxBWT->LF( p );

      // {
      //       InstrumentationTimer timer("array_C");
			       interval = SFMI.arrayC(c);
      // }

			// cout << "LETTER =" << c << ", INTERVAL FIRST = " << interval.first << ", SA=" << idxBWT->convertToTextCoordinate(interval.first) << "| SECOND = " << interval.second << ", SA=" << idxBWT->convertToTextCoordinate(interval.second-1) << endl;
			// cout << "RMQ[sp,ep] = " << RMQ.query(interval.first, interval.second-1) << endl;

			sp = interval.first;
			ep = interval.second;

      // {
      //       InstrumentationTimer timer("RMQ");

            // cout << "1) RMQ(" << sp << ", " << ep-1 << ") traditional = " <<  RMQ.query(sp, ep-1) ;
            if (sp <= ep-1){
              if (ep-1 >= nR){
                  if (sp >= nR){
                    // sp = vec.size()-1;
                    max_Ar_idx = RMQsdsl(nR-1, nR-1);
                    max_Ar = idxBWT->convertToTextCoordinate(nR);

                  } else {
                    // ep = vec.size();
                    max_Ar_idx = RMQsdsl(sp, nR-1);
                    max_Ar = idxBWT->convertToTextCoordinate(max_Ar_idx);
                  }
                // cout << " and RMQ-sdsl = " << max_Ar << endl;
              } else {
                max_Ar_idx = RMQsdsl(sp, ep-1);
                max_Ar = idxBWT->convertToTextCoordinate(max_Ar_idx);
                // cout << " and RMQ-sdsl = " << max_Ar << endl;
              }

          } else {
            max_Ar = 0;
            // cout << " and RMQ-sdsl = 0" << endl;
          }

          // assert(RMQ.query(sp, ep-1) == max_Ar);
      // }
			iprime = 0;

			heap = priority_queue<aln_entry_t, vector<aln_entry_t>, CompareAlnScores>();

      // {
      //       InstrumentationTimer timer("heap");
			         heap.push(aln_entry_t(i, sp, ep, 0, 0, p, 0, j+1, 0, c, c));
      // }

      best_score = aln_score(2, 2);
      best_diff = 2;
      max_diff = atoi(options.get("nummis"));
      if (max_diff == 0) {
        k_seed = std::numeric_limits<int>::max();
      } else {
      max_diff = max_diff-1;}
      // cout << "MAAAAX DIFF = " << max_diff << endl;
      num_best = 0;
      max_entries = 0;

      total_entries = 0;
      last_num_entries = 1;


      alns_get_num_entries = 0;

			// cout << "Max_Ar = " << max_Ar << endl;
			// cout << "i= " << i << " j= " << j << " c= " << c << " p= " << p << " sp=" << sp << " ep=" << ep << endl;

			// while ((j + 1 < n) && (sp <= ep) && (max_Ar >= (n-(i+1)))){
      aln_entry_t a;

      // allow_diff = 0;
      while (!heap.empty()){
        // cout << "LEN HEAP = " << heap.size() << endl;
        a = heap.top();
        // {
        //       InstrumentationTimer timer("heap");
              heap.pop();
        // }

        // allow_diff counts the numer of exact matchings
        allow_diff = a.state;

        diff_left = max_diff - a.num_mm;
        sp = a.L;
        ep = a.U;
        matches = a.aln_length;

        c = a.c;
        // {
        //       InstrumentationTimer timer("idxBWTS_LF");
              p = idxBWTS->LF( a.next_p );
        // }

        // {
        //       InstrumentationTimer timer("RMQ");
              // max_Ar = RMQ.query(sp, ep-1);
              // cout << "WTF) RMQ(" << sp << ", " << ep-1 << ") traditional = " <<  RMQ.query(sp, ep-1) ;
              if (sp <= ep-1){
                if (ep-1 >= nR){
                  if (sp >= nR){
                    // sp = vec.size()-1;
                    // max_Ar_idx = RMQsdsl(vec.size()-1, vec.size()-1);
                    max_Ar = idxBWT->convertToTextCoordinate(nR);
                  } else {
                    // ep = vec.size();
                    max_Ar_idx = RMQsdsl(sp, nR-1);
                    max_Ar = idxBWT->convertToTextCoordinate(max_Ar_idx);
                  }
                  // cout << " and RMQ-sdsl = " << max_Ar << endl;
                } else {
                  max_Ar_idx = RMQsdsl(sp, ep-1);
                  max_Ar = idxBWT->convertToTextCoordinate(max_Ar_idx);
                  // cout << " and RMQ-sdsl = " << max_Ar << endl;
                }

            } else {
              max_Ar = 0;
              // cout << " and RMQ-sdsl = 0" << endl;
            }
            // assert(RMQ.query(sp, ep-1) == max_Ar);
        // }

        if (debug){
          cout << "Pop i = " << a.i << " c=" << a.c << " and c =" << c << "int(c) =" << int(a.c) << " , actual base = "<< a.b << ", max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps <<", j new=" << a.i+a.aln_length << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << " num_mm = " << a.num_mm << " n-(a.i+1) =" << n-(a.i+1) << endl;

          // cout << "DUPLICATED Pop i + allow_diff = " << a.i + a.state << " allow_diff = " << a.state << " c=" << a.c << " and c =" << c << " , actual base = "<< a.b << ", max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps <<", j new=" << a.i+a.aln_length+allow_diff << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << " num_mm = " << a.num_mm << endl;

          cout << "Allow diff = " << allow_diff << endl;
        }

        if (diff_left < 0){
          if (debug){
            cout << " INSIDE DIFF LEFT !!!" << diff_left << endl;
          }
	  // j = j - 1;
	  // a.aln_length = a.aln_length - 1;
          continue;}

        if (sp > (ep-1)){
          if (debug){
            cout << " INSIDE SP > (EP-1) !!!" << endl;
          }
          continue;}

        // LZ
        // if (max_Ar < n-(a.i+1)){
        // RLZ
        if (max_Ar < 0){
          if (debug){
          cout << " INSIDE THE MAX AR !!!" << endl;
        }
          continue;
        }

        // LZ
        // if (a.num_snps >= (n-1)){
        // RLZ
        // if (a.num_snps >= (n)){
        //   if (debug){
        //   cout << " INSIDE THE NUM SNPS !!! " << a.num_snps << endl;
        // }
        //   // a.num_snps = j + a.num_mm;
        //
        //   // if (allow_diff >= k_seed){
        //   //   j = a.num_snps;
        //   // }
        //   // j = j + 1;
        //   continue;
        // }

        if ((int(a.c) == 36) || (int(a.c) == 0)){
          if (debug){
            cout << " INSIDE $ BASE DETECTION !!! " << a.c << endl;
          }
          continue;
        }

        // in RLZ we comment all the above code, up to the empty line
        // if (a.next_p == 1) {
        //   if (debug){
        //     cout << " CYCLING THROUGH BWT !! STOP BEFORE GETTING TO FIRST CHAR !!" << endl;
        //   }
        //   break;
        // }

        // cout << "## Pop i = " << a.i << " c=" << a.c << " and c =" << c << " , actual base = "<< a.b << " max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << endl;
        nhashtags ++;
        if (debug){
          // cout << repeat(s, a.aln_length+2) << " ALLOWING INNER MATCHES" << endl;
          cout << "## " << " ALLOWING INNER MATCHES" << endl;
        }
        // LZ - start
				// iprime = n - max_Ar;
        // LZ - end

        // RLZ - begin
        iprime = nR - max_Ar;
        // RLZ - end

        a.num_gapo = iprime;
        j = a.num_snps;
        // c = idxBWT->at( p );
        // {
        //       InstrumentationTimer timer("idxBWTS_");
              c = idxBWTS->at( p );
        // }

        if ((int(c) == 0) || (int(c) == 36)){
          continue;
        }

        if (debug){
        // cout << repeat(s, a.aln_length+2) << " IN c=" << c << " with iprime=" << iprime << " and i=" << a.i+allow_diff << ", p=" << p << endl;
        cout << "### " << " IN c=" << c << " with iprime=" << iprime << " and i=" << a.i+allow_diff << ", p=" << p << endl;
        }

        next_p =  a.next_p;
        // p = a.next_p;

        uchar alphabet[] = {'A', 'C', 'G', 'T'};
        struct query_bases_t qbases[4];

				// j = j + 1;
        int is_mm, i_ref;
        // int allow_diff = 1;
        // allow_diff = allow_diff - a.num_mm;
        // cout << "PROXIMO ALLOW DIFF = " << allow_diff << endl;

        if (allow_diff > k_seed){


          for(int i=0; i < 4; i++){
            // {
            //       InstrumentationTimer timer("idxBWT_exact");
                  interval = idxBWT->exact_match( alphabet[i], a.L, a.U);
            // }
    				sp = interval.first;
    				ep = interval.second;

            qbases[i].b = alphabet[i];
            if (qbases[i].b == c){
              i_ref = i;
            }
            qbases[i].L = sp;
            qbases[i].U = ep;
            // {
            //       InstrumentationTimer timer("RMQ");
                  // qbases[i].maxAr =RMQ.query(sp, ep-1);
                  // cout << "2) RMQ(" << sp << ", " << ep-1 << ") traditional = " <<  RMQ.query(sp, ep-1) ;
                  if (sp <= ep-1){
                    if (ep-1 >= nR){
                      if (sp >= nR){
                        // sp = vec.size()-1;
                        // max_Ar = vec[vec.size()-1];
                        // cout << "%%%%%%%%% HERE WE MAKE A MISTAKE!! " << idxBWT->convertToTextCoordinate(nR-1) << ", " << idxBWT->convertToTextCoordinate(nR) << endl;
                        qbases[i].maxAr = idxBWT->convertToTextCoordinate(nR);

                      } else {
                        // ep = vec.size();
                        max_Ar_idx = RMQsdsl(sp, nR-1);
                        qbases[i].maxAr = idxBWT->convertToTextCoordinate(max_Ar_idx);
                      }
                      // cout << " and RMQ-sdsl = " << max_Ar << endl;
                    } else {
                      max_Ar_idx = RMQsdsl(sp, ep-1);
                      qbases[i].maxAr = idxBWT->convertToTextCoordinate(max_Ar_idx);
                      // cout << " and RMQ-sdsl = " << max_Ar << endl;
                    }

                } else {
                  qbases[i].maxAr = 0;
                  // cout << " and RMQ-sdsl = 0" << endl;
                }
                // assert(RMQ.query(sp, ep-1) == max_Ar);
                // qbases[i].maxAr = max_Ar;
            // }
            // cout << "letter = " << alphabet[i] << " with L = " << sp << " and U = " << ep << " and max Ar=" <<   qbases[i].maxAr << endl;
          }


          if (debug){
            // cout << repeat(s, a.aln_length+2) << "&&& BEGIN MISMATCH STYLE ALLOWED !!!" << endl;
            // cout << repeat(s, a.aln_length+2) << " max Ar = " << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << " allow_diff =" << allow_diff << endl;

            cout << "## "<< "&&& BEGIN MISMATCH STYLE ALLOWED !!!" << endl;
            cout << "## "<< " max Ar = " << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << " allow_diff =" << allow_diff << endl;

          }

            for(int i=0; i < 4; i++){

              if (debug) {
            			cout << "BEFORE CONDITION BEGIN" << endl;
            			// cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " << allow_diff << endl;
                  cout << "## " << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " << allow_diff << endl;
            			cout << "BEFORE CONDITION END" << endl;
            		  }

              // LZ
              // if ((qbases[i].L < qbases[i].U) && (max_Ar >= n-(a.i+1)) ){
              // RLZ
              if ( (qbases[i].L < qbases[i].U) ){

                // cout << "AAAAQUI " << n-(a.i+1)<< endl;
                is_mm = 0;

                if (c != qbases[i].b){
                  is_mm = 1;
                }

                if (is_mm == 0){

                  if (debug){
                    // cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " << allow_diff << endl;
                    cout << "## " << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " <<  allow_diff << " max(n-(i+1), n-(ref_len+1)) = " << n-(a.i+1) << endl;
                  }

                  // LZ
                  // if (qbases[i].maxAr >= n-(a.i+1))
                  // RLZ
                  if (qbases[i].maxAr >= 0)
                  {
                  //       InstrumentationTimer timer("heap");
                        // cout << "PUSH MATCH = " << i << " ep-sp=" << qbases[i].U-qbases[i].L << endl;
                        heap.push(aln_entry_t(a.i, qbases[i].L, qbases[i].U, a.num_mm, (nR-qbases[i].maxAr), p, a.state+1, j+1, a.aln_length+1, c, qbases[i].b, a.mis_pos, a.mis_bases));
                      // }

                }} else {

                  if (debug){
                    // cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm+1 << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << endl;
                    cout << "## " << " Pushed base " << qbases[i].b << " with " << a.num_mm+1 << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << " max(n-(i+1), n-(ref_len+1)) = " << std::max(n-(a.i+1), n-(13+1)) << endl;
                  }

                  // LZ
                  // if (qbases[i].maxAr >= n-(a.i+1))
                  // RLZ
                  if (qbases[i].maxAr >= 0)
                  {
                        // InstrumentationTimer timer("heap");

                    std::vector<int> new_;
                    new_ = a.mis_pos;
                    new_.push_back(a.aln_length+1);
                    // new_.push_back(a.num_gapo);


                    std::vector<char> newbase_;
                    newbase_ = a.mis_bases;
                    newbase_.push_back(c);

                    // cout << "PUSH MISMATCH = " << i << " ep-sp=" << qbases[i].U - qbases[i].L << endl;
                    heap.push(aln_entry_t(a.i, qbases[i].L, qbases[i].U, a.num_mm+1, (nR-qbases[i].maxAr), p, 0, j+1, a.aln_length+1, c, qbases[i].b, new_, newbase_ ));
                    // for(int r=0; r < new_.size(); r++)
                    //   cout << "-- IN CODE RETRIEVING " << new_[r] << " j = " << a.num_snps << endl;
                  }
                }


            }

        if (debug){
          // cout << repeat(s, a.aln_length+2) << "&&& END MISMATCH STYLE ALLOWED !!!" << endl;
          cout << "## " << "&&& END MISMATCH STYLE ALLOWED !!!" << endl;
        }
      }} else {

          // {
          //       InstrumentationTimer timer("idxBWT_exact");
                interval = idxBWT->exact_match( c, a.L, a.U);
          // }
          sp = interval.first;
          ep = interval.second;
          ulint maxArm;
          // {
          //       InstrumentationTimer timer("RMQ");
                // maxArm = RMQ.query(sp, ep-1);
                // cout << "3) RMQ(" << sp << ", " << ep-1 << ") traditional = " <<  RMQ.query(sp, ep-1) ;
                if (sp <= ep-1){
                  if (ep-1 >= nR){
                    if (sp >= nR){
                      // sp = vec.size()-1;
                      max_Ar_idx = RMQsdsl(nR-1, nR-1);
                      maxArm = idxBWT->convertToTextCoordinate(nR);

                    } else {
                      // ep = vec.size();
                      max_Ar_idx = RMQsdsl(sp, nR-1);
                      maxArm = idxBWT->convertToTextCoordinate(max_Ar_idx);
                    }
                    // cout << " and RMQ-sdsl = " << maxArm << endl;
                  } else {
                    max_Ar_idx = RMQsdsl(sp, ep-1);
                    maxArm = idxBWT->convertToTextCoordinate(max_Ar_idx);
                    // cout << " and RMQ-sdsl = " << maxArm << endl;
                  }

              } else {
                maxArm = 0;
                // cout << " and RMQ-sdsl = 0" << endl;
              }
              // assert(RMQ.query(sp, ep-1) == maxArm);


          // }
          // maxAr = RMQ.query(sp, ep-1);
          if (debug){
            // cout << repeat(s, a.aln_length+2) << "BEFORE Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << " max_Ar_ori=" << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << endl;
            // cout << repeat(s, a.aln_length+2) << "@@ TESTING FOR n-(i+1) = " <<  n-(a.i+allow_diff+1) << " BUT " << allow_diff << endl;

            cout << "## " << "BEFORE Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << " max_Ar_ori=" << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << endl;
            cout << "## " << "@@ TESTING FOR n-(i+1) = " <<  n-(a.i+allow_diff+1) << " BUT " << allow_diff << endl;
          }
          // LZ
          // if ( (sp < ep) && (maxArm >= n-(a.i+1)) ){
          // RLZ
          if ( (sp < ep) && (maxArm >= 0) ){

          // if ( (sp < ep) && (maxArm >= (a.num_gapo+1)) ){
            if (debug){
              // cout << repeat(s, a.aln_length+2) << " Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << ", allow_diff=" << a.state+1 << " iprime = " << a.num_gapo << " n-(a.i+1) = " << n << endl;
              cout << "## " << " Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << ", allow_diff=" << a.state+1 << " iprime = " << a.num_gapo << " n-(a.i+1) = " << n << endl;
            }
              // cout << "PUSH GOING EXACT !!! ep-sp=" << ep-sp<< endl;
            // {
            //       InstrumentationTimer timer("heap");
                  heap.push(aln_entry_t(a.i, sp, ep, a.num_mm, (nR-maxArm), p, a.state+1, j+1, a.aln_length+1, c, c , a.mis_pos,a.mis_bases));
            // }
            // allow_diff ++;
            i ++;
          }
        }


        if (debug){
          // cout << repeat(s, a.aln_length) << " Heap SIZE = " << heap.size() << endl;
            cout << "## " << " Heap SIZE = " << heap.size() << endl;
        }

				// perc = (100*j)/n;
				// if (perc > last_perc and (perc%5)) {
				// 	cout << " " << perc << "% done." << endl;
				// 	last_perc=perc;
				// }
        //
				// c = idxBWT->at( p );
				// p = idxBWT->LF( p );
        //
				// interval = idxBWT->exact_match( c, sp, ep );
				// sp = interval.first;
				// ep = interval.second;
        //
				// // cout << "#### sp = " << sp << ", ep = " << ep << endl;
				// max_Ar = RMQ.query(sp, ep - 1);
				// // cout << "Max_Ar = " << max_Ar << endl;
				// // cout << "i= " << i << " j= " << j << " c= " << c << " p= " << p << " sp=" << sp << " ep=" << ep << endl;
        // aln = a;
			}

      // j = j;
      if (debug){
        cout << "J final=" << j << endl;
      }

			if (a.num_gapo - 1 == -1) {
				// outfile << j << "\t" << c << "\t-\t-" <<endl;
        cout << j << "\t" << c << "\t-\t-" <<endl;
        if (int(c)==36){
          cout << "HERE $" << endl;
          break;
        }
			} else {
				// cout << "(" << iprime-(j-i) << "," << j-i-1 << ", " << c << ")" << endl;
        if (debug){
          cout << "c=" << a.b << " num_gapo (iprime) = " << a.num_gapo << ", num_snps = " << a.num_snps << " a.i = " << a.i << endl;
        }

        // if ( (a.num_gapo-(a.num_snps-a.i)) == 0) {
        //   cout << "BIRL" << endl;
        //   cout << "(" << (a.num_gapo-(a.num_snps-a.i)) << "," << a.aln_length + 1 << ")" <<  endl;
        // }
        // else {
          if (a.num_mm > 0){
                pos_factor = ((a.num_gapo)-(a.num_snps-a.i));
                which_chr = 0;
                while(pos_factor >= chrom_lens[which_chr]){
                  which_chr ++;
                }

                outfile << j << "\t" << ((a.num_gapo)-(a.num_snps-a.i)) << "\t" << a.aln_length + 1 << "\tM" << "\t" << Chref_rank(pos_factor+1) << endl;
                cout << j << "\t" << ((a.num_gapo)-(a.num_snps-a.i)) << "\t" << a.aln_length + 1 << "\tM" << "\t" << Chref_rank(pos_factor+1) << endl;
                for(int r=0; r < a.mis_pos.size(); r++){
                    cout << "-- RETRIEVING " << a.mis_pos[r] << " = " << a.mis_bases[r] << " j = " << a.mis_pos[r] << " track_len = " << track_len << endl;
                    // B[a.mis_pos[r]+((a.num_gapo)-(a.num_snps-a.i))] = 1;
                    cout << "ERROR HERE!!" << endl;
                    B[track_len+a.mis_pos[r]] = 1;
                    M.push_back(a.mis_bases[r]);
                }
                track_len += a.aln_length+1;
                I[track_len] = 1;

                // cout << j << "\t" << ((a.num_gapo)-(a.aln_length)) << "\t" << a.aln_length + 1 << "\tM" <<  endl;
          } else {
            pos_factor = ((a.num_gapo)-(a.num_snps-a.i));
            which_chr = 0;
            while(pos_factor >= chrom_lens[which_chr]){
              which_chr ++;
            }
            outfile << j<< "\t" << ((a.num_gapo)-(a.num_snps-a.i)) << "\t" << a.aln_length + 1 << "\tN" << "\t" << Chref_rank(pos_factor+1) << endl;
            cout << j<< "\t" << ((a.num_gapo)-(a.num_snps-a.i))<< "\t" << a.aln_length + 1 << "\tN" << "\t" << Chref_rank(pos_factor+1) << endl;
            track_len += a.aln_length+1;
            I[track_len] = 1;
            // cout << j<< "\t" << ((a.num_gapo)-(a.aln_length)) << "\t" << a.aln_length + 1 << "\tN" << endl;
          }
        // }
			}

      a.state = 0;
			i = j;

		}

    outfile.close();
    I[n-1] = 0;
    // cout << I << endl;
    // cout << B << endl;
    if (readablefiles){
      ofstream outfileM;
      outfileM.open(pref+"_outM.txt");
      outfileM << M;
      outfileM.close();
    }


    sd_vector<> sdbI(I);
    store_to_file(sdbI, pref+"_sdbI.sdsl");
    if (readablefiles){
      ofstream outfileI;
      outfileI.open(pref+"_outI.txt");
      outfileI << I;
      outfileI.close();
    }


    sd_vector<> sdbB(B);
    store_to_file(sdbB, pref+"_sdbB.sdsl");
    if (readablefiles){
      ofstream outfileB;
      outfileB.open(pref+"_outB.txt");
      outfileB << B;
      outfileB.close();
    }
	}

  Instrumentor::Get().EndSession();
	printRSSstat();
	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();

	if(total>=3600){

		uint h = total/3600;
		uint m = (total%3600)/60;
		uint s = (total%3600)%60;

		cout << "Total time: " << h << "h " << m << "m " << s << "s" << endl;

	}else if (total>=60){

		uint m = total/60;
		uint s = total%60;

		cout << "Total time: " << m << "m " << s << "s" << endl;

	}else{

		cout << "Total time: " << total << "s" << endl;

	}
}
