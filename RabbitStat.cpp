#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include <string>
#include <iostream>
#include "thirdparty/CLI11.hpp"
#include "io/DataQueue.h"
#include <thread>
#include "io/Formater.h"
#include <sys/time.h>
#include <cstdint>
#include <vector>
#include "thirdparty/robin_hood.h"
#include <unordered_map>

//typedef robin_hood::unordered_map unordered_map;

double get_time(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

typedef int pfunc(int, char**);
struct Counter{
	uint64_t A, T, G, C;
	std::unordered_map<uint64_t, uint64_t> gid2len;
	std::unordered_map<uint64_t, std::string> gid2name;
};

typedef rabbit::core::TDataQueue<rabbit::fa::FastaChunk> FaChunkQueue;

int producer_fasta_task(std::string file, rabbit::fa::FastaDataPool* fastaPool, FaChunkQueue &dq) {
  //rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1 << 22);
  rabbit::fa::FastaFileReader *faFileReader;
  faFileReader = new rabbit::fa::FastaFileReader(file, *fastaPool, false);
  int n_chunks = 0;
  while (true) {
    rabbit::fa::FastaChunk *fachunk = new rabbit::fa::FastaChunk;
    //fachunk = faFileReader->readNextChunkList();
    fachunk = faFileReader->readNextChunk();
    if (fachunk == NULL) break;
    n_chunks++;
		dq.Push(n_chunks, fachunk);
  }
	dq.SetCompleted();
  return 0;
}

//typedef core::TDataQueue<FastaDataChunk> FastaDataQueue;
void consumer_fasta_task(rabbit::fa::FastaDataPool *fastaPool,  FaChunkQueue &dq, Counter *counter) {
  long line_sum = 0;
  rabbit::int64 id = 0;
  //rabbit::fa::FastaDataChunk *fachunk = new rabbit::fa::FastaDataChunk;
	std::vector<uint64_t> lengths;
	std::unordered_map<uint64_t, uint64_t> &gid2len = counter->gid2len;
	std::unordered_map<uint64_t, std::string> &gid2name = counter->gid2name;

  rabbit::fa::FastaChunk *fachunk;// = new rabbit::fa::FastaChunk;
	while (dq.Pop(id, fachunk)) {
		//rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
    std::vector<Reference> data;
    int ref_num = rabbit::fa::chunkFormat(*fachunk, data);
		for(Reference &r: data){
			//count length
			//std::cout << "ref num: " << ref_num << " reference global id: " << r.gid << " name: " << r.name
			//					<< " length: " << r.length << " seq size: " << r.seq.length() << std::endl;
			gid2len[r.gid] += r.length;
			gid2name[r.gid] += r.name;
      //count ATGC
			const char* base = r.seq.c_str();
			for(int i = 0; i < r.length; i++){
				switch(base[i]){
				case 'A':
					counter->A++; break;
				case 'T':
					counter->T++; break;
				case 'G':
					counter->G++; break;
				case 'C':
					counter->C++; break;
				default:
					break;
				}
			}
		}
    //-----relaease
		rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
    do {
      fastaPool->Release(tmp);
      tmp = tmp->next;
    } while (tmp != NULL);
  }
}

int main(int argc, char **argv) {
  CLI::App app{"Wellcome to RabbitIO"};
  CLI::Option* opt;
  //std::string file1 = "/home/old_home/haoz/workspace/data/hg38/hg38.fa";
  std::string filename;
  int th;  // thread number
  app.add_option("file", filename, "input file name")->required();
		//app.add_option("-f, --file", filename, "input file name")
		//	->required();
  app.add_option("-t, --threads", th, "worktreads")
    ->default_val(4);
  //----------------------------------------
  CLI11_PARSE(app, argc, argv);
  rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1 << 24);
  FaChunkQueue queue1(128, 1);
  std::thread producer(producer_fasta_task, filename, fastaPool, std::ref(queue1));
  std::thread **threads = new std::thread *[th];
	Counter* counters[th];
  for (int t = 0; t < th; t++) {
		counters[t] = new Counter{0,0,0,0};
    threads[t] = new std::thread(std::bind(consumer_fasta_task, fastaPool, std::ref(queue1), counters[t]));
  }
  producer.join();
  for (int t = 0; t < th; t++) {
    threads[t]->join();
  }
	uint64_t ca = 0, ct = 0, cc = 0, cg = 0;
	std::unordered_map<uint64_t, uint64_t> global_g2l;
	std::unordered_map<uint64_t, std::string> global_g2n;
	uint64_t sum_length = 0;
	for(int t = 0; t < th; t++){
		ca += counters[t]->A;
		ct += counters[t]->T;
		cc += counters[t]->C;
		cg += counters[t]->G;
		for(auto& x : counters[t]->gid2len){
			global_g2l[x.first] += x.second;
		}
		for(auto& x : counters[t]->gid2name){
			global_g2n[x.first] += x.second;
		}
	}
	for(int i = 0; i < global_g2l.size(); i++){
		std::cout << global_g2n[i]<< " - " << global_g2l[i] << std::endl;
	}
	//std::cout << "total length: " << sum_length << std::endl;
	//std::cout << "ATCG infromatic: " << ca << " " << ct << " " << cc << " " << cg << std::endl;
	printf("A:%ld, T:%ld, C:%ld, G:%ld\n", ca, ct, cc, cg);
  delete fastaPool;
  for (int t = 0; t < th; t++) {
    delete threads[t];
  }
  return 0;
}
