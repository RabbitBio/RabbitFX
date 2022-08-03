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

double get_time(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FqChunkQueue;
typedef int pfunc(int, char**);
struct Counter{
	uint64_t A, T, G, C;
};


int producer_pe_fastq_task(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool, FqChunkQueue &dq) {
  rabbit::fq::FastqFileReader *fqFileReader;
  fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, file2, false);
  int n_chunks = 0;
  int line_sum = 0;
  while (true) {
    rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk;
    fqchunk->chunk = fqFileReader->readNextPairChunk1();
    if (fqchunk->chunk == NULL) break;
    n_chunks++;
    //std::cout << "readed chunk: " << n_chunks << std::endl;
    dq.Push(n_chunks, fqchunk->chunk);
  }

  dq.SetCompleted();
  delete fqFileReader;
  std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
  return 0;
}

void consumer_pe_fastq_task(rabbit::fq::FastqDataPool *fastqPool, FqChunkQueue &dq, Counter *counter) {
  long line_sum = 0;
  rabbit::int64 id = 0;
  rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk;
  while (dq.Pop(id, fqchunk->chunk)) {
		std::vector<neoReference> data1;
		data1.resize(10000);
		rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk*)(fqchunk->chunk->left_part), data1);
		for(neoReference &read : data1){
			for(int i = 0; i < read.lseq; i++){
				switch(read.base[read.pseq + i]){
				case 'A':
					counter->A++; break;
				case 'T':
					counter->T++; break;
				case 'G':
					counter->G++; break;
				case 'C':
					counter->C++; break;
				}
			}
		}
		std::vector<neoReference> data2;
		data2.resize(10000);
		rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk*)(fqchunk->chunk->right_part), data2);
		for(neoReference &read : data2){
			for(int i = 0; i < read.lseq; i++){
				switch(read.base[read.pseq + i]){
				case 'A':
					counter->A++; break;
				case 'T':
					counter->T++; break;
				case 'G':
					counter->G++; break;
				case 'C':
					counter->C++; break;
				}
			}
		}
    fastqPool->Release(fqchunk->chunk->left_part);
    fastqPool->Release(fqchunk->chunk->right_part);
  }
}

int main(int argc, char **argv) {
  std::string file1 = std::string(argv[1]); //"/home/user_home/haoz/data/fastv_experiment_data/SRR1030141_1.fastq";
  //std::string file2 = "/home/old_home/haoz/ncbi/public/sra/mashscreen_test/ERR1711677_2.fastq";
  std::string file2 = std::string(argv[2]); // "/home/user_home/haoz/data/fastv_experiment_data/SRR1030141_2.fastq";
  int th = std::stoi(argv[3]);  // thread number
  rabbit::fq::FastqDataPool *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
  FqChunkQueue queue1(128, 1);
  std::thread producer(producer_pe_fastq_task, file1, file2, fastqPool, std::ref(queue1));
  std::thread **threads = new std::thread *[th];
	Counter* counters[th];
  for (int t = 0; t < th; t++) {
		counters[t] = new Counter{0,0,0,0};
    threads[t] = new std::thread(std::bind(consumer_pe_fastq_task, fastqPool, std::ref(queue1), counters[t]));
  }
  producer.join();
  for (int t = 0; t < th; t++) {
    threads[t]->join();
  }
	uint64_t ca = 0, ct = 0, cc = 0, cg = 0;
	for(int t = 0; t < th; t++){
		ca += counters[t]->A;
		ct += counters[t]->T;
		cc += counters[t]->C;
		cg += counters[t]->G;
	}
	std::cout << ca << " " << ct << " " << cc << " " << cg << std::endl;
 
  delete fastqPool;
  for (int t = 0; t < th; t++) {
    delete threads[t];
  }
  return 0;
}
// pcmode(){
//   rabbit::fq::FastqDataPool *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
//   FqChunkQueue queue1(128, 1);
//   std::thread producer(producer_pe_fastq_task, file1, file2, fastqPool, std::ref(queue1));
//   std::thread **threads = new std::thread *[th];
// 	read_chunks();
// 	Counter* counters[th];
//   for (int t = 0; t < th; t++) {
// 		counters[t] = new Counter{0,0,0,0};
//     threads[t] = new std::thread(std::bind(consumer_pe_fastq_task, fastqPool, std::ref(queue1), counters[t]));
//   }
//   producer.join();
//   for (int t = 0; t < th; t++) {
//     threads[t]->join();
//   }
// }
// int main(){
// 	readfile();
// 	process
// }
