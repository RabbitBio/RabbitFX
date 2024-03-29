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
#include "io/RabbitFX.h"

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


void consumer_pe_fastq_task(FXReader<FQ_PE> &m_reader, Counter *counter) {
  long line_sum = 0;
  rabbit::int64 id = 0;
  rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk;
  while (true) {
    auto data = m_reader.get_formated_reads_nocp();
    if(data.first.chunk_p == NULL || data.second.chunk_p == NULL) break;
		for(neoReference &read : data.first.vec){
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
		for(neoReference &read : data.second.vec){
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
		m_reader.release_chunk(data);
  }
}

int main(int argc, char **argv) {
  std::string file1 = std::string(argv[1]); //"/home/user_home/haoz/data/fastv_experiment_data/SRR1030141_1.fastq";
  //std::string file2 = "/home/old_home/haoz/ncbi/public/sra/mashscreen_test/ERR1711677_2.fastq";
  std::string file2 = std::string(argv[2]); // "/home/user_home/haoz/data/fastv_experiment_data/SRR1030141_2.fastq";
  int th = std::stoi(argv[3]);  // thread number
  
  FXReader<FQ_PE> m_reader(file1, file2);
	Counter* counters[th];
  std::thread** threads = new thread*[th];
  for (int t = 0; t < th; t++) {
		counters[t] = new Counter{0,0,0,0};
    threads[t] = new std::thread(std::bind(consumer_pe_fastq_task, std::ref(m_reader), counters[t]));
  }
  m_reader.join_producer();
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
