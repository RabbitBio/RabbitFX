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


void consumer_se_fastq_task(FXReader<FQ_SE> &m_reader, Counter *counter) {
  long line_sum = 0;
  rabbit::int64 id = 0;
  while (true) {
    auto data = m_reader.get_formated_reads();
    if(data.size() == 0) break;
		for(Reference &read : data){
			for(int i = 0; i < read.seq.size(); i++){
				switch(read.seq[i]){
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
  }
}

int main(int argc, char **argv) {
  std::string file1 = std::string(argv[1]); //"/home/user_home/haoz/data/fastv_experiment_data/SRR1030141_1.fastq";
  int th = std::stoi(argv[2]);  // thread number
  
  FXReader<FQ_SE> m_reader(file1);
	Counter* counters[th];
  std::thread** threads = new thread*[th];
  for (int t = 0; t < th; t++) {
		counters[t] = new Counter{0,0,0,0};
    threads[t] = new std::thread(std::bind(consumer_se_fastq_task, std::ref(m_reader), counters[t]));
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
