#include <string>
#include <iostream>
#include "thirdparty/CLI11.hpp"
#include <thread>
#include "io/Formater.h"
#include <sys/time.h>
#include <cstdint>
#include <vector>
#include "thirdparty/robin_hood.h"
#include <unordered_map>
#include "io/RabbitFX.h"

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


void consumer_fasta_task(FXReader<FA> &m_reader, Counter *counter){
  std::vector<Reference> data;
  int ref_num;
  std::unordered_map<uint64_t, uint64_t> &gid2len = counter->gid2len;
  std::unordered_map<uint64_t, std::string> &gid2name = counter->gid2name;

  while ( true ) {
		auto data = m_reader.get_formated_reads();
		if(data.size() == 0) break;
		//rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
		for(Reference &r: data){
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
  //FA faReader(filename);
  FXReader<FA> m_reader(filename);
  std::thread **threads = new std::thread *[th];
	Counter* counters[th];
  for (int t = 0; t < th; t++) {
		counters[t] = new Counter{0,0,0,0};
    threads[t] = new std::thread(std::bind(consumer_fasta_task, std::ref(m_reader), counters[t]));
  }
  m_reader.join_producer();
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
  return 0;
}
