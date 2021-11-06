#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include <string>
#include <iostream>
#include "thirdparty/CLI11.hpp"
#include "io/DataQueue.h"
#include <thread>
#include "io/Formater.h"
#include <sys/time.h>

double get_time(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FqChunkQueue;
typedef int pfunc(int, char**);

int count_line(rabbit::fq::FastqChunk *fqchunk) { return 1000; }

int producer_pe_fastq_task(std::string file, std::string file2, rabbit::fq::FastqDataPool &fastqPool, FqChunkQueue &dq) {
  rabbit::fq::FastqFileReader *fqFileReader;
  fqFileReader = new rabbit::fq::FastqFileReader(file, fastqPool, file2, false);
  int n_chunks = 0;
  int line_sum = 0;
  while (true) {
    rabbit::fq::FastqDataPairChunk *fqdatachunk = new rabbit::fq::FastqDataPairChunk;
    fqdatachunk = fqFileReader->readNextPairChunk();
    if (fqdatachunk == NULL) break;
    n_chunks++;
    //std::cout << "readed chunk: " << n_chunks << std::endl;
    dq.Push(n_chunks, fqdatachunk);
  }

  dq.SetCompleted();
  delete fqFileReader;
  std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
  return 0;
}

void consumer_pe_fastq_task(rabbit::fq::FastqDataPool &fastqPool, FqChunkQueue &dq) {
  long line_sum = 0;
  rabbit::int64 id = 0;
  std::vector<neoReference> data;
  rabbit::fq::FastqDataPairChunk *fqdatachunk = new rabbit::fq::FastqDataPairChunk;
  data.resize(10000);
  while (dq.Pop(id, fqdatachunk)) {
		line_sum += rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk*)(fqdatachunk->left_part), data, true);
		line_sum += rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk*)(fqdatachunk->right_part), data, true);
    fastqPool.Release(fqdatachunk->left_part);
    fastqPool.Release(fqdatachunk->right_part);
  }
}

int producer_fastq_task(std::string file, rabbit::fq::FastqDataPool& fastqPool, rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq){
  rabbit::fq::FastqFileReader fqFileReader(file, fastqPool);
  rabbit::int64 n_chunks = 0; 
  while(true){ 
		rabbit::fq::FastqDataChunk* fqdatachunk;// = new rabbit::fq::FastqDataChunk;
    fqdatachunk = fqFileReader.readNextChunk(); 
    if (fqdatachunk == NULL) break;
    n_chunks++;
    //std::cout << "readed chunk: " << n_chunks << std::endl;
    dq.Push(n_chunks, fqdatachunk);
  }
  dq.SetCompleted();
  std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
  return 0;
}

void consumer_fastq_task(rabbit::fq::FastqDataPool& fastqPool, rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq){
    long line_sum = 0;
    rabbit::int64 id = 0;
    std::vector<neoReference> data;
		rabbit::fq::FastqDataChunk* fqdatachunk;// = new rabbit::fq::FastqDataChunk;
    data.resize(10000);
    while(dq.Pop(id, fqdatachunk)){
        line_sum += rabbit::fq::chunkFormat(fqdatachunk, data, true);
        fastqPool.Release(fqdatachunk);
    }
    std::cout << "line_sum: " << line_sum << std::endl;
}

void print_chunk(rabbit::fa::FastaDataChunk *chunk) {
  std::cout << "chunk size: " << chunk->size << std::endl;
  std::cout << "chunk head: " << std::string((char *)chunk->data.Pointer(), 100) << std::endl;
}

void print_fachunkpart_info(rabbit::fa::FastaChunk *fachunk) {
  std::cout << "------------------chunk info-----------------" << std::endl;
  print_chunk(fachunk->chunk);
  rabbit::fa::FastaDataChunk *current_chunk = fachunk->chunk;
  while (current_chunk->next != NULL) {
    std::cout << "next" << std::endl;
    current_chunk = current_chunk->next;
    print_chunk(current_chunk);
  }
}

int producer_fasta_task(std::string file) {
  rabbit::fa::FastaDataPool fastaPool(256, 1 << 22);
  rabbit::fa::FastaFileReader faFileReader(file, fastaPool, false);
  int n_chunks = 0;
  int line_sum = 0;
  while (true) {
    rabbit::fa::FastaChunk *fachunk = new rabbit::fa::FastaChunk;
    fachunk = faFileReader.readNextChunkList();
    //fachunk = faFileReader->readNextChunk();
    if (fachunk == NULL) break;
    n_chunks++;
		//std::cout << n_chunks << std::endl;
    // line_sum += count_line(fqchunk);
    //std::vector<Reference> data;
    //print_fachunkpart_info(fachunk);
    //-----relaease
    rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
    do {
      fastaPool.Release(tmp);
      tmp = tmp->next;
    } while (tmp != NULL);
    // line_sum += rabbit::fa::chunkFormat(*fachunk, data);
  }
  std::cout << "file " << file << " has " << line_sum << " lines" << std::endl;
  return 0;
}

int test_fastq_pe(int argc, char **argv) {
  std::string file1 = "/home/old_home/haoz/ncbi/public/sra/fastv_test/SRX7918726.sra_1.fastq";
  std::string file2 = "/home/old_home/haoz/ncbi/public/sra/fastv_test/SRX7918726.sra_2.fastq";
  int th = 20;  // thread number
  rabbit::fq::FastqDataPool fastqPool(256, 1 << 22);
  FqChunkQueue queue1(128, 1);
  std::thread producer(producer_pe_fastq_task, file1, file2, std::ref(fastqPool), std::ref(queue1));
  std::thread **threads = new std::thread *[th];
  for (int t = 0; t < th; t++) {
    threads[t] = new std::thread(std::bind(consumer_pe_fastq_task, std::ref(fastqPool), std::ref(queue1)));
  }
  producer.join();
  for (int t = 0; t < th; t++) {
    threads[t]->join();
  }
  //delete fastqPool;
  for (int t = 0; t < th; t++) {
    delete threads[t];
  }
  return 0;
}

// test SE
int test_fastq_se(int argc, char** argv){
  std::string file = "/home/old_home/haoz/workspace/QC/fastp_dsrc/out_1.fq";
  //std::string file = "/home/old_home/haoz/workspace/data/hg19/hg19.fa";
  //---------------cmd parser----------------
  CLI::App app{"Wellcome to RabbitIO"};
  CLI::Option* opt;
  std::string filename ;
  int th;
  app.add_option("-f, --file", filename, "input file name")
    ->default_val(file);
  app.add_option("-t, --threads", th, "worktreads")
    ->default_val(2);
  //----------------------------------------
  CLI11_PARSE(app, argc, argv);
  if(app.count("-f"))
    std::cout << "filename: " << filename << std::endl;
  else{
    std::cout << "-f not find, use default: " << filename << std::endl;
  }
  rabbit::fq::FastqDataPool fastqPool(32, 1<<22);
  rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);
  std::thread producer(producer_fastq_task, filename, std::ref(fastqPool), std::ref(queue1));
  vector<std::thread> threads;
  for(int t = 0; t < th; t++){
    threads.emplace_back(thread(consumer_fastq_task, std::ref(fastqPool), std::ref(queue1)));
  }
  producer.join();
  for(int t = 0; t < th; t++){
    threads[t].join();
  }
  return 0;
}

int test_fasta(int argc, char** argv){
	producer_fasta_task("/home/old_home/haoz/workspace/data/hg19/hg19.fa");
	return 0;
}

inline void check_sucess(pfunc func, int argc, char** argv, std::string desc){
  std::cout << "runing " << desc << std::endl;
	double start_time = get_time();
  if(func(argc, argv) == 0){
    std::cout << "["<< desc << "] runing sucess!";
  }else{
    std::cout << "["<< desc << "] not sucess!";
  }
	double end_time = get_time();
	std::cout << " time: " << end_time - start_time << " s" << std::endl;
}

int main(int argc, char** argv){
  check_sucess(test_fastq_pe, argc, argv, "test fastq pair end");
  check_sucess(test_fastq_se, argc, argv, "test fastq single end");
  check_sucess(test_fasta,    argc, argv, "test fasta process");
}

/*
FastqReader processer( );
io::data::chunk<FastqChunk> fqchunk = processer.get_chunk(file);
//or
io::data::seq<FastqSeq>  fqseqs = processer.get_formated_seq(file);

processer.process_seq(worker_num, func, param);

thread_pool.process(producer_task, file);
thread_pool.process(consumer_task, file);

return 0;
}*/
