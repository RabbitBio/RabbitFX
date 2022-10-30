#include "./FastxStream.h"
#include "./FastxChunk.h"
#include "assert.h"

class FQ{
typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FqChunkQueue;  
typedef rabbit::fq::FastqDataPool FqDataPool;
private:
  rabbit::fq::FastqFileReader* fqFileReader;
  FqDataPool* dp_;
  FqChunkQueue* dq_;
public:
  FQ(std::string file1, std::string file2){
    dq_ = new  FqChunkQueue(128, 1); // data queue
    dp_ = new FqDataPool(256, 1 << 22); // data pool
    fqFileReader = new rabbit::fq::FastqFileReader(file1, *dp_, file2, false);
    //start_producer();
  }
  void start_producer(){
    std::thread producer([&] {
      int n_chunks = 0;
      while (true) {
        rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk;
        fqchunk->chunk = fqFileReader->readNextPairChunk1();
        if (fqchunk->chunk == NULL) break;
        n_chunks++;
        dq_->Push(n_chunks, fqchunk->chunk);
      }

      dq_->SetCompleted();
    }); 
  }
  template <typename REFTYPE>
  int get_formated_reads(vector<REFTYPE> &data) {
    rabbit::fq::FastqChunk *fqchunk;  // = new rabbit::fa::FastaChunk;
    rabbit::int64 id = 0;
    int ref_num;
    if (dq_->Pop(id, fachunk)) {
      // rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
      ref_num = rabbit::fa::chunkFormat(*fachunk, data);
      //------------------relaease-----------------//
      rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
      do {
        dp_->Release(tmp);
        tmp = tmp->next;
      } while (tmp != NULL);
    }else{
      ref_num = 0;
    }
    return ref_num;
  }

  ~FQ(){
    delete fqFileReader;
    delete dp_;
    delete dq_;
  }
};

class FA{
typedef rabbit::core::TDataQueue<rabbit::fa::FastaChunk> FaChunkQueue;
typedef rabbit::fa::FastaDataPool FaDataPool;
public:
  rabbit::fa::FastaFileReader *faFileReader;
  FaDataPool* dp_;
  FaChunkQueue* dq_;
  std::thread* producer_;
  std::thread** consumers_;
  int consumer_th_;
public:
  FA(std::string &file){
    printf("in constuctor of FA!, fname: %s\n", file.c_str());
    dq_ = new  FaChunkQueue(128, 1); // data queue
    dp_ = new FaDataPool(256, 1 << 24); // data pool
    faFileReader = new rabbit::fa::FastaFileReader(file, dp_, false);
  }

  void start_producer() {
    producer_ = new std::thread([&]{
      int n_chunks = 0;
      while (true) {
        rabbit::fa::FastaChunk *fachunk;// = new rabbit::fa::FastaChunk;
        // fachunk = faFileReader->readNextChunkList();
        fachunk = faFileReader->readNextChunk();
        if (fachunk == NULL) break;
        n_chunks++;
        dq_->Push(n_chunks, fachunk);
      }
      dq_->SetCompleted();
    });
  }
  template<typename T>
	void process_data_mt(const int tn, T (*func_ptr)(Reference &), vector<vector<T> > &v_res_data){ //tn means thread number
		consumers_ = new std::thread*[tn];
    assert(v_res_data.size() == tn);
    this->consumer_th_ = tn;
    for (int i = 0; i < tn; i++) {
      rabbit::int64 id = 0;
      vector<T> &res_data = v_res_data[i];

      consumers_[i] = new thread([&]() {
        cout << "starting worker: " << i << endl;
        rabbit::fa::FastaChunk *fachunk;  // = new rabbit::fa::FastaChunk;
        while (dq_->Pop(id, fachunk)) {
          // rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
          std::vector<Reference> data;
          int ref_num = rabbit::fa::chunkFormat(*fachunk, data);
          for (Reference &r : data) {
            res_data.emplace_back(func_ptr(r));
          }
          //------------------relaease-----------------//
          rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
          do {
            dp_->Release(tmp);
            tmp = tmp->next;
          } while (tmp != NULL);
        }
      });
    }
  }

  template<typename REFTYPE>
  int get_formated_reads(vector<REFTYPE>& data){
    rabbit::fa::FastaChunk *fachunk;  // = new rabbit::fa::FastaChunk;
    rabbit::int64 id = 0;
    int ref_num;
    if (dq_->Pop(id, fachunk)) {
      // rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
      ref_num = rabbit::fa::chunkFormat(*fachunk, data);
      //------------------relaease-----------------//
      rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
      do {
        dp_->Release(tmp);
        tmp = tmp->next;
      } while (tmp != NULL);
    }else{
      ref_num = 0;
    }
    return ref_num;
  }

  ~FA(){
    delete faFileReader;
    delete dq_;
    delete dp_;
    delete producer_;
  }
};


// Just a wapper
template<class T>
class FXReader{
public:
  T reader_;
  FXReader(std::string fn)
    : reader_(fn)
  {
    reader_.start_producer();
  }
  FXReader(const std::string &fn1, const std::string &fn2)
    : reader_(fn1, fn2)
  {
    reader_.start_producer();
  }
 // bool pop_dq(rabbit::int64 &id, ){
 //   return reader_.dq_.Pop(id, chunk);
 // }
  template<typename REFTYPE>
  int get_formated_reads(vector<REFTYPE> &data){
    return reader_.get_formated_reads(data);
  }
  void join_producer(){
    reader_.producer_->join();
  }
  void join_consumers(){
    for(int i = 0; i < reader_.consumer_th_; i++){
      reader_.consumers_[i]->join();
    }
  }
};
