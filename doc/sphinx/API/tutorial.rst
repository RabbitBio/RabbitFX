.. figure:: rabbitio.png
   :alt: RabbitIO

RabbitIO Tutorial
=================

Installation
------------

Dependency
~~~~~~~~~~

1. c++11
2. `zlib <https://zlib.net/>`__ 
   
Using Cmake (recommend) 
~~~~~~~~~~~~~~~~~~~~~~~

Copy folder `io` to your program, and then you can integrate RabbitIO in your ``CMakeLists.txt``:

.. code:: cmake

  AUX_SOURCE_DIRECTORY(. SOURCE_LIST)
  ADD_LIBRARY(io_lib ${SOURCE_LIST})
  TARGET_LINK_LIBRAIES(io_lib z)

Highlight
---------

-  RabbitIO highly support multi-core paltform
-  RabbitIO efficiency processing FASTQ/FASTA files
-  Concurrency data pool and data queue.
-  Non-copy read format
-  Bit-based sequencing processing

Illustration
------------

.. figure:: pipeline.png
   :alt: Pipeline

Case study
----------

-  `RabbitIO-Ktrim <https://github.com/RabbitBio/RabbitIO-Casestudy/tree/master/RabbitIO-Ktrim>`__
-  `RabbitIO-fastp <https://github.com/RabbitBio/RabbitIO-Casestudy/tree/master/RabbitIO-fastp>`__
-  `RabbitIO-Mash <https://github.com/RabbitBio/RabbitIO-Casestudy/tree/master/RabbitIO-Mash>`__

Runing Example in main.cpp and TestCount.cpp
--------------------------------------------

.. code:: bash

    cd RabbitIO
    mkdir build && cd build
    cmake ..
    make
    #then there is an test file in build file
    time ./test
    time ./testcount

**Note:** We integtated `CLI <https://github.com/CLIUtils/CLI11>`__ as
default command line patser for user convenience

FASTQ data example
------------------

Single-end data processing example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Example of define a mult-threading task

.. code-block:: cpp

  int test_fastq_se(int argc, char** argv){
    std::string file = "/home/old_home/haoz/workspace/QC/out_1.fq";
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
    rabbit::fq::FastqDataPool *fastqPool = new rabbit::fq::FastqDataPool(32, 1<<22);
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);
    std::thread producer(producer_fastq_task, filename, fastqPool, std::ref(queue1));
    std::thread** threads = new std::thread*[th];
    for(int t = 0; t < th; t++){
      threads[t] = new std::thread(std::bind(consumer_fastq_task, fastqPool, std::ref(queue1)));
    }
    producer.join();
    for(int t = 0; t < th; t++){
      threads[t]->join();
    }
    //-----free
    delete fastqPool;
    delete[] threads;
    return 0;
  }

2. example of define producer and consumer task 

.. code-block:: cpp

  int producer_fastq_task(std::string file, rabbit::fq::FastqDataPool* fastqPool, rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq){
    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool);
    rabbit::int64 n_chunks = 0;
    while(true){
      rabbit::fq::FastqDataChunk* fqdatachunk;// = new rabbit::fq::FastqDataChunk;
      fqdatachunk = fqFileReader->readNextChunk();
      if (fqdatachunk == NULL) break;
      n_chunks++;
      //std::cout << "readed chunk: " << n_chunks << std::endl;
      dq.Push(n_chunks, fqdatachunk);
    }
    dq.SetCompleted();
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    return 0;
  }
  
  void consumer_fastq_task(rabbit::fq::FastqDataPool* fastqPool, rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq){
      long line_sum = 0;
      rabbit::int64 id = 0;
      std::vector<neoReference> data;
  	rabbit::fq::FastqDataChunk* fqdatachunk;// = new rabbit::fq::FastqDataChunk;
      data.resize(10000);
      while(dq.Pop(id, fqdatachunk)){
        line_sum += rabbit::fq::chunkFormat(fqdatachunk, data, true);
        fastqPool->Release(fqdatachunk);
      }
      std::cout << "line_sum: " << line_sum << std::endl;
  }
  
Pair-end data processing example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example of processing Pair-end sequencing data is showed in file
`TestCount.cpp <https://github.com/RabbitBio/RabbitIO/blob/master/TestCount.cpp>`__. It is tested that compared to
`FQReader <https://github.com/rob-p/FQFeeder>`__, in the task of
counting ATCG of pair-end data, RabbitIO is 2 times faster in 20 thread.

RabbitIO is about 2G/s I/O speed now

FASTA data example
------------------

this is an example of reading and processing FASTA file

-  example code of using only one thread (count chunk number of input
   file): 
   
.. code-block:: cpp

  int proces_fasta_task(std::string file) {
    rabbit::fa::FastaDataPool *fastaPool = new rabbit::fa::FastaDataPool(256, 1 << 22);
    rabbit::fa::FastaFileReader *faFileReader;
    faFileReader = new rabbit::fa::FastaFileReader(file, *fastaPool, false);
    int n_chunks = 0;
    int line_sum = 0;
    while (true) {
      rabbit::fa::FastaChunk *fachunk = new rabbit::fa::FastaChunk;
      fachunk = faFileReader->readNextChunkList();
      if (fachunk == NULL) break;
      n_chunks++;
      //-----relaease
      rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
      do {
        fastaPool->Release(tmp);
        tmp = tmp->next;
      } while (tmp != NULL);
      // line_sum += rabbit::fa::chunkFormat(*fachunk, data);
    }
    std::cout << "file " << file << " has " << line_sum << " lines" << std::endl;
    return 0;
  }
  
  int test_fasta(int argc, char** argv){
    producer_fasta_task("/home/old_home/haoz/workspace/data/hg19/hg19.fa");
    return 0;
  }

Cite
----

RabbitIO paper is under review now.
