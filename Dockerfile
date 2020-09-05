FROM ubuntu:16.04


#
# Install general tools
#
RUN apt-get update && apt-get install -y \
    git \
    python \
    make \
    curl \
    bzip2 \
    libbz2-dev \
    libz-dev \
    liblzma-dev \
    g++ \
    wget \
    libncurses-dev \
    libcurl4-openssl-dev \
    unzip \
    ftp \
    vim



#-#-#-#-#-#-#-#-#-#-#-#-#
# BAM processing tools
#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/htslib-1.10.2 && make && make install && \
    rm -rf /tmp/htslib-1.10.2
RUN curl -kL https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/samtools-1.10 && make && make install && \
    rm -rf /tmp/samtools-1.10
RUN curl -kL https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/bcftools-1.10.2 && make && make install && \
    rm -rf /tmp/bcftools-1.10.2



#-#-#-#-#-#-#-#-#-#-#-#-#
# DNA analysis tools
#-#-#-#-#-#-#-#-#-#-#-#-#
#bwa
RUN curl -kL https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/bwa-0.7.17 && make && find /tmp/bwa-0.7.17/ -type f -executable -exec mv '{}' /usr/local/bin/ ';' && \
    rm -rf /tmp/bwa-0.7.17/
#minimap2
RUN cd /tmp && git clone https://github.com/lh3/minimap2 && \
    cd /tmp/minimap2 && make && \
    find /tmp/minimap2/ -type f -executable -exec mv '{}' /usr/local/bin/ ';'




#-#-#-#-#-#-#-#-#-#-#-#-#
# RNA analysis tools 
#-#-#-#-#-#-#-#-#-#-#-#-#
#STAR
RUN curl -kL https://github.com/alexdobin/STAR/archive/2.7.5c.tar.gz | tar -C /tmp -zxf - && \ 
    mv /tmp/STAR-2.7.5c/bin/Linux_x86_64_static/* /usr/local/bin/ && \
    rm -rf /tmp/STAR-2.7.5c/

#bowtie2 + tophat + cufflinks
#RUN curl -kL http://netix.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip -o /tmp/bowtie2-2.2.9-linux-x86_64.zip && \
#    unzip /tmp/bowtie2-2.2.9-linux-x86_64.zip -d /tmp && \
#    find /tmp/bowtie2-2.2.9/ -maxdepth 1 -type f -executable -exec mv '{}' /usr/local/bin/ ';' && \
#    rm -rf /tmp/bowtie2-2.2.9-linux-x86_64.zip /tmp/bowtie2-2.2.9
#RUN curl -kL https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz | tar -C /tmp -zxf - && \
#    find /tmp/tophat-2.1.1.Linux_x86_64/ -maxdepth 1 -executable -type f -exec mv '{}' /usr/local/bin/ ';' && \
#    rm -rf /tmp/tophat-2.1.1.Linux_x86_64
#RUN curl -kL http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar -C /tmp -zxf - && \
#    find /tmp/cufflinks-2.2.1.Linux_x86_64/ -maxdepth 1 -executable -type f -exec mv '{}' /usr/local/bin/ ';' && \
#    rm -rf /tmp/cufflinks-2.2.1.Linux_x86_64    


#-#-#-#-#-#-#-#-#-#-#-#-#
# de novo assembly tools
#-#-#-#-#-#-#-#-#-#-#-#-#
RUN apt-get update && apt-get install -y \
    libsparsehash-dev \
    bamtools \
    libbamtools-dev \
    libz-dev \
    autoconf
RUN curl -kL https://github.com/jts/sga/archive/v0.10.15.tar.gz | tar -C /tmp -zxf - && \
    mkdir /usr/lib/bamtools && ln -s /usr/lib/x86_64-linux-gnu/libbamtools.* /usr/lib/bamtools/ && \
    cd /tmp/sga-0.10.15/src && ./autogen.sh && ./configure --with-bamtools=/usr && make && make install && \
    rm -rf /tmp/sga-0.10.15




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Additional tools (using JAVA)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
RUN apt-get update && apt-get install -y picard-tools



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Other tools to install
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# SRA toolkit
RUN curl -kL http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.7.0/sratoolkit.2.7.0-ubuntu64.tar.gz | tar -C /tmp -zxf - 




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Install IntaRNA
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
COPY ViennaRNA-1.8.5.tar.gz /tmp/
RUN tar -C /tmp/ -zxf /tmp/ViennaRNA-1.8.5.tar.gz 
WORKDIR /tmp/ViennaRNA-1.8.5
RUN ./configure CFLAGS='-D inline=' --without-perl --without-forester --without-kinfold && make && make install

COPY intarna-1.2.5.tar.gz /tmp/
RUN tar -C /tmp/ -zxf /tmp/intarna-1.2.5.tar.gz
WORKDIR /tmp/intarna-1.2.5
RUN ./configure && make



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Install my script to perform basic SNP calling
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
COPY mplp2vcf.cpp /tmp/
RUN cd /tmp;g++ -O3 -Wall -std=c++0x -o /usr/local/bin/mplp2vcf -I/usr/local/include/htslib mplp2vcf.cpp -Lsrc/usr/htslib -Wl,-Bstatic -lhts -Wl,-Bdynamic -lz -llzma -lbz2 -lcurl -lpthread



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Bowtie
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL https://freefr.dl.sourceforge.net/project/bowtie-bio/bowtie/1.2.0/bowtie-1.2-linux-legacy-x86_64.zip -o /tmp/bowtie-1.2-linux-legacy-x86_64.zip && \
    unzip /tmp/bowtie-1.2-linux-legacy-x86_64.zip -d /tmp && \
    find /tmp/bowtie-1.2-legacy -maxdepth 1 -type f -executable -exec mv '{}' /usr/local/bin/ ';' && \
    rm -rf /tmp/bowtie-1.2-linux-legacy-x86_64.zip /tmp/bowtie-1.2-legacy



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Jellyfish
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz | tar -C /tmp -zxf - && \
    cd /tmp/jellyfish-2.3.0 && ./configure && make -j4 && make install && \
    rm -rf /tmp/jellyfish-2.3.0


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Fastx toolkit
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
RUN apt-get update && apt-get install -y fastx-toolkit

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Install Python
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
RUN apt-get update && apt-get install -y \
    build-essential \
    python2.7-dev \
    python-numpy \
    python-matplotlib \
    python-pip \
    libbz2-dev \
    liblzma-dev


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Cutadapt HTSeq-count, MACS2, HTSeq, umi_tools
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
RUN pip install --upgrade cutadapt==1.9
RUN pip install --upgrade pip && easy_install -U setuptools
#RUN pip install pysam HTSeq MACS2==2.1.3
#RUN pip install cython pandas future umi_tools

RUN curl -kL http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz | tar -C /tmp -xzf - && find /tmp/SPAdes-3.14.1-Linux/ -type f -executable -exec mv '{}' /usr/local/bin/ ';'


# GATK
# snpEff
# Mfold/UNAfold
# http://data.broadinstitute.org/igv/projects/downloads/igvtools_2.3.79.zip
# my Makefiles with NGS rules
#  minimap


# TOOLS WITH GUI
# * IGV
# * Cytoscape
# * SeqMonk
# * FastQC








VOLUME ["/export/"]
WORKDIR /export
ENTRYPOINT ["/bin/bash","-c"]
#EXPOSE :80



