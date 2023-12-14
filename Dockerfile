FROM ubuntu:16.04

RUN apt-get update && apt-get -y upgrade && apt-get -y install vim sqlite3 python2.7-dev python-pip wget libatlas3-base zlib1g-dev agrep libcurl3-gnutls libncurses5-dev libncursesw5-dev libagg-dev ssh openssh-client\
    && rm -rf /var/lib/apt/lists/\
    && wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh\
    && mkdir /root/.conda\
    && bash Miniconda2-latest-Linux-x86_64.sh -b\
    && rm -f Miniconda2-latest-Linux-x86_64.sh


ENV PATH="/root/miniconda2/bin:${PATH}"
RUN conda install -y -c bioconda samtools=1.6 \
    && conda install -c anaconda biopython=1.72 \
    && pip install Pillow==5.0.0 regex==2019.3.12 joblib==0.13.2 pandas==0.24.2 
    

RUN apt-get update && \
    apt-get install -y perl \
                       bioperl \
                       default-jre

RUN mkdir -p /media
RUN mkdir -p /media/Pipelines
RUN mkdir -p /media/Data

COPY jscb_src /media/Pipelines/jscb_src
WORKDIR /media/Pipelines/jscb_src

ENTRYPOINT ["python","jscb.py"]
CMD ["0"]


