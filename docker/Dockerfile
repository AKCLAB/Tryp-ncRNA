FROM continuumio/miniconda3

#Work directory /
WORKDIR /

RUN apt-get update && apt-get install -y \
    perl \
    gcc \
    gfortran \
    wget \
    unzip \
    build-essential \
    && apt-get clean 

# Create and activate tryp-ncrna enviroment
RUN conda create -y -n tryp-ncrna python=3.12 r-base=4.3.1 r-essentials

SHELL ["conda", "run", "-n", "tryp-ncrna", "/bin/bash", "-c"]

RUN conda install -c conda-forge -c bioconda bowtie2=2.5.4 samtools=1.21 bedtools=2.31.1 diamond=2.1.11 igvtools=2.3.93 snoscan openjdk=11 
RUN conda install -c conda-forge cmake dos2unix
RUN conda install bioconda::picard

RUN ln -s /opt/conda/envs/tryp-ncrna/bin/dos2unix /usr/bin/dos2unix && \
    ln -s /opt/conda/envs/tryp-ncrna/bin/Rscript /usr/bin/Rscript

RUN pip install regex xlsxwriter matplotlib seaborn scipy pysam fuzzywuzzy biopython==1.77 cutadapt

#Download Tryp-ncrna pipeline
#RUN git clone https://github.com/AKCLAB/Tryp-ncRNA.git

#Install PORTRAIT
RUN mkdir programs && \
    cd programs && \
    wget https://www.bioinformatics.org/portrait/download/portrait.tar.gz && \
    tar -zxvf portrait.tar.gz && \
    wget https://www.bioinformatics.org/portrait/download/angle.tar.gz && \
    tar -zxvf angle.tar.gz && \
    wget https://www.bioinformatics.org/portrait/download/cast-linux.tar.gz && \
    tar -zxvf cast-linux.tar.gz && \
    chmod +x cast-linux && \    
    wget https://www.bioinformatics.org/portrait/download/libsvm-2.84.tar.gz && \
    dpkg --add-architecture i386 && \
    apt update && \
    apt install -y libc6:i386 libstdc++6:i386 && \
    apt install -y libstdc++5 && \
    tar -zxvf libsvm-2.84.tar.gz  && \
    rm *.gz && \
    cd libsvm-2.84  && \
    make 
#Install tRNAscan-SE
RUN cd programs && \
    git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git && \
    cd tRNAscan-SE && \
    ./configure && \
    make install && \
    apt install -y infernal && \
    ln -s /usr/bin/cmsearch /usr/local/bin/cmsearch && \
    ln -s /usr/bin/cmscan /usr/local/bin/cmscan
#Install ptRNApred
RUN cd programs && \
    wget http://www.ptrnapred.org/ptRNApred1.0.zip && \
    unzip ptRNApred1.0.zip && \
    wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz && \
    tar -zxvf ViennaRNA-2.7.0.tar.gz && \
    cd ViennaRNA-2.7.0 && \
    ./configure --without-perl --without-python --without-python2 --without-forester --without-rnalocmin && \
    make && make install && \
    cd .. && rm -rf ViennaRNA-2.7.0 ViennaRNA-2.7.0.tar.gz
# Install RNAcon
RUN cd programs && \
    wget http://webs.iiitd.edu.in/raghava/rnacon/RNAcon_v1.0.tar.gz && \
    tar -zxvf RNAcon_v1.0.tar.gz && \
    wget https://osmot.cs.cornell.edu/svm_light/current/svm_light_linux64.tar.gz && \
    tar -zxvf svm_light_linux64.tar.gz && \
    mv svm_classify RNAcon_standalone/progs && \
    wget https://prdownloads.sourceforge.net/weka/weka-3-8-6-azul-zulu-linux.zip && \
    unzip weka-3-8-6-azul-zulu-linux.zip && \
    chmod +x weka-3-8-6/weka.jar && \
    mv weka-3-8-6/weka.jar RNAcon_standalone/progs && \
    wget https://github.com/satoken/ipknot/releases/download/v1.1.0/ipknot-1.1.0-x86_64-linux.zip && \
    unzip ipknot-1.1.0-x86_64-linux.zip && \
    mv ipknot-1.1.0-x86_64-linux/ipknot RNAcon_standalone/progs && \
    wget https://varna.lisn.upsaclay.fr/bin/VARNAv3-8-src.jar && \
    chmod +x VARNAv3-8-src.jar && \
    mv VARNAv3-8-src.jar RNAcon_standalone/progs

RUN Rscript -e "install.packages('igraph', repos='https://cloud.r-project.org')" && \
    Rscript -e "install.packages('cpp11', repos='https://cloud.r-project.org')" && \
    Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install('Biostrings')"
 
#Define the conda enviroment to start 
#CMD ["conda", "run", "-n", "tryp-ncrna", "/bin/bash"]
