# base image : ubuntu
FROM ubuntu:22.04
RUN apt-get update --fix-missing \
&& apt-get install -y wget \
&& cd /usr/local/ \
&& wget https://bioconductor.org/packages/release/bioc/src/contrib/DESeq2_1.36.0.tar.gz \
&& tar xzvf DESeq2_1.36.0.tar.gz  \
&& cp -r DESeq2 /usr/local/bin \
&& chmod +x /usr/local/bin/DESeq2 \
&& rm DESeq2_1.36.0.tar.gz   \
&& rm -rf DESeq2

ENTRYPOINT ["/usr/local/bin/DESeq2"]