# base image : ubuntu
FROM ubuntu:22.04
RUN apt-get update --fix-missing \
&& apt-get install -y wget \
&& cd /usr/local/ \
&& wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz/download\
&& tar xzvf 2.7.10a.tar.gz \
&& cp -r STAR-2.7.10a /usr/local/bin \
&& chmod +x /usr/local/bin/STAR-2.7.10a \
&& rm 2.7.10a.tar.gz   \
&& rm -rf STAR-2.7.10a

ENTRYPOINT ["/usr/local/bin/STAR-2.7.10a"]