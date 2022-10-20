FROM ubuntu:22.04
	
	
RUN apt-get update --fix-missing \
&& apt-get install -y wget build-essential bzip2 cmake default-jdk git libnss-sss libtbb2 libtbb-dev ncurses-dev unzip zlib1g zlib1g-dev \
&& cd /usr/local/ \
&& wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz \
&& tar -xzvf 2.7.10a.tar.gz \
&& rm -rf 2.7.10a.tar.gz  \
&& cp -r STAR-2.7.10a /usr/local/bin \
&& chmod +x /usr/local/bin/STAR-2.7.10a \
&& rm -rf 2.7.0f.tar.gz  \
&& apt-get autoremove -y  \
&& apt-get autoclean -y \
&& apt-get clean

ENTRYPOINT ["/usr/local/bin/STAR-2.7.10a"]
	
	