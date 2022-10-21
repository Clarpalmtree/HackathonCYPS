FROM ubuntu:22.04
	
	ARG STAR_VERSION=2.7.10a
	
	ENV PACKAGES gcc g++ make wget zlib1g-dev 
	
	RUN set -ex
	
	RUN apt-get update && \
	    apt-get install -y --no-install-recommends ${PACKAGES} && \
	    apt-get clean && \
	    g++ --version && \
	    cd /home && \
	    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz && \
	    tar -xzvf ${STAR_VERSION}.tar.gz && \
	    cd STAR-${STAR_VERSION}/source && \
	    make STARstatic && \
	    mkdir /home/bin && \
	    cp STAR /home/bin && \
	    cd /home && \
	    'rm' -rf STAR-${STAR_VERSION} && \
	    apt-get --purge autoremove -y  ${PACKAGES}
	
	ENV PATH /home/bin:${PATH}

	

