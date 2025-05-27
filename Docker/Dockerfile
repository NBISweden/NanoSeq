# Build stage 1: src, samtools and dependencies

FROM ubuntu:18.04 AS buildsrc
ENV DEBIAN_FRONTEND=noninteractive

## Install apt dependencies

RUN apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	bzip2 \
	ca-certificates \
	gcc \
	g++ \
	gfortran \
	libblas-dev \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgsl-dev \
	liblapack-dev \
	liblzma-dev \
	libncurses5-dev \
	libperl-dev \
	libpng-dev \
	libssl-dev \
	libxml2-dev \
	locales \
	lsb-release \
	make \
	parallel \
	pkg-config \
	software-properties-common \
	wget \
	zlib1g-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Apt continued, install cmake to build libdeflate

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | \
	gpg --dearmor - | \
	tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null && \
	apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" && \
	apt-get install -yq --no-install-recommends \
	cmake=3.25.2-0kitware1ubuntu18.04.1 \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Environment

ENV OPT=/usr/local
ENV PATH=$OPT/bin:$PATH
ENV LD_LIBRARY_PATH=$OPT/lib

## Locale

RUN locale-gen en_US.UTF-8 && \
	update-locale LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

## Software versions

ENV SAMTOOLS_ECOSYSTEM="1.18"

## Build libdeflate (for htslib)

WORKDIR /tmp
RUN wget https://github.com/ebiggers/libdeflate/archive/refs/tags/v$SAMTOOLS_ECOSYSTEM.tar.gz -O libdeflate.tar.gz && \
	mkdir -p libdeflate && \
	tar --strip-components=1 -C libdeflate -zxf libdeflate.tar.gz && \
	cd libdeflate && \
	cmake -B build && \
	cmake --build build && \
	cmake --install build && \
	cmake --install build --prefix $OPT && \
	cd /tmp && \
	rm -rf libdeflate libdeflate.tar.gz

## Build htslib (for source code makefile)

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/$SAMTOOLS_ECOSYSTEM/htslib-$SAMTOOLS_ECOSYSTEM.tar.bz2 && \
	mkdir -p htslib && \
	tar --strip-components=1 -C htslib -jxf htslib-$SAMTOOLS_ECOSYSTEM.tar.bz2 && \
	cd htslib && \
	export CFLAGS="-I$OPT/include" && \
	export LDFLAGS="-L$OPT/lib" && \
	./configure --disable-plugins --enable-libcurl --with-libdeflate --prefix=$OPT && \
	make -j$(nproc) && \
	make install && \
	mkdir -p $OPT/include/cram && \
	cp ./cram/*.h $OPT/include/cram/ && \
	cp header.h $OPT/include && \
	cd /tmp && \
	rm -rf htslib htslib-$SAMTOOLS_ECOSYSTEM.tar.bz2

## Build samtools (used alongside compiled code)

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_ECOSYSTEM/samtools-$SAMTOOLS_ECOSYSTEM.tar.bz2 -O samtools.tar.bz2 && \
	mkdir -p samtools && \
	tar --strip-components=1 -C samtools -xjf samtools.tar.bz2 && \
	cd samtools && \
	./configure --with-htslib=$OPT --enable-plugins --enable-libcurl --prefix=$OPT && \
	make -j$(nproc) && \
	make install && \
	cd /tmp && \
	rm -rf samtools samtools.tar.bz2

## Build gzstream (for compiled code)

WORKDIR /tmp
RUN wget https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz -O gzstream.tgz && \
	mkdir -p gzstream && \
	tar --strip-components=1 -C gzstream -xzf gzstream.tgz && \
	cd gzstream && \
	cp *.h $OPT/include && \
	make && \
	cp libgzstream.a $OPT/lib && \
	cd /tmp && \
	rm -rf gzstream gzstream.tgz

## Compile source code

WORKDIR /tmp
COPY Makefile .
COPY src ./src/
RUN make PREFIX=$OPT && \
	make install PREFIX=$OPT && \
	cd /tmp && \
	rm -rf src Makefile

# Build stage 2: R packages

FROM ubuntu:18.04 AS buildr
ENV DEBIAN_FRONTEND=noninteractive

## Install apt dependencies

RUN apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	bzip2 \
	ca-certificates \
	libblas-dev \
	libbz2-dev \
	libcurl4-openssl-dev \
	liblapack-dev \
	liblzma-dev \
	libpng-dev \
	libssl-dev \
	libxml2-dev \
	locales \
	make \
	pkg-config \
	software-properties-common \
	wget \
	zlib1g-dev \
	&& add-apt-repository ppa:ubuntu-toolchain-r/test && \
	apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	gcc-10 \
	gfortran-10 \
	g++-10 \
	&& update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 \
	&& update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100 \
	&& update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-10 100 \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Apt continued, install R

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
	tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
	add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
	apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	r-base-core=4.1.3-1.1804.0 && \
	apt-mark hold r-base-core && \
	apt-get install -yq --no-install-recommends \
	r-cran-mass=7.3-51.5-2bionic0 \
	r-cran-class=7.3-16-1bionic0 \
	r-cran-nnet=7.3-13-1bionic0 \
	r-recommended=4.1.3-1.1804.0 \
	r-base=4.1.3-1.1804.0 && \
	apt-mark hold r-base r-recommended \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Environment

ENV OPT=/usr/local
ENV PATH=$OPT/bin:$PATH
ENV LD_LIBRARY_PATH=$OPT/lib
ENV R_LIBS=$OPT/R-lib
ENV R_LIBS_USER=$R_LIBS

## Locale

RUN locale-gen en_US.UTF-8 && \
	update-locale LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# Build R packages

RUN mkdir -p $R_LIBS && \
	R -e "options(repos = list(CRAN = 'https://cloud.r-project.org')); install.packages(c('BiocManager'), lib = Sys.getenv('R_LIBS'))"
RUN R -e "BiocManager::install(c('deepSNV', 'vcfR', 'VGAM'), ask = FALSE, lib = Sys.getenv('R_LIBS'))"
RUN R -e "options(repos = list(CRAN = 'https://cloud.r-project.org')); install.packages(c('ggplot2', 'data.table', 'epitools', 'gridExtra', 'seqinr'), lib = Sys.getenv('R_LIBS'))"
RUN R -e 'pkgs <- c("ggplot2", "data.table", "epitools", "gridExtra", "seqinr", "deepSNV", "vcfR", "VGAM"); sapply(pkgs, function(pkg) library(pkg, character.only=TRUE, lib.loc=Sys.getenv("R_LIBS")))'

# Runtime image

FROM ubuntu:18.04 AS runtime
ENV DEBIAN_FRONTEND=noninteractive

## Metadata

LABEL org.opencontainers.image.description="Docker image to run NanoSeq source codebase" \
	org.opencontainers.image.version="0.94" \
	org.opencontainers.image.authors="cormac.kinsella@nbis.se"

## Install apt dependencies

RUN apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	dirmngr \
	libcurl4-openssl-dev \
	libperl5.26 \
	libcapture-tiny-perl \
	libfile-which-perl \
	locales \
	software-properties-common \
	wget \
	&& add-apt-repository ppa:ubuntu-toolchain-r/test && \
	apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	libstdc++6 \
	&& wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
	tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
	add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
	apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	r-base-core=4.1.3-1.1804.0 && \
	apt-mark hold r-base-core && \
	apt-get install -yq --no-install-recommends \
	r-cran-mass=7.3-51.5-2bionic0 \
	r-cran-class=7.3-16-1bionic0 \
	r-cran-nnet=7.3-13-1bionic0 \
	r-recommended=4.1.3-1.1804.0 \
	r-base=4.1.3-1.1804.0 && \
	apt-mark hold r-base r-recommended \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
	&& apt-get purge -y --auto-remove software-properties-common dirmngr wget

## Environment

ENV OPT=/usr/local
ENV PATH=$OPT/bin:$PATH
ENV LD_LIBRARY_PATH=$OPT/lib
ENV R_LIBS=$OPT/R-lib
ENV R_LIBS_USER=$R_LIBS

## Locale

RUN locale-gen en_US.UTF-8 && \
update-locale LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# Copy binaries

COPY --from=buildsrc $OPT $OPT
COPY --from=buildr $R_LIBS $R_LIBS

## Configure user entrypoint

RUN adduser --disabled-password ubuntu && \
	mkdir -p /home/ubuntu && \
	chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
