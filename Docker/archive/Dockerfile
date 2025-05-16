FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# Dependencies

RUN apt-get update -yq && \
	apt-get install -yq --no-install-recommends \
	wget \
	curl \
	g++ \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev\
	libncurses5-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	libblas-dev \
	liblapack-dev \
	libgsl-dev \
	gfortran \
	libxml2-dev \
	libperl-dev \
	libpng-dev \
	locales \
	ca-certificates \
	parallel \
	software-properties-common \
	lsb-release && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install cmake

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | \
	gpg --dearmor - | \
	tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null && \
	apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" && \
	apt-get install -yq --no-install-recommends \
	cmake=3.25.2-0kitware1ubuntu18.04.1 \
	make \
	pkg-config && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set locale

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# Set paths

ENV OPT=/opt/wtsi-cgp
ENV PATH=$OPT/bin:$PATH
ENV LD_LIBRARY_PATH=$OPT/lib

# External tools: set tool versions for build script, add & run build script

ENV VER_SAMTOOLS="1.21"
ENV VER_HTSLIB="1.21"
ENV VER_BCFTOOLS="1.21"
ENV VER_VERIFYBAMID="2.0.1"
ENV VER_LIBDEFLATE="v1.23"

ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# Internal tools:
COPY Makefile .
COPY src /src
ADD build/opt-build-local.sh build/
RUN bash build/opt-build-local.sh $OPT

# User config

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu

# Entry point
CMD ["/bin/bash"]
