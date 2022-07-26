FROM python:3.7

RUN  	apt-get update && apt-get upgrade -y && \ 
	apt-get install -y software-properties-common && \ 
	add-apt-repository ppa:team-xbmc/ppa && \
	apt-get install -y git &&  \
	apt-get install -y git-lfs && \
	apt-get install -y nlohmann-json3-dev && \
	apt-get install -y ninja-build && \
	apt-get autoremove && \
	python -m pip install  --no-cache-dir  numpy==1.20.1 pandas==1.2.2 && \
	python -m pip install  --no-cache-dir  PyYAML plotly==5.9.0 SALib==1.4.5 &&\
	python -m pip install --no-cache-dir  kaleido==0.2.1 seaborn==0.11.1 
	

COPY cmake-3.14.0-Linux-x86_64.sh cmake.sh 

RUN chmod +x cmake.sh && ./cmake.sh --skip-license --prefix=/usr/ 

RUN git clone https://github.com/koparasy/HPAC.git -b Puppeteer && cd HPAC && git lfs pull && rm -rf .git

COPY puppet_env.sh /HPAC/puppet_env.sh

COPY smallRelease.cmake /HPAC/smallRelease.cmake
COPY setup.sh /HPAC/setup.sh  

RUN cd /HPAC && chmod +x setup.sh && ./setup.sh /PUPPET_INSTALL 4 && rm -rf /HPAC/build_compiler/ && rm -rf /HPAC/build_hpac 

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]

