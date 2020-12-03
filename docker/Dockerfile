FROM continuumio/miniconda3:4.8.2-alpine AS builder

ENV PATH /opt/conda/bin:$PATH

USER root

RUN apk update && apk upgrade && apk add bash cmake g++ gfortran perl make vim nano && \
	mkdir /workspace && \
	mkdir /software && \
	cd /software && \
    wget http://ftp.gromacs.org/pub/gromacs/gromacs-2020.4.tar.gz && \
	tar xzf gromacs-2020.4.tar.gz && cd gromacs-2020.4 && \
	mkdir build && cd build && \
	cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON && \
	make -j 4 && make check && make install && \
	echo "source /usr/local/gromacs/bin/GMXRC" >> /home/anaconda/.bashrc && \
	conda update conda -y && \
	. /opt/conda/etc/profile.d/conda.sh && \
    conda create --name mc mosdef_cassandra foyer matplotlib seaborn pycifrw garnett \
    -c conda-forge -c mosdef -c omnia && \
	echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/anaconda/.bashrc && \
	echo "conda activate base" >> /home/anaconda/.bashrc && \
	echo "conda activate mc" >> /home/anaconda/.bashrc && \
	conda clean -afy

RUN apk add git && \
	. /opt/conda/etc/profile.d/conda.sh && \
	cd /software && \
    git clone https://github.com/rsdefever/mc_examples.git && \
	cd mc_examples && \
	conda activate mc && \
	pip install .


WORKDIR /workspace

COPY entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh
RUN chown -R anaconda /workspace
RUN chmod 755 /workspace

USER anaconda

ENTRYPOINT ["/entrypoint.sh"]
CMD ["none"]
