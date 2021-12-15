### author: Sebastian Leifeld

FROM nvcr.io/nvidia/nvhpc:20.11-devel-cuda11.1-ubuntu20.04

RUN apt update
RUN apt install -y make cmake
RUN apt install -y python3 python3-numpy
RUN apt install -y bc libxml2-utils xsltproc tidy vim bash-completion
RUN apt install -y libtinfo5 libncursesw5 valgrind

RUN echo "source /etc/profile.d/bash_completion.sh" >> /root/.bashrc
RUN echo "source /etc/profile.d/bash_completion.sh" >> /.bashrc
