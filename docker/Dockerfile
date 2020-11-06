### author: Sebastian Leifeld

FROM nvcr.io/nvidia/nvhpc:20.9-devel-centos7

RUN yum update -y
RUN yum remove -y cmake
RUN yum install -y make cmake3 python python-numpy -y bc libxml2-utils xsltproc tidy python3 python3-numpy vim
RUN ln -s /usr/bin/cmake3 /usr/bin/cmake
