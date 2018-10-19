FROM ubuntu:bionic

LABEL maintainer="Bertrand Neron"
LABEL maintainer.email="bneron@pasteur.fr"
LABEL authors="Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha"
LABEL package.version="2.0"
LABEL package.license="GPLv3"
LABEL package.homepage="https://github.com/gem-pasteur/Integron_Finder/"
LABEL publications="https://doi.org/10.1093/nar/gkw319"


ARG DEBIAN_FRONTEND=noninteractive
ARG IF_VERS="2.0rc4"

RUN apt-get update -y &&\
    apt-get install -y --no-install-recommends python3.6  python3-pip git hmmer infernal prodigal &&\
    apt-get install -y python3-numpy python3-pandas &&\
    apt-get clean -y

RUN pip3 install setuptools &&\
    pip3 download --index-url https://test.pypi.org/simple/ --no-deps integron_finder==${IF_VERS} &&\
    pip3 install "integron_finder-${IF_VERS}.tar.gz"

ENV DEBIAN_FRONTEND teletype
ENV PYTHONIOENCODING UTF-8

ENTRYPOINT ["/usr/local/bin/integron_finder"]

CMD ['-h']