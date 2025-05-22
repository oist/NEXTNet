FROM ghcr.io/phusion/holy-build-box/hbb-64:rel-3.0.5
RUN mkdir -p /ext/boost/tmp && cd /ext/boost/tmp && \
    curl -O -L https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz && \
    tar xzf boost_1_86_0.tar.gz && cd boost_1_86_0 && \
    source /hbb_exe/activate && \
    ./bootstrap.sh && ./b2 headers && \
    mv boost/* ../../ && cd ../../ && rm -r tmp
