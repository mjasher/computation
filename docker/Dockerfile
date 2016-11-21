# -----------------------------------------------
# A container with the best of Python (OpenCV and GDAL)
# 
# build your own image:
# docker build -t py_cv_gdal .
#
# open bash in image
# docker run -it docker_image bash
#
# run a server on http://localhost:5000
# docker run -it -p 5000:5000 -v ~/share_with_docker:/share_with_docker py_cv_gdal  python /share_with_docker/flask_demo.py 
# deploy
# docker run -e PYTHONPATH=/share_with_docker -it -p 5000:5000 -v ~/share_with_docker:/share_with_docker py_cv_gdal /usr/local/bin/gunicorn flask_demo:app --bind 0.0.0.0:5000 -w 2 -k gevent --access-logfile -
# -----------------------------------------------

# -----------------------------------------------
# Install OpenCV 3.1
# thanks http://www.pyimagesearch.com/2016/10/24/ubuntu-16-04-how-to-install-opencv/
# -----------------------------------------------

FROM ubuntu:16.04

RUN apt-get -y update && \
    apt-get -y upgrade

RUN apt-get -y install \
    build-essential \
    cmake \
    pkg-config \
    libjpeg8-dev \
    libtiff5-dev \
    libjasper-dev \
    libpng12-dev \
    libavcodec-dev \
    libavformat-dev \
    libswscale-dev \
    libv4l-dev \
    libxvidcore-dev \
    libx264-dev \
    libgtk-3-dev \
    libatlas-base-dev \
    gfortran \
    python2.7-dev \ 
    python3.5-dev \
    wget \
    unzip

RUN cd ~ && \
    wget -O opencv.zip https://github.com/Itseez/opencv/archive/3.1.0.zip && \
    unzip opencv.zip && \
    wget -O opencv_contrib.zip https://github.com/Itseez/opencv_contrib/archive/3.1.0.zip && \
    unzip opencv_contrib.zip

RUN cd ~ && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py

RUN pip install numpy==1.11.1

RUN cd ~/opencv-3.1.0/ && \
    mkdir build && \
    cd build && \
    cmake -D CMAKE_BUILD_TYPE=RELEASE \
    -D CMAKE_INSTALL_PREFIX=/usr/local \
    -D INSTALL_PYTHON_EXAMPLES=ON \
    -D INSTALL_C_EXAMPLES=OFF \
    -D OPENCV_EXTRA_MODULES_PATH=~/opencv_contrib-3.1.0/modules \
    -D PYTHON_EXECUTABLE=/usr/bin/python \
    -D BUILD_EXAMPLES=ON .. && \
    make -j4 && \
    make install && \
    ldconfig

# -----------------------------------------------
# Install GDAL 2.1
# thanks https://github.com/GeographicaGS/Docker-GDAL2/blob/master/2.1.2/Dockerfile
# -----------------------------------------------

ENV ROOTDIR /usr/local/
ENV GDAL_VERSION 2.1.2

WORKDIR $ROOTDIR/

# RUN apt-get update -y && apt-get install -y python-numpy \
RUN apt-get install -y \
    python-dev \
    software-properties-common \
    python-software-properties \
    build-essential \
    libspatialite-dev \
    sqlite3 \
    libpq-dev \
    libcurl4-gnutls-dev \
    libproj-dev \
    libxml2-dev \
    libgeos-dev \
    libnetcdf-dev \
    libpoppler-dev \
    libspatialite-dev \
    libhdf4-alt-dev \
    libhdf5-serial-dev \
    wget

RUN wget http://download.osgeo.org/gdal/${GDAL_VERSION}/gdal-${GDAL_VERSION}.tar.gz --directory-prefix=$ROOTDIR/src/ && \
    cd src && tar -xvf gdal-${GDAL_VERSION}.tar.gz && cd gdal-${GDAL_VERSION} \
    && ./configure --with-python --with-spatialite --with-pg --with-curl \
    && make && make install && ldconfig \
    && apt-get update -y \
    && apt-get remove -y --purge build-essential wget \
    && rm -Rf $ROOTDIR/src/*

# Output version and capabilities by default.
# CMD gdalinfo --version && gdalinfo --formats && ogrinfo --formats

# -----------------------------------------------
# Install other python libraries
# -----------------------------------------------

RUN pip install \
    numpy==1.11.1 \
    scipy==0.18.0 \
    matplotlib==1.5.3 \
    Pillow==3.2.0 \
    requests==2.11.1 \
    scikit-image==0.12.3 \
    scikit-learn==0.18.0 \
    gevent==1.1.2 \
    gunicorn==19.6.0 \
    Flask==0.11.1 \
    boto3==1.4.1
