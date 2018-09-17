FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("destiny")'

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'devtools::install_cran("FateID")'

RUN R -e 'devtools::install_cran("RaceID")'

LABEL version 0.1.2

ADD . /code

ENTRYPOINT Rscript /code/run.R
