FROM fredhutch/r-shiny-base:4.4.2

RUN apt-get --allow-releaseinfo-change update -y

RUN apt-get install -y curl

RUN R -e "install.packages(c('tidyr', 'shiny', 'dplyr', 'ggplot2', 'readr', 'gridExtra', 'DT', 'plotly'))"

ADD ./app.R /app/
ADD ./data_folder /app/data_folder
ADD ./www /app/www
ADD ./check.R /tmp


WORKDIR /app

# make sure all packages are installed
# because R does not fail when there's an error installing a package.
RUN R -f /tmp/check.R --args tidyr shiny dplyr ggplot2 readr gridExtra DT plotly

EXPOSE 3838

CMD R -f app.R
