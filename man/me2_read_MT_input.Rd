% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/me2_read_MT_input.R
\name{me2_read_MT_input}
\alias{me2_read_MT_input}
\title{Get the multi time data}
\usage{
me2_read_MT_input(file, tz = "Etc/GMT-1")
}
\arguments{
\item{file}{File name of the multi time data file.}

\item{tz}{Time zone for date time, defaults to \code{Etc/GMT-1}}
}
\value{
list containing full data, both concentrations and uncertainties.
}
\description{
The data from the multi time data file can be read using this function. It
also performs some clean up, allowing the use of data in R.
}
\section{Data format}{

The format of the multi time data has to consists of the following columns,
in the same order as given in the table.

\tabular{ll}{
\strong{Name}    \tab \strong{Description}\cr
Datestart:       \tab Start date of the measurement in the 'yyyy-mm-dd' or 'dd-mm-yyyy' format\cr
Time             \tab Start time of the measurement in the 'hh:mm' format\cr
Dateend:         \tab End date of the measurement in the 'yyyy-mm-dd' or 'dd-mm-yyyy' format\cr
Time             \tab End time of the measurement in the 'hh:mm' format\cr
Tzone            \tab Timezone associated with the start and end date times\cr
Begin            \tab Unique number, expressing the beginning of the sample in time-units used in the data\cr
Length           \tab Unique number, expressing the length of the sampling time in time-units used in the data\cr
End              \tab Unique number, expressing the end of the sample in time-units used in the data\cr
...              \tab Columns of species with their appropriate names\cr
}

For more information on the multi time data format see, for example,
\href{https://doi.org/10.1016/j.scitotenv.2022.157981}{Mooibroek et al.(2022)}
}

