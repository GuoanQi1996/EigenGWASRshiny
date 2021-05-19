# EigenGWAS Shiny
The core algorithm of this online platform is implemented in efficient C language and provides a simple and interactive user-friendly interface using R shiny framework.  

Online tool is now available at www.eigengwas.com.  

You can easily deploy EigenGWAS directly at you PC/sever by the followings.  
###   -Download the source codes of EigenGWAS Shiny
Find the download link at the homepage **Code** button or **Releases** panel, download the zipped source codes.
###   -Unzip the source codes locally
###   -Open R/Rstudio
Run the commands below to initialize EigenGWAS Shiny.
~~~
# This is an R console

# R packages shiny, bsplus and rmarkdown is required in EigenGWAS
# Install them if they are not all set in your system
install.packages(c("shiny","bsplus","rmarkdown"))

# Before start the app, be sure to set the directory of zipped source codes, specifically the file contains the "app.R", as the working directory
# For example if you unzipped the source code at /home/your_name/EigenGWAS/, then you should run
setwd("/home/your_name/EigenGWAS/")

library(shiny)
runApp()
~~~
Normally for most of the users working with Windows/MacOS/Unbuntu etc., there should be a window or browser tab pops up, EigenWGAS is then ready for the analysis.
For some of the linux OS without displays, LAN remote access can be used for the EigenGWAS Shiny.
~~~
# This is an R console

# Assume your linux IP is 100.100.100.1, and port 1234 is accessible
runApp(host='100.100.100.1',port=1234)
~~~
EigenGWAS is then availble at computer in LAN, by visit 100.100.100.1:1234 in the browser.

The way to open the specific port (for example 1234) in linux:
~~~
# This is a linux bash/terminal

# make sure the firewall is active
systemctl status firewalld
systemctl start firewalld

# open port 1234
firewall-cmd --zone=public --add-port=1234/tcp --permanent

# note
# --zone: scope zone
# --permanent: set options permanently, a change will only be part of the runtime configuration without this option
~~~

All platforms support.