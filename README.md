reportRx installation and tutorial
=

reportRx ~ report + R(language) + R(yan) + x(docx) + Rx (prescription (clinical)) + more
-

[If you have any concerns comments or suggestions please edit the wiki here by clicking 'Edit Page'](https://github.com/rdelbel/reportRx/wiki/Suggestions)

[Email rdelbel at gmail dot com to report any bugs with the corresponding .Rnw file and data files as attachments](mailto:rdelbel@gmail.com)

reportRx relies heavily on knitr. First we must install knitr and all its dependencies.

First we install knitr in R
```github
install.packages('knitr')
```

Knitr does not need to use LaTeX. Markdown might actually be a better choice as it is much easier to use. [(In fact this document was written in markdown.)](https://raw.github.com/rdelbel/reportRx/master/README.md) We will use LaTeX though as it works better in our ecosystem

[Download MikTeX here](http://miktex.org/download)

Now we need to setup RStudio to use knitr. Tools->Options->Sweve->Weve Rnw files-> knitr. Also make sure 'Typeset LaTeX into pdf using' is set to pdfLaTeX.

In order to use LaTeX & knitr in RStudio we need to use .Rnw files. We can open them in RStudio by going to File->New->R Sweve. Note that selecting R Markdown would produce a knitr + markdown (Rmd) file.


We are now ready to use knitr to knit LaTeX and R. Although knitr is not interactive in the sense that we have to compile the pdf, we can still interact in the R session inside of the R chunks. The output will appear in the console as normal. Note that we can create a new chunk using ```ctrl+alt+i```. Also after placing an initial comma in the ```<<>>=``` to get ```<<,>>=``` knitr will autocomplte all of the options. You can move around with the arrows and select with the enter key.  [You should go though my short knitr tutorial before continuing. Just copy and paste this file into a blank .Rnw file in RSudio.](https://raw.github.com/rdelbel/reportRx/master/Vignettes/ClinicalReportVignette.Rnw) There is a big button on top of the RStudio editor to compile the .Rnw file to pdf. You should read the source and view the resulting pdf to see what is going on. Note that the compiling process creates a lot of files. It is best to make a new folder for each .Rnw file you use. 

[The knitr website had surprisingly good documentation](http://yihui.name/knitr/options)

reportRx
-

reportRx is just a set of functions that do (hopefully) useful things in a .Rnw file. The 'new R Sweve (.Rnw)' file does not contain everything we will need to use reportRx. Unfortunately there is no way to change the file you get when you click 'New R Sweve'. I recommend saving this file opening it every time you want to start a new document that will use reportRx. You can then 'save as' to your desired file name and location. The very last part should be modified for your personal situation which will be discussed later.

```TeX
\documentclass{article}
\usepackage{multirow}
\setlength\parindent{0pt}
\usepackage{geometry}
\usepackage{longtable}
\usepackage{float}

%This is word 'normal' margins
\geometry{left=1.25in,right=1in,top=1.25in,bottom=1.25in}
%word 'moderate' margins
%\geometry{left=.75in,right=.75in,top=1in,bottom=1in}
%word 'narrow' margins
%\geometry{left=.5in,right=.5in,top=.5in,bottom=.5in}



%Change title here
\title{reportRx}
%Change Author here
\author{Ryan Del Bel}

\begin{document}
\maketitle

<<,include=FALSE>>=
require(reportRx)
@

%Do what you want here


<<,include=F>>=
#makedocx("Directory of file",
#          "Name of file with no extension",
#          "Directory of pandoc",
#          "Directory of imagemagick")
@

\end{document}
```
Installing reportRx is slightly more complicated than a normal package. It is not currently on CRAN (the main repository for R packages), but it is on github. We first need to install the package devtools, and then we can use a function from devtools to install my package off github. We use the :: notation so we don't have to require the devtools package.

```R
install.packages('devtools')
devtools::install_github('reportRx','rdelbel')
```
`devtools::install_github('reportRx','rdelbel')` can also be used to update the package.

This is all we need to use reportRx and create .pdf documents. If you want to create a .docx document then we will need to install more software to work in the command line. These three programs are only used in the `reportRx::makedocx` function.

* [pandoc: convert .tex to .docx](https://code.google.com/p/pandoc/downloads/list)
* [ghostscrip: pandoc dependency](http://www.ghostscript.com/download/gsdnld.html)
* [imagemagick: convert image formats](http://www.imagemagick.org/script/binary-releases.php)

**Install pandoc and imagemagick in a directory with no spaces** as we will call them from the command line and spaces in the command line are always a pain. I recommend downloading the ```Win64 dynamic at 16 bits-per-pixel``` version of imagemagick. ghostscript is not needed on Ubuntu and probably not needed on OSX. Installing pandoc on Ubuntu 13.04 is somewhat complicated. Contact me if you want to do it.


[You should go through the reportRx short tutorial now.](https://raw.github.com/rdelbel/reportRx/master/Vignettes/ClinicalReportVignette.Rnw)

reportRx has several functions not in the tutorial, and the functions in the tutorial also may have more options not specified. To view all of the functions in reportRx, in RStudio write `reportRX::<tab>` where ```<tab>``` is pressing the tab button. This will give a drop-down list of all of the functions. You can go to the function you want by using the arrow keys. After hovering over the function you want press `<enter>(<tab>` to view all of the parameters and the documentation for each parameter. While in this menu you can press ```f1``` at any time to have the full documentation of the function pop up.

You can spell-check your document in RStudio by pressing ```f7```. Once you are happy with your .pdf file you might want to convert it to .docx. If the above software is installed properly we can do it as follows **on windows only** (for now)
```R
makedocx("Directory of file",
          "Name of file with no extension",
          "Directory of pandoc",
          "Directory of imagemagick")
```
The 4th argument is optional and only required when there is at least one image file in the report. For example on my computer I might use the following command

```R
makedocx("C:\\pm\\presentations\\reportRx\\",
           "ClinicalReportVignette",
           "C:\\pm\\Pandoc",
           "C:\\pm\\ImageMagick")
```

It is important to note that we can not run this function while we are knitting the pdf. This will try and convert a .tex file that does not exist yet (or is an old version). You may have noticed in the template that the bottom had the following code

```R
<<,include=F>>=
#makedocx("Directory of file",
#          "Name of file with no extension",
#          "Directory of pandoc",
#          "Directory of imagemagick")
@
```

Since there we have specified ```include=F``` this code will never show up in our documents, and since we have commented out the code nothing will ever run. After we have created the .tex file if we can uncomment and run this code to convert it to .docx. Remember to comment it again so that it does not run the next time you compile a pdf. Note that usually your directory of pandoc and imagemagick will not move, and likely the directory of the file is your working directory. Thus we only really need to specify the name of the file each time if we use the following code.

```R
<<,include=F>>=
#makedocx( getwd(),
#          "Name of file with no extension",
#          "your filled in directory for pandoc",
#          "your filled in directory for imagemagick")
@
```
Of course getwd() will only work if the working directory is the directory of the .Rnw file. 

As explained in the tutorial, we have to make scarifies when creating a document that we will eventually convert to .docx. In addition to needing to use an inferior package to make tables, we also can not add color (hence having bold pvalues), and **we can not automatically number our ```\section{}```.** We have to use ```\section*{}``` and number it ourselves.

After we have created the .docx document we will notice that there are some formatting problems. In particular there are [H] everywhere, our tables take up too much space, and the images are not aligned. We can fix all of these problems with the [microsoft word macro here](https://raw.github.com/rdelbel/reportRx/master/wordmacro.txt). Note that there are some instructions at the top.

reportRx functions
-

There are three main sets of functions:

* Functions that output LaTeX 
* Helper functions that make working with LaTeX inside of knitr easier
* GWAS QC functions. 

The functions that output to LaTeX will automatically parse the colnames of your dataframe to turn ```.``` and ```_``` into spaces. Usually this is what we want, however sometimes it is not. At the moment there is no way to override this. Some of the more useful helper functions are ```sanitizestr```, ```lpvalue```, ```nicenamee```, and ```lbld```. You can view the help files for their documentation. I suggest talking to me before using the GWAS QC functions, and being careful when using the ```citime``` and ```plotci``` functions as they are not that well documented and may do things you don't expect.

Citation
-
To cite the 'reportRx' package in publications use:

Ryan Del Bel, Wei Xu (2013). reportRx: A package for clinical report generation in R. R package version 1.0.0

Disclaimer
-


reportRx has been tested by several people and there are no known cases of it producing the wrong estimates, CI, or pvalues. It is possible that there is a bug for some edge-cases that have not been encountered yet. If something looks funny it is best to fit a model by hand and see if the answer is the same. To ease your heart it might be good to fit the main model of the paper by hand to just make sure it is correct. Like R, this package comes with no warranty and I take no responsibility if it breaks science.

If you do encounter a bug [email rdelbel at gmail dot com. Attach the corresponding .Rnw and data files.](mailto:rdelbel@gmail.com)


Known bugs
-

* We should not even have to use imagemagick to convert the image files because knitr can output any extension we want. The images will only show up in the .docx if we convert them via imagemagick. So strange.
* It is possible that the pmvsum function will not work in some cases with ambiguous factor variable and level names. In this case the function will give you a cryptic warning and end gracefully producing no output. This has never occurred in practice but will be fixed eventually.
* When converting to .docx the images will not be that sharp
* The pcovsum function right justify columns of the table
* petsum will produce NA if you ask for a survival probability past the latest recorded time. Although not totally wrong it would be nice if the function handled it more gracefully
* You may sometimes notice that if you update your plots they will not update in the .pdf file. This will only happen if you have already converted to .docx. If you convert the new .pdf to .docx the new plots will show up. If you want to see the new plots in the .pdf then you will have to delete the plots in the /figures folder. This happens because we need to convert the .pdf images to .png images when we convert to .docx. Unfortunately when we subsequently compile the next .pdf it will take the old .png images instead of the new .pdf images. I do not know enough about LaTeX to fix this at the moment.