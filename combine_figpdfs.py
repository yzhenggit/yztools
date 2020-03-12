from PyPDF2 import PdfFileWriter, PdfFileReader
import os

def combine_figpdfs(figdir, outpdfname='new.pdf', outputdir='.'):
    '''
    Combine all the pdf in the figdir into one pdf file. 
    '''

    # look for the pdfiles in the figdir
    pdffiles=os.listdir(figdir)
    pdffiles.sort()
    if len(pdffiles) == 0:
        print('Empty ', figdir)
        return 0

    # Creating an object where pdf pages are appended to
    outpdf = PdfFileWriter()

    for ifile in pdffiles:
        if ifile[-3:] != 'pdf': continue
        inpdf = PdfFileReader(open(figdir+'/'+ifile,"rb"))
        [outpdf.addPage(inpdf.getPage(ipage)) for ipage in range(inpdf.numPages)]

    # Writing all the collected pages to a file
    outpdf.write(open(outputdir+'/'+outpdfname,"wb"))

    # delete the single files
    for ifile in pdffiles:
        os.remove(figdir+'/'+ifile)

    return outputdir+'/'+outpdfname
