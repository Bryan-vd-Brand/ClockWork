import argparse
import os.path
import pandas as pd
import glob
from PyPDF2 import PdfWriter, PdfReader, Transformation
from PIL import Image
import fpdf

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-dir",
        "--directory",
        dest = "directory",
        nargs="+",
        required=True,
        help="directory containing bam files"
    )

    requiredArgs.add_argument(
        "-outDir",
        "--outDir",
        dest = "outDir",
        nargs="+",
        required=True,
        help="place to store screenshots"
    )

    requiredArgs.add_argument(
        "-KO",
        "--KnockoutReport",
        dest = "KO",
        nargs="+",
        required=True,
        help="data file containing KO calls"
    )

    return parser.parse_args()

#method composes PDFs for each well containing the different graphs/report data
def generatePDFs():

    args = parse_args()
    inputDir = args.directory[0]
    outDir = os.path.abspath(args.outDir[0])
    KOFile = args.KO[0]

    if not os.path.exists(inputDir):
        print("ERROR: directory supplied does not exist.")

    if not os.path.exists(KOFile):
        print("ERROR: file supplied does not exist.")
    

    KOTable = pd.read_table(KOFile, sep = '\t', header=0)
    Samples = KOTable['Sample'].unique()
    
    for sampleGroup in Samples:
        
        #Add well plot for sampleGroup's supersummary
        superSummary = PdfWriter()
        wellplot = PdfReader(f"{sampleGroup}_wellPlot.pdf")
        superSummary.append(wellplot)

        pageCounter = 1

        for directory in glob.iglob(inputDir + f"CRISPResso_on_*_{sampleGroup}_*"):
            
            #Dirty Hack for bad iglobs
            if ".html" in directory:
                continue
            
            lastDir = os.path.basename(os.path.normpath(directory))
            sampleName = lastDir.split('_')[2] + "_" + lastDir.split('_')[3]
            #determine if the report is for a fully frameshift well, if so add it to a super_summary
            alleles = KOTable.loc[(KOTable['SampleName'] == sampleName)]
            KOAlleles = alleles.loc[alleles['Reason'] == "Frameshift"]

            print(f"Processing PDF's for {sampleName}")
            #grab pdfs and merge: 1a.Read_barplot.pdf 1b.Alignment_pie_chart.pdf 2a.*.pdf (flip on same page?) 2b.*.pdf 5.*.Frameshift_in-frame_mutations_pie_chart.pdf 9.*.pdf
            #grab IGV screenshot and add
            merger = PdfWriter()
            merger.add_blank_page(width = 2480, height = 3508) #A4 page with 300 dpi
            basePage = merger.pages[0]
            oriMediabox = basePage.mediabox
            
            barplotPath = os.path.join(directory,"1a.Read_barplot.pdf")
            barplot = PdfReader(barplotPath).pages[0]
            barplot.mediabox.lower_left = oriMediabox.lower_left
            barplot.mediabox.upper_right = oriMediabox.upper_right
            barplot.add_transformation(Transformation().translate(tx = 25, ty = 25))
            basePage.merge_page(barplot, expand = True)

            piechartPath = os.path.join(directory,"1b.Alignment_pie_chart.pdf")
            piechart = PdfReader(piechartPath).pages[0]
            piechart.mediabox.lower_left = oriMediabox.lower_left
            piechart.mediabox.upper_right = oriMediabox.upper_right
            piechart.add_transformation(Transformation().translate(tx=1000, ty = 25))
            basePage.merge_page(piechart, expand = True)

            i = 0
            l = 0
            n = 0
            for smallNQPath in glob.iglob(os.path.join(directory,"2b.*.pdf")):
                smallNQ = PdfReader(smallNQPath).pages[0]
                smallNQ.mediabox.lower_left = oriMediabox.lower_left
                smallNQ.mediabox.upper_right = oriMediabox.upper_right
                smallNQ.add_transformation(Transformation().translate(tx= 25, ty = 950 + 350 * i))
                basePage.merge_page(smallNQ)
                i += 1
        
            for framePiePath in glob.iglob(os.path.join(directory, "5.*.pdf")):
                framePie = PdfReader(framePiePath).pages[0]
                framePie.mediabox.lower_left = oriMediabox.lower_left
                framePie.mediabox.upper_right = oriMediabox.upper_right
                framePie.add_transformation(Transformation().translate(tx= 1250, ty = 1500 + 1000 * l))
                basePage.merge_page(framePie)
                l += 1

            for freqTablePath in glob.iglob(os.path.join(directory, "9.*.pdf")):
                freqTable = PdfReader(freqTablePath).pages[0]
                freqTable.mediabox.lower_left = oriMediabox.lower_left
                freqTable.mediabox.upper_right = oriMediabox.upper_right
                freqTable.add_transformation(Transformation().translate(tx= 25 , ty = 1650 + 900 * n))
                basePage.merge_page(freqTable, expand = True)
                n += 1
                    
            #save pdf a4 to output
            writer = PdfWriter()
            writer.add_page(basePage)

            #add extra a4 with IGV screenshot and save as pdf
            sampleGroup = sampleName.split('_')[1]
            path = f'results/4_quantifyMutation/{sampleGroup}/{sampleName}'
            Image.open(f'{path}.png').convert('RGB').save(f'{path}.pdf')
            IGVPage = PdfReader(f'{path}.pdf').pages[0]
            writer.add_page(IGVPage)
            
            #Want to merge FS into supersummary
            if not alleles.empty and alleles.size == KOAlleles.size:
                
                # Create new page and merge it
                namePDF = fpdf.FPDF(format='letter')
                namePDF.add_page()
                namePDF.set_font("Arial", size = 18)
                namePDF.cell(200, 10, txt = f"{sampleGroup} \n {sampleName}", ln = 1, align="L")
                namePDF.output(f"{sampleName}_name.pdf")

                name = PdfReader(f"{sampleName}_name.pdf").pages[0]
                name.mediabox.lower_left = oriMediabox.lower_left
                name.mediabox.upper_right = oriMediabox.upper_right
                name.add_transformation(Transformation().translate(tx= 25 , ty = 2500))
                basePage.merge_page(name, expand = True)
                #Save page to large summary pdf with name indicator
                superSummary.add_page(basePage)
                superSummary.add_page(IGVPage)


            #save
            with open(os.path.join(outDir,f"{sampleName}_Summary.pdf"), "wb") as report:
                writer.write(report)

        #save super
        with open(os.path.join(outDir,f"All_Frameshifts_{sampleGroup}_Summary.pdf"), "wb") as fsreport:
                superSummary.write(fsreport)


   
if __name__ == "__main__":
	generatePDFs()