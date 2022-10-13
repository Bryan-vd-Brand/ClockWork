import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-b",
        "--barcode_tsv",
        dest = "barcodes",
        nargs="+",
        required=True,
        help="TSV file containing barcodes"
    )
    
    return parser.parse_args()

#Script opens FASTA format barcode file, returns all barcode names for input filename generation
def getBarcodes():
    args = parse_args()
    barcodeFile = args.barcodes[0]
    print(F"Opening {args.barcodes} and constructing barcodeIDs")
    BarcodeIDs = []
    lines = [line.strip() for line in open(barcodeFile).readlines()]
    for line in lines:
        if '>' in line:
            barcodeID = line.replace('>','')
            BarcodeIDs.append(barcodeID)
    return BarcodeIDs  

if __name__ == "__main__":
	getBarcodes()
