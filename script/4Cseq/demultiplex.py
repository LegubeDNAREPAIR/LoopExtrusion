import itertools
import gzip
import HTSeq
import argparse

def main():
  """Script to demultiplex a fastq single end library.
  The fastq file specified in --fastq is demultiplexed 
  using the sequences provided in the barcode fasta file (--barcode)."""
  parser = argparse.ArgumentParser(description='Demultiplex a fastq file with the given barcodes provided in a fasta file.')
  parser.add_argument('--fastq', type=str, help='fastq file to be demultiplexed', required=True)
  parser.add_argument('--barcode', type=str, help='fasta file of barcodes for demultiplexing', required=True)
  parser.add_argument('--startBase', type=int, help='Base position from which the matching of the barcode is started (Default=Start of read).', required=False, default=0) 
  parser.add_argument('--endBase', type=int, help='Base position until which the barcode is used for demultiplexing (Default=minimum lenght of all barcode sequences provided).', required=False, default=None)
  parser.add_argument('--keep', type=int, help='Bases to keep at the end of the primer sequence for trimming, e.g. cutter sequence (Default=0).', required=False, default=0)
  parser.add_argument('--maxMismatches', type=int, help='Maximum number of missmatches allowed for assigning a read to a barcode (Default=0).', required=False, default=0)    
  parser.add_argument('--out', type=str, help='output directory for the demultilexed files (Default = current working directory).', required=False, default=".")
  parser.add_argument('--zip', type=str, help='define whether the output should be compressed. (Default is "True")', required=False, default="True") 

  ## get arguments
  args = parser.parse_args()
  barcodeFile = args.barcode
  startBase = args.startBase
  endBase = args.endBase
  keep = args.keep
  maxMismatch = args.maxMismatches
  outdir = args.out
  fastqFile = args.fastq
  barcodeFile = args.barcode
  outdir = args.out
  zipFiles = args.zip
  
  if outdir[-1] != "/":
     outdir = outdir + "/"

  ## read the barcodes from the fasta file into a dictionary 
  with open( barcodeFile, "r" ) as f: 
     barcodes = dict( (s.name, s.seq.rstrip().upper()) for s in HTSeq.FastaReader(f) )
     
  ## if endBase is not set the shortest length of the barcodes is used
  if endBase is None:
     endBase = min([len(barcode) for barcode in barcodes.values()])
          
  ## read the barcodes from the fasta file into a dictionary 
  with open( barcodeFile, "r" ) as f:      
     seqTable = dict( (s[startBase:endBase].seq.rstrip().upper(), s.name) for s in HTSeq.FastaReader(f) )

  ## create connections for file output.
  outfile = outdir + fastqFile.split("/")[-1].split(".")[0] + "_barcode.fastq"
  if zipFiles == "True":
    outfile=outfile+".gz"
  outFiles = {}
  if zipFiles == "False":
    for bc in barcodes:
      outFiles[bc] = open( outfile.replace("barcode", bc), "w" )
    outFileUnknown = open( outfile.replace("barcode", "unmatched"), "w" )
  else:
    for bc in barcodes:
      outFiles[bc] = gzip.GzipFile( outfile.replace("barcode", bc), "w" )
    outFileUnknown = gzip.GzipFile( outfile.replace("barcode", "unmatched"), "w" ) 

  ## if endBase is not set the shortest length of the barcodes is used
  if endBase is None:
      endBase = min([len(barcode) for barcode in barcodes.values()])  

  ## count the number of missmatches of the read against each barcode
  with  gzip.open(fastqFile, "rb") as f:
    count = 0  
    for seq in HTSeq.FastqReader( f ):
        count+=1
        if count % 100000 == 0:
            print str(count) + " reads processed"
        try:
            bc = seqTable[seq.seq[startBase:endBase]]
        except KeyError:
            aligned = False
            if maxMismatch > 0:
                match = []
                for bc in barcodes:
                    mismatch = 0
                    for b1, b2 in itertools.izip(barcodes[bc][startBase:endBase], seq.seq[startBase:endBase]):
                        if b1 != b2:
                            mismatch+=1
                        if mismatch > maxMismatch:
                            break
                    if mismatch > maxMismatch:
                        continue
                    else:
                        match.append(bc)
                if len(match) == 1:
                    seq[len(barcodes[match[0]]):].write_to_fastq_file( outFiles[ match[0] ] )
                    aligned = True
            if not aligned:
                seq.write_to_fastq_file( outFileUnknown )        
        else:
            seq[(len(barcodes[bc])-keep):].write_to_fastq_file( outFiles[ bc ] )

  for bc in barcodes:
     outFiles[bc].close()
  outFileUnknown.close()

if __name__ == "__main__":
   main()


